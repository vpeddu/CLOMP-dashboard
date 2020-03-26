#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(c3)
library(formattable)

# Define UI for application that draws a histogram

ui <- fluidPage(

    # Application title
    titlePanel("RPM values across all samples"),

    
    numericInput("RPM_threshold", "Minmum RPM threshold for filtering", 200, min = 1),
    
    
    checkboxInput('show_subchart', 'Show subchart', value = FALSE, width = NULL), 
    # Sidebar with a slider input for number of bins 
   
           c3Output("RPM_plot",width = "100%", height = 'auto'), 
    
    titlePanel("Comparison"),
    numericInput("Comparison_threshold", "Minmum RPM threshold for filtering", 10, min = 1),
    
            box(formattableOutput("table"))
    
    )


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    prep_data<-function(path, taxa_class){
        setwd(path)
        
        files<-list.files(pattern = '*.tsv')
        
        taxa_detect<-function(df, taxid){ 
            temp_rpm<-df$RPM[which(df$taxid == taxid)]
            if(identical(temp_rpm, numeric(0))){ 
                temp_rpm <- 0
            }
            return(temp_rpm)
        }
            
        
        
        for(i in 1:length(files)){
            #print(files[i])
            temp_tsv<-read.csv(files[i], sep = "\t", col.names = c('percent_clade_reads', 'number_clade_reads_rooted_at_taxon','number_clade_reads_this_taxon', 'taxa', 'taxid', 'name'), header = FALSE)
            total_reads = temp_tsv$number_clade_reads_rooted_at_taxon[2] + temp_tsv$number_clade_reads_rooted_at_taxon[1]
            temp_tsv$RPM = temp_tsv$number_clade_reads_this_taxon  / (total_reads / 1e6)
            temp_tsv$cumulative_RPM<-temp_tsv$number_clade_reads_rooted_at_taxon / (total_reads / 1e6)
            temp_tsv$taxa<-trimws(temp_tsv$taxa , which = "both", whitespace = "\t")
            temp_tsv$name<-as.character(temp_tsv$name)
            file_name = strsplit(files[i],"_L001")[[1]][1]
            
            #initialize dataframe for below
            if( i == 1 ){ 
                final_tsv<-temp_tsv[,c(4,5,6,8)] 
                colnames(final_tsv)[4]<-file_name
            }
            
            else{ 
                final_tsv[,(i+3)]<-0
                new_list<-c()
                for(j in 1:nrow(temp_tsv)){ 
                    #index is which row of the current pavian tsv equals the current taxid
                    index<-which(temp_tsv[j,5] == final_tsv[,2])
                    
                    if(identical(index, integer(0))){ 
                        final_tsv[(nrow(final_tsv)+ 2), ]<- 0
                        final_tsv[nrow(final_tsv),1]<-temp_tsv[j,4]
                        final_tsv[nrow(final_tsv),2]<-temp_tsv[j,5]
                        final_tsv[nrow(final_tsv),3]<-as.character(temp_tsv[j,6])
                        final_tsv[nrow(final_tsv),ncol(final_tsv)]<-temp_tsv[j,8]
                        
                        next
                    }
                    else{ 
                        final_tsv[index,i+3]<-temp_tsv[j,8]
                    }
                    colnames(final_tsv)[i+3]<-file_name
                }
            }
        }
        
        
        final_tsv<-final_tsv[complete.cases(final_tsv),]
        
        to_remove<-c()
        
        for(i in 1:nrow(final_tsv)){ 
            if( all(final_tsv[i,3:ncol(final_tsv)] == 0 )){ 
                to_remove<-append(to_remove, i)
            }
        }
        if(length(to_remove > 0)){
            final_tsv_zero_removed<-final_tsv[-to_remove,]
        }
        else{
            final_tsv_zero_removed<-final_tsv
        }
        
        taxa_filtered<-final_tsv_zero_removed[which(grepl(taxa_class,final_tsv_zero_removed$taxa)),]
        # 
        # sample_list<-c()
        # taxa_list<-c()
        # rpm_list<-c()
        # taxid_list<-c()
        # name_list<-c()
        # for(i in 1:nrow(final_tsv_zero_removed)){ 
        #     for(j in 4:ncol(final_tsv_zero_removed)){ 
        #         sample_list<-append(sample_list,colnames(final_tsv_zero_removed)[j])
        #         taxa_list<-append(taxa_list,final_tsv_zero_removed[i,1])
        #         rpm_list<-append(rpm_list,final_tsv_zero_removed[i,j])
        #         taxid_list<-append(taxid_list,final_tsv_zero_removed[i,2])
        #         name_list<-append(name_list, final_tsv_zero_removed[i,3])
        #     }  
        # }
        # accumulated_df<-as.data.frame(cbind(sample_list,taxa_list,rpm_list,taxid_list,name_list))
        # graph_df<-accumulated_df[which(accumulated_df$taxa_list == 'G'),]
        # 
        # graph_df$rpm_list<-as.numeric(as.character(graph_df$rpm_list))
        # graph_df$rpm_list<-round(graph_df$rpm_list, digits = 0)
        # 
        # graph_df$taxid_list<-as.character(graph_df$taxid_list)
        # graph_df<-graph_df[-which(graph_df$taxid_list == '9605'),]
        
        return(taxa_filtered)
    }
    
    filepath  = '/Users/vikas/Documents/UW/Greninger_lab/clomp_viz/scratch/remove_euk_test'
  
        df<-(prep_data(filepath, 'S'))
        colnames(df)<-gsub("\t", "", colnames(df))
        df$name<-trimws(df$name, which = c("left"), whitespace = "[ \t\n]")
        
    output$RPM_plot <- renderC3({
            
        
        to_remove<-c()
        
        for(i in 1:nrow(df)){ 
            if(all(df[i,(4:ncol(df))] < input$RPM_threshold)){ 
                to_remove<-append(to_remove, i)
            }
        }
        if(length(to_remove > 0)){
            df<-df[-to_remove,]
        }
        
        
        
        #comment out when done with test_set
        
        
        df<-df[,c(2:ncol(df))]
        data_df<-t(df[,c(3:ncol(df))])
        colnames(data_df)<-as.character(df[,2])
        #colnames(data_df)<-gsub("[[:space:]]", "", colnames(data_df))
        sample_names<-as.character(rownames(data_df))
        #taxa_names<-gsub("[[:space:]]", "", colnames(data_df))
        data_df<-log10(data_df + 1)
        
        
        
        if(input$show_subchart){ 
        data.frame(data_df) %>%
            c3() %>% 
            c3_bar(stacked = FALSE, rotated = FALSE, zerobased = TRUE)%>% 
            tooltip( grouped = FALSE) %>%
            xAxis( type = 'category', categories = sample_names, rotate = 45)  %>%
            legend(position = 'right') %>%
            zoom(type = 'scroll') %>% 
            subchart(height = 60) 
        }
        else {        data.frame(data_df) %>%
                c3() %>% 
                c3_bar(stacked = FALSE, rotated = FALSE, zerobased = TRUE)%>% 
                tooltip( grouped = FALSE) %>%
                xAxis( type = 'category', categories = sample_names, rotate = 45)  %>%
                yAxis(label = 'log10(RPM)') %>%
                legend(position = 'right') %>%
                zoom(type = 'scroll')
            }
    })
    
    
    output$table <- renderFormattable({
        
        comparison_df<-(prep_data(filepath, '*'))
        #comparison_df$taxa 
        to_remove<-c()
        
        for(i in 1:nrow(comparison_df)){ 
            if(all(comparison_df[i,(4:ncol(comparison_df))] < input$Comparison_threshold)){ 
                to_remove<-append(to_remove, i)
            }
        }
        if(length(to_remove > 0)){
            comparison_df<-comparison_df[-to_remove,]
        }
        
        #works but can't sort. Probably need to add buttons to patch it together (formattable is static and datatable is uggo)
        formattable(as.data.frame(comparison_df))
        
        })
}

# Run the application 
shinyApp(ui = ui, server = server)

