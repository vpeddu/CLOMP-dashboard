#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(c3)
library(formattable)
library(DT)
library(reactable)
#library(pavian)
library(shinyHeatmaply)
library(Cairo)
library(collapsibleTree) 
library(tidyverse)
library(viridis)
library(visNetwork)
library(sunburstR)
library(plotly)

# Read and merge pavian TSVs 
# Returns dafaframe of dataframe of RPM values 
prep_data<-function(path, taxa_class){
    path <- '/Users/uwvirongs/Documents/clomp-dashboard-data/test_data/'
    setwd(path)
    
    files<-list.files(path, pattern = '*.tsv',full.names = TRUE)
    
    reports<-pavian::read_reports(files)
    
    filenames<-strsplit(basename(files), '[.]')
    
    reports_merged<-as.data.frame(pavian::merge_reports2(reports,fix_taxnames = TRUE, col_names = filenames))
    
    keep<-which(grepl('clade',colnames(reports_merged)))
    
    only_clade<-reports_merged[,c(1:4, keep)]
    
    #colnames(only_clade)[keep] <- strsplit(colnames(only_clade), 'f')
    
    new_colnames<-c()
    for(i in 5:length(colnames(only_clade))){ 
        #print(colnames(only_clade[i]))
        first<-strsplit(colnames(only_clade)[i],'[..]')[[1]][4]
        new_colnames<-append(new_colnames,first)
    }
    colnames(only_clade)[keep]<-new_colnames
    colnames(only_clade)[1:4]<-c('name','rank','taxID','lineage')
    
    
    final_tsv<-only_clade
    
    for(i in 5:ncol(only_clade)){ 
        final_tsv[,i]<-only_clade[,i] / ((only_clade[1,i] + only_clade[2,i] ) / 1e6) 
    }
    
    
    final_tsv<-final_tsv[complete.cases(final_tsv),]
    
    to_remove<-c()
    
    for(i in 1:nrow(final_tsv)){ 
        if( all(final_tsv[i,5:ncol(final_tsv)] == 0 )){ 
            to_remove<-append(to_remove, i)
        }
    }
    if(length(to_remove > 0)){
        final_tsv_zero_removed<-final_tsv[-to_remove,]
    } else{
        final_tsv_zero_removed<-final_tsv
    }
    
    taxa_filtered<-final_tsv_zero_removed[which(grepl(paste(taxa_class,collapse="|"),final_tsv_zero_removed$rank)),]
    
    final<-clean_order(taxa_filtered)
    
    return(final)
}

# Returns dataframe of RPMr values (water excluded)
# Samples must have "DNA" or "RNA" in sample names 
# Water samples must be labeled "H2O_RNA" and "H2O_RNA" 
RPMr_create<-function(df){ 
    DNA_cols<-which(grepl('DNA',colnames(df)))
    RNA_cols<-which(grepl('RNA',colnames(df)))
    
    RNA_df<-df[,RNA_cols]
    DNA_df<-df[,DNA_cols]
    
    water_RNA_col<-which(grepl('H2O_RNA',colnames(RNA_df)))
    water_DNA_col<-which(grepl('H2O_DNA',colnames(DNA_df)))
    
    # Adding pseudocount of 1 to account for RPM of 0 in water controls
    # Avoids divide by 0 error 
    RNA_df<-RNA_df + 1 
    DNA_df<-DNA_df + 1 
    
    RNA_r_df<- RNA_df / RNA_df[,water_RNA_col]
    DNA_r_df<- DNA_df / DNA_df[,water_DNA_col]
    
    # getting rid of water columns 
    RNA_r_df[,water_RNA_col]<-NULL
    DNA_r_df[,water_DNA_col]<-NULL
    
    # Binding back into same dataframe
    final<-cbind(df[,c(1:4)], RNA_r_df, DNA_r_df)
    
    # Reordering columns
    sample_order<-order(colnames(final[5:ncol(final)])) + 4
    
    final<-final[,c(1:4,sample_order)]
    }


# Clean column names and reorder columns 
clean_order<-function(df){ 
    for(i in 5:ncol(df)){ 
        #print(colnames(df)[i])
        #colnames(df)[i]<-
        colnames(df)[i]<-(paste0(strsplit(colnames(df)[i], '_')[[1]][3], "_" ,strsplit(colnames(df)[i], '_')[[1]][4]))
    }
    samples<-df[5:ncol(df)]
    samples<-samples[,order(colnames(samples))]
    
    sample_order<-order(colnames(df[5:ncol(df)])) + 4
    
    df<-df[,c(1:4,sample_order)]
    
    
    return(df)
}


filepath  = '/Users/uwvirongs/Documents/clomp-dashboard-data/test_data/'

df<-(prep_data(filepath, 'G'))
colnames(df)<-gsub("\t", "", colnames(df))
df$name<-trimws(df$name, which = c("left"), whitespace = "[ \t\n]")

# global data frame of RPM values with name, lineage, rank, and sample names
comparison_df<-(prep_data(filepath, "*"))

# global data frame of RPMr values (RPM_sample / RPM_water)
RPM_r_df<-RPMr_create(comparison_df)



# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    # output$distPlot <- renderPlot({
    # 
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white')
    # 
    # })

    show_list<-c()
    
    uniques<-unique(c(sapply(strsplit(as.character(colnames(RPM_r_df)[5:ncol(RPM_r_df)]), "_"), "[[", 1)))
    print(uniques)
    output$selected_samples <- renderUI({
      #checkboxGroupInput("selected_samples", "select samples:", choices = colnames(RPM_r_df[5:ncol(RPM_r_df)]), selected = colnames(RPM_r_df[5:ncol(RPM_r_df)]))
      checkboxGroupInput("selected_samples", "select samples:", choices = uniques, selected = uniques)
    })
  
    
    x <- reactive({
        z<-RPM_r_df[RPM_r_df$rank == input$Heatmap_rank,]
        #print(input$Heatmap_rank)

        for(i in 1:nrow(z)){
            if( any(z[i,5:ncol(z)] > input$RPM_threshold)){
                show_list<-append(i, show_list)
            }
            else{ 
                #print('failed')
                }
        }
        show_list<-unique(show_list)
        
        selected_colnames<-(c(sapply(strsplit(as.character(input$selected_samples), "\t"), "[[", 1)))
        to_include<-which(grepl(paste(selected_colnames, collapse = "|"),colnames(RPM_r_df)))
        to_include<-append(unlist(to_include),c(1:4) )
        print(to_include)
        final<-z[show_list,to_include]
        final

    })
    
    
    
   
    # 
    # filter_df<- reactive({
    #   print('here')
    #   selected_colnames<-(c(sapply(strsplit(as.character(input$selected_samples), "\t"), "[[", 1)))
    #   to_include<-which(colnames(x()) %in% selected_colnames)
    #   print(to_include)
    #   out<-x()[,to_include]
    #   out
    # })
    # 
    
    #selected_colnames<-(c(sapply(strsplit(as.character(input$selected_samples), "\t"), "[[", 1)))
    output$heatmap<-renderPlotly({
      print('here')
      #print(length(filter_df()))
      
    #to_include<-c(5:ncol(x()))
    #print(class(input$selected_samples))
    #print(strsplit(input$selected_samples, '\t') )
    #print(asdfasdf[1])
    
      #heatmap_df<-x()
      #to_include<-c('5')
      #to_include<-append(1,to_include)
        #x<-x()[,to_include]
        #print(colnames(x))
        #original<-x()
        heatmap_df<-x()
        #print(colnames(heatmap_df))
      
      
        #heatmap_df<-heatmap_df[,c(input$selected_samples)]
        
        
        #print(input$selected_samples)
        #heatmap_names<-colnames(heatmap_df)

        
        
        #print(dim(heatmap_df))
        rownames(heatmap_df)<-heatmap_df$name
        #heatmap_df$name<-NULL
        to_remove<-which(grepl(paste(c('name','rank','taxID','lineage'), collapse = "|"),colnames(heatmap_df)))
        heatmap_df[,to_remove]<-NULL
        #heatmap_df[,c('name','rank','taxID','lineage')]<-NULL

        print(colnames(heatmap_df))
        
        
        #heatmap_df<-heatmap_df[,c(input$selected_samples)]
        heatmap_mat<-as.matrix(heatmap_df)
        heatmap_mat<-log10(heatmap_mat)
        
        heatmap_min<-min(heatmap_mat)
        heatmap_max<-max(heatmap_mat)
        
        
        heatmaply(heatmap_mat,
                  dendrogram = "none",
                  xlab = "", ylab = "",
                  main = "Heatmap of RPMr values",
                  scale = "row",
                  #margins = c(60,100,40,20),
                  grid_color = "black",
                  grid_width = 0.00001,
                  titleX = FALSE,
                  hide_colorbar = FALSE,
                  branches_lwd = 0.1,
                  label_names = c("Taxa", "Sample", "Log(RPM)"),
                  fontsize_row = 5, 
                  #limits = c(heatmap_min,heatmap_max),
                  #fontsize_col = 5,
                  #subplot_heights=c(0.05, 0.95),
                  #margins = (c(10,10,10,10)),
                  
                  heatmap_layers = theme(axis.line=element_blank())
        ) %>% layout(width=1500, height = 500)
    })
    
    
})
