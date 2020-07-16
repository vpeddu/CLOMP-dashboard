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
    path <- '/Users/gerbix/Documents/vikas/CLOMP/clomp_view/test_data'
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




filepath  = '/Users/gerbix/Documents/vikas/CLOMP/clomp_view/test_data'

df<-(prep_data(filepath, 'G'))
colnames(df)<-gsub("\t", "", colnames(df))
df$name<-trimws(df$name, which = c("left"), whitespace = "[ \t\n]")

# global data frame of RPM values with name, lineage, rank, and sample names
comparison_df<-(prep_data(filepath, "*"))

# global data frame of RPMr values (RPM_sample / RPM_water)
RPM_r_df<-RPMr_create(comparison_df)



# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')

    })

    show_list<-c()
    x <- reactive({
        z<-RPM_r_df[RPM_r_df$rank == input$Heatmap_rank,]
        print(input$Heatmap_rank)

        for(i in 1:nrow(z)){
            if( any(z[i,5:ncol(z)] > input$RPM_threshold)){
                show_list<-append(i, show_list)
            }
            else{ 
                #print('failed')
                }
        }
        show_list<-unique(show_list)
        final<-z[show_list,]
        final

    })


    output$heatmap<-renderPlotly({
        heatmap_df<-x()
        
     
        rownames(heatmap_df)<-heatmap_df$name
        heatmap_df[,1:4]<-NULL
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
