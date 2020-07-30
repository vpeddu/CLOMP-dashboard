#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

list.of.packages <- c("shiny","c3","formattable","DT","reactable","shinyHeatmaply","Cairo","collapsibleTree","tidyverse",
                      "viridis","visNetwork","sunburstR","plotly","dplyr")
lapply(list.of.packages,library,character.only = TRUE)

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
    
    
    # MOVE HEATMAP XAXIS TO TOP!
    
    #selected_colnames<-(c(sapply(strsplit(as.character(input$selected_samples), "\t"), "[[", 1)))
    output$heatmap<-renderPlotly({
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

        #print(colnames(heatmap_df))
      
        #heatmap_df<-heatmap_df[,c(input$selected_samples)]
        heatmap_mat<-as.matrix(heatmap_df)
        heatmap_mat<-log10(heatmap_mat)
        
        heatmap_min<-min(heatmap_mat)
        heatmap_max<-max(heatmap_mat)
        
        heatmaply(heatmap_mat,
                  dendrogram = "none",
                  heatmap_layers = theme(plot.background = element_rect(fill="#272B30"),
                                         panel.background = element_rect(fill="#272B30"),
                                         legend.background = element_rect(fill="#272B30"),
                                         text = element_text(colour = "#C8C8C8", family="Arial", size=10),
                                         #plot.title = element_text(colour = "#C8C8C8", face="bold", size=12),
                                         axis.text = element_text(colour = "#C8C8C8"),
                                         axis.line=element_blank()),
                  xlab = "", ylab = "",
                  #main = "Heatmap of RPMr values",
                  scale = "row",
                  #margins = c(60,100,40,20),
                  grid_color = "#272B30",
                  grid_width = 0.00001,
                  titleX = FALSE,
                  hide_colorbar = FALSE,
                  branches_lwd = 0.1,
                  label_names = c("Taxa", "Sample", "Log(RPM)"),
                  fontsize_row = 10,
                  fontsize_col = 10
                  #limits = c(heatmap_min,heatmap_max),
                  #subplot_heights=c(0.05, 0.95),
                  #margins = (c(10,10,10,10)),
                  
  
        ) %>% layout(autosize = TRUE)          
        #) %>% layout(width=1500, height = 500)
    })
    
    
    table_theme <- reactableTheme(
      backgroundColor <- "#272B30",
      borderColor <- "hsl(213, 10%, 17%)"
      #color <- "hsl(0, 100%, 100%)"
    )
    
    
    spotify_theme <- function() {
      search_icon <- function(fill = "none") {
        # Icon from https://boxicons.com
        svg <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24"><path fill="%s" d="M10 18c1.85 0 3.54-.64 4.9-1.69l4.4 4.4 1.4-1.42-4.39-4.4A8 8 0 102 10a8 8 0 008 8.01zm0-14a6 6 0 11-.01 12.01A6 6 0 0110 4z"/></svg>', fill)
        sprintf("url('data:image/svg+xml;base64,%s')", jsonlite::base64_enc(svg))
      }
      
      text_color <- "hsl(0, 0%, 95%)"
      text_color_light <- "hsl(0, 0%, 70%)"
      text_color_lighter <- "hsl(0, 0%, 55%)"
      #bg_color <- "hsl(0, 0%, 10%)"
      bg_color <- "#1C1E22"
      
      reactableTheme(
        color = "#272B30",
        backgroundColor = bg_color,
        borderColor = "hsl(0, 0%, 16%)",
        borderWidth = "1px",
        highlightColor = "rgba(255, 255, 255, 0.1)",
        cellPadding = "10px 8px",
        style = list(
          fontFamily = "Work Sans, Helvetica Neue, Helvetica, Arial, sans-serif",
          fontSize = "14px",
          "a" = list(
            color = text_color_light,
            "&:hover, &:focus" = list(
              textDecoration = "none",
              borderBottom = "1px solid currentColor"
            )
          ),
          ".number" = list(
            color = text_color_light,
            fontFamily = "Source Code Pro, Consolas, Monaco, monospace"
          ),
          ".tag" = list(
            padding = "2px 4px",
            color = "hsl(0, 0%, 40%)",
            fontSize = "12px",
            border = "1px solid hsl(0, 0%, 24%)",
            borderRadius = "2px"
          )
        ),
        headerStyle = list(
          color = text_color_light,
          fontWeight = 400,
          fontSize = "12px",
          letterSpacing = "1px",
          textTransform = "uppercase",
          "&:hover, &:focus" = list(color = text_color, backgroundColor = "hsl(0, 0%, 24%)")
        ),
        rowHighlightStyle = list(
          ".tag" = list(color = text_color, borderColor = text_color_lighter)
        ),
        # Full-width search bar with search icon
        searchInputStyle = list(
          paddingLeft = "30px",
          paddingTop = "8px",
          paddingBottom = "8px",
          width = "100%",
          border = "none",
          backgroundColor = bg_color,
          backgroundImage = search_icon(text_color_light),
          backgroundSize = "16px",
          backgroundPosition = "left 8px center",
          backgroundRepeat = "no-repeat",
          "&:focus" = list(backgroundColor = "rgba(255, 255, 255, 0.1)", border = "none"),
          "&:hover, &:focus" = list(backgroundImage = search_icon(text_color)),
          "::placeholder" = list(color = text_color_lighter),
          "&:hover::placeholder, &:focus::placeholder" = list(color = text_color)
        ),
        paginationStyle = list(color = text_color_light),
        pageButtonHoverStyle = list(backgroundColor = "hsl(0, 0%, 20%)"),
        pageButtonActiveStyle = list(backgroundColor = "hsl(0, 0%, 24%)")
      )
    }
    
    
    comparison_df<-comparison_df[comparison_df$name != 'root',]
    output$table<-renderReactable({
      #comparison_df$taxa 
      #comparison_df<-filter_all(any_vars(abs(.) > 0.4))
      # 
      # comparison_df<-comparison_df %>% 
      #   select_if(is.numeric) %>% 
      #   cor() %>% 
      #   as.data.frame() %>%
      #   select_if(funs(any(abs(.) > 0.4)))
      to_keep<-function(df, phylogeny, threshold){ 
        
        keep_rows<-which(df$rank %in% as.character(input$phyloRank))
        for(i in 1:nrow(df)){
          if( any(df[i,5:ncol(df)] > threshold)){
            show_list<-append(i, keep_rows)
          }
        return(keep_rows)
        }
      }
      
      
      if(input$normalizeWater==1) {
        #print(to_keep(comparison_df, input$phyloRank ,2 ))
        keep_values<-to_keep(comparison_df, input$phyloRank, input$Comparison_threshold)
        comparison_df<-comparison_df[keep_values,]
        #comparison_df<-comparison_df[comparison_df$rank %in% as.character(input$phyloRank),]
        graph_df<-comparison_df
        graph_df[,c(2,3,4)]<-NULL
      } else if (input$normalizeWater==2) {
        graph_df <- x()
        rownames(graph_df)<-graph_df$name
        to_remove<-which(grepl(paste(c('name','rank','taxID','lineage'), collapse = "|"),colnames(graph_df)))
        graph_df[,to_remove]<-NULL
        graph_df<-rownames_to_column(graph_df, var = 'name')
      }
      
      colorpal <- viridis(255,begin=0,end=1)

      #GnYlRd <- function(x) rgb(colorRamp(c("#63be7b", "#ffeb84", "#f8696b"))(x), maxColorValue = 255)
      GnYlRd <- function(x) rgb(colorRamp(colorpal[c(15:255)])(x), maxColorValue = 255)
      reactable(graph_df, 
                defaultColDef = colDef(
                  style = function(value) {
                    if (!is.numeric(value))
                    {return()}
                    else 
                      normalized <- (value - min(graph_df[,2:ncol(graph_df)])) / (max(graph_df[,2:ncol(graph_df)]) - min(graph_df[,2:ncol(graph_df)]))
                    color <- GnYlRd(normalized)
                    list(background = color) 
                    #list(background = 1) 
                  },
                  format = colFormat(digits = 0),
                  #minWidth = 50,
                  #headerStyle = list(background = "#f7f7f8")
                ),
                searchable = TRUE, 
                height = "auto",
                width = "auto",
                bordered = TRUE,
                outlined = FALSE,
                resizable = TRUE,
                highlight = TRUE,
                fullWidth = TRUE,
                wrap = FALSE,
                theme = spotify_theme,
                columns = list(
                  name = colDef(name = "Taxa", sortable = TRUE, align = "left", resizable = TRUE, minWidth=1600,
                                cell = function(value, index) {
                                  # Render as a link
                                  url <- sprintf("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s", comparison_df$taxID[index])
                                  htmltools::tags$a(href = url, target = "_blank", as.character(value))
                                }, 
                                style = list(fontFamily = "Helvetica", color = "#8ec8c8",whiteSpace = "pre")
                  )
                ) )
      
    })
})
