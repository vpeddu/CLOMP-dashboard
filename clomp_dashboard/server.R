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


#filepath  = '/Users/gerbix/Documents/vikas/CLOMP/clomp_view/test_data'
#filepath = '/Users/gerbix/Documents/vikas/CLOMP/clomp_view/test_data_new'
filepath = '/Users/uwvirongs/Downloads/CLOMP-dashboard/test_data_new'

# Read and merge pavian TSVs 
# Returns dafaframe of dataframe of RPM values 
prep_data<-function(path, taxa_class){
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
    
    #print(only_clade[,c(1,5:6)])
    final_tsv<-only_clade
    #print('colname')
    #print(colnames(final_tsv))
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
    

    for(k in 1:length(files)){
      if( k == 1 ){
      temp<-read.csv(file = files[k], sep = "\t", header = FALSE)[ ,4:5]
      }
      else {
        temp2 = read.csv(file = files[k], sep = "\t", header = FALSE)[ ,4:5]
        new_taxids = temp2[,2][temp2[,2] %in% temp[,2] == FALSE]
        temp2<-temp2[new_taxids,]
        temp<-rbind(temp,temp2)
        }
    }
    #print(temp)
    #print(unique(temp[,2]))
    taxa_filtered<-final_tsv_zero_removed[which(grepl(paste(taxa_class,collapse="|"),final_tsv_zero_removed$rank)),]
    #print(temp)


    
    final<-clean_order(taxa_filtered)
    
    for(f in 1:nrow(final)){
      #print(f)
    #print(temp[,1][which(temp[,2] == final$taxID[f])[1]])  
      final$rank[f]<-as.character(temp[,1][which(temp[,2] == final$taxID[f])[1]])
      #print(final$rank[f])
      #print(final$taxID[f])
    }
    #print('here')
    
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
    #print(uniques)
    output$selected_samples <- renderUI({
      #checkboxGroupInput("selected_samples", "select samples:", choices = colnames(RPM_r_df[5:ncol(RPM_r_df)]), selected = colnames(RPM_r_df[5:ncol(RPM_r_df)]))
      checkboxGroupInput("selected_samples", "select samples:", choices = uniques, selected = uniques)
    })
  
    
    x <- reactive({
        z<-RPM_r_df[RPM_r_df$rank %in% input$Heatmap_rank,]
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
        #print(to_include)
        final<-z[show_list,to_include]
        final

    })
   


    
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
      
      text_color <- "hsl(0, 0%, 100%)"
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

      to_keep<-function(df, phylogeny, threshold,ranges){ 
        #print(colnames(df))
        #print(unique(df$rank))
        keep_rows<-which(df$rank %in% as.character(phylogeny))
        threshold_filtered<-c()
        
        for(i in 1:length(keep_rows)){
          
          if( any(df[keep_rows[i],ranges] > threshold)){
            
            threshold_filtered<-append(keep_rows[i], threshold_filtered)
          }
        }
        out_list<-unique((threshold_filtered))
        #print(out_list)
        return(out_list)
      }
      
      
      if(input$normalizeWater==1) {
        #print(to_keep(comparison_df, input$phyloRank ,2 ))
        keep_values<-to_keep(comparison_df, input$phyloRank, input$Comparison_threshold,5:ncol(comparison_df))
        #print(keep_values)
        comparison_df<-comparison_df[keep_values,]
        #print(keep_values)
        #comparison_df<-comparison_df[comparison_df$rank %in% as.character(input$phyloRank),]
        graph_df<-comparison_df
        graph_df[,c(2,3,4)]<-NULL
      } else if (input$normalizeWater==2) {
        graph_df <- x()
        rownames(graph_df)<-graph_df$name
        keep_values<-to_keep(graph_df, input$phyloRank, input$Comparison_threshold, 1:(ncol(graph_df) - 5))
        #print(keep_values)
        graph_df<-graph_df[keep_values,]
        to_remove<-which(grepl(paste(c('name','rank','taxID','lineage'), collapse = "|"),colnames(graph_df)))
        graph_df[,to_remove]<-NULL
        graph_df<-rownames_to_column(graph_df, var = 'name')
      }
      
      colorpal <- viridis(255,begin=0,end=1,option="D")

      
      #print(colnames(graph_df))
      GnYlRd <- function(x) rgb(colorRamp(colorpal[c(15:255)])(x), maxColorValue = 255)
      #GnYlRd <- function(x) rgb(colorRamp(c("#63be7b", "#ffeb84", "#f8696b"))(x), maxColorValue = 255)
      reactable(graph_df, 
                defaultColDef = colDef(
                  style = function(value) {
                    if (!is.numeric(value))
                    {return()}
                    else {

                      if (value <= 1) {
                        normalized = 0
                      } else {
                      normalized <- log10(value)/(log10(max(graph_df[,2:ncol(graph_df)])))
                      
                      }
                      
                      #print("that was value")
                      #print(normalized)
                    }
                    #print(min(graph_df[,2:ncol(graph_df)]))
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
                  name = colDef(name = "Taxa", sortable = TRUE, align = "left", resizable = TRUE, minWidth=250,
                                cell = function(value, index) {
                                  # Render as a link
                                  url <- sprintf("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s", comparison_df$taxID[index])
                                  htmltools::tags$a(href = url, target = "_blank", as.character(value))
                                }, 
                                style = list(fontFamily = "Helvetica", color = "#8ec8c8",whiteSpace = "pre")
                  )
                ) )
      
    })
    
    
    # 
    # output$sunburst<-renderSunburst({
    #   #col_index = which(colnames(comparison_df) == input$select)
    #   col_index = 1
    #   test<-comparison_df[,c(4,col_index)]
    #   rough<-sapply(test$lineage, make_relation)
    #   
    #   all_values<-do.call(rbind, rough)
    #   all_values<-rownames_to_column(all_values, "lineage")
    #   all_values<-distinct(all_values,from, to, .keep_all = TRUE)
    #   
    #   all_values$lineage<-sapply(strsplit(as.character(all_values$lineage), "\\."), "[[", 1)
    #   
    #   merged<-left_join(all_values, test ,by = 'lineage' )
    #   colnames(merged)[length(colnames(merged))]<-'size'
    #   
    #   merged$taxa<-sapply(strsplit(as.character(merged$from), "_"), "[[", 1)
    #   merged$from<-sapply(strsplit(as.character(merged$from), "_"), "[[", 2)
    #   merged$to<-sapply(strsplit(as.character(merged$to), "_"), "[[", 2)
    #   map_colors<-colorRampPalette(viridis(12))
    #   colors = map_colors(length(unique(merged$taxa)))
    #   merged$color= NA
    #   for(i in 1:length(unique(merged$taxa))){ 
    #     indexes <- which(merged$taxa == unique(merged$taxa)[i])
    #     merged$color[indexes]<-colors[i]
    #     }
    #   
    #   merged$lineage<-gsub("-_root", "root", merged$lineage)
    #   #merged$lineage<-gsub('-_','r_', merged$lineage)
    #   #merged$lineage<-gsub('-_','', merged$lineage)
    #   #merged$lineage<-gsub('\\|','-', merged$lineage)
    #   
    #   sunburst_df<-data.frame(merged$lineage, merged$size)
    #   sunburst_df<-sunburst_df[order(sunburst_df$merged.lineage),]
    #   
    #   sund2b(sunburst_df,
    #          #count = TRUE,
    #          #percent = TRUE,
    #          #height = "2000px",
    #          #width = "100%"
    #          # explanation would put the % and count out
    #          #,explanation = "function(d){return d.data.name}"
    #   )
    #   
    #   #merged<-merged %>% mutate(trimmed = str_replace_all( merged$lineage, '\\|', ""))  
    #   #merged<-merged %>% mutate(trimmed_final = str_replace_all( merged$trimmed, '_', ""))  
    #   
    # #})
    #   })
    # 
    # 
    

    make_relation<-function(sample){ 
      #sample<-data.frame(sample, stringsAsFactors = FALSE)
      temp<-unlist(strsplit(sample, '|', fixed = TRUE))
      #print(temp)
      count = 1
      to<-c()
      from<-c()
      if(length(temp) > 10)
        end = 10
      else
        end = length(temp)
      while(count < end ){ 
        from<-append(from, temp[count])
        to<-append(to, temp[count + 1])
        count = count + 1
      }
      df<-data.frame(from = from, to = to)
      return(df)
    }
    
    
    output$visnetworkPlot<-renderVisNetwork({
      col_index = which(colnames(comparison_df) == input$visnetworkSelect)
      test<-comparison_df[,c(4,col_index)]
      rough<-sapply(test$lineage, make_relation)
      
      all_values<-do.call(rbind, rough)
      all_values<-rownames_to_column(all_values, "lineage")
      all_values<-distinct(all_values,from, to, .keep_all = TRUE)
      
      all_values$lineage<-sapply(strsplit(as.character(all_values$lineage), "\\."), "[[", 1)
      
      merged<-left_join(all_values, test ,by = 'lineage' )
      colnames(merged)[length(colnames(merged))]<-'size'
      
      merged$taxa<-sapply(strsplit(as.character(merged$from), "_"), "[[", 1)
      merged$from<-sapply(strsplit(as.character(merged$from), "_"), "[[", 2)
      merged$to<-sapply(strsplit(as.character(merged$to), "_"), "[[", 2)
      
      merged<-merged[merged$size > 1,]
      
      
      colorpal2 <- viridis(200,begin=0,end=1,option="B")
      colorpal2 <- colorpal2[20:200]
      map_colors<-colorRampPalette(colorpal2)
      
      colors = map_colors(length(unique(merged$taxa)))
      merged$color= NA
      for(i in 1:length(unique(merged$taxa))){ 
        indexes <- which(merged$taxa == unique(merged$taxa)[i])
        merged$color[indexes]<-colors[i]
        
      }
      
      nodes<-data.frame(id = (merged$to), label = (merged$to), size = merged$size, title = merged$size, color = merged$color)
      nodes$level <- str_count(as.character(merged$lineage), '\\|') 
      #print(nodes)
      #print(merged$lineage)
      
      nodes$title <- paste(nodes$id," :RPM = ", round(nodes$title))
      
      edges<-data.frame(from = merged$from, to = merged$to, size = merged$size)
      
      nodes$id<-as.character(nodes$id)
      
      
      uniques<-unique(nodes$id  )
      keep<-c()
      for(i in 1:length(uniques)){ 
        indexes<-which(nodes$id == uniques[i])
        keep<-append(keep,(indexes[(which.max(nodes$size[indexes]))]))
        #print(which(nodes$size == max(nodes$size[nodes$id == i])))
        #keep<-append(keep, which(nodes$size == max(nodes$size[nodes$id == i])))
      }
      #print(nodes)
      
      nodes<-nodes[keep,]
      edges<-edges[complete.cases(edges),]
      #edges<-which(edges$from[edges$from == edges$to])
      #edges<-edges[-which(edges$from == 'Actinobacteria'),]
      #print(edges)
      nodes$size[nodes$title <= 1] <- NA
      #nodes$size[nodes$title >= 1] <- nodes$size[nodes$title >= 1] + 3 
      nodes$size<-log10(nodes$size) ^2
      
      nodes<-nodes[complete.cases(nodes$title),]
      edges$size <- log10(edges$size ^ 2 )
      #print(edges$size)
      edges$value <- edges$size
      #print(nodes[nodes$size == 0])
      
      
      
      visNetwork(nodes = nodes, edges = edges, width = "100%") %>%
        #visHierarchicalLayout(levelSeparation = 500) %>%
        visHierarchicalLayout(sortMethod = 'directed', direction="UD") %>%
        visNodes(label = nodes, font = list(color = 'white')) %>% 
        #visEdges(value= edges$size , scaling = list(edges$size)) %>%
        #visNodes(width)
        visPhysics(stabilization = FALSE) %>%
        #  visEdges(smooth = FALSE) %>%
        visPhysics(solver = "barnesHut")
      #            hierarchicalRepulsion = list(gravitationalConstant = -200))
      
    })
})

