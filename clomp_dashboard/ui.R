#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(shinythemes)
library(reactable)

ranks <-
  c(
    #"root",
    #"no_rank",
    "superkingdom",
    #"clade",
    "kingdom",
    "phylum",
    #"subphylum",
    #"superclass",
    "class",
    #"superorder",
    "order",
    #"suborder",
    #"infraorder",
    #"parvorder",
    #"superfamily",
    "family",
    #"subfamily",
    "genus",
    "species"
    #"subgenus",
    #"subspecies",
    #"infraclass",
    #"subclass",
    #"cohort",
    #"subcohort",
    #"subkingdom",
    #"strain",
    #"tribe",
    #"subtribe",
    #"species_group",
    #"species_subgroup"
  )

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("slate"),

    # # Application title
    titlePanel("ClompViz"),
    # 
    # # Sidebar with a slider input for number of bins
     sidebarLayout( fluid = TRUE,
         sidebarPanel(
           width = 3,
           numericInput("RPM_threshold", "Minimum RPMr threshold for filtering", 10, min = 1),
           checkboxGroupInput('Heatmap_rank', 'Rank', ranks, selected = ranks, inline = FALSE),
           uiOutput('selected_samples')
         ),

        # Show a plot of the generated distribution
        mainPanel(
            #plotOutput("distPlot"),
            plotlyOutput(outputId = "heatmap")
        )
        
    ),
    sidebarLayout( fluid = TRUE,
      sidebarPanel( 
        width = 3,
        numericInput("Comparison_threshold", "Minimum RPM threshold for filtering", 10, min = 1),
        #uiOutput('phyloRank'),
        #checkboxGroupInput('phyloRank', 'Rank', c('D','P','class','O','F','G','S','-'), selected = c('class','G','S','-'), inline = TRUE),
        checkboxGroupInput('phyloRank', 'Rank',ranks, selected = ranks, inline = FALSE),
        radioButtons("normalizeWater", h3(""), choices = list("Include water samples" = 1, "Normalize to water samples" = 2), selected = 1)
        ),
      mainPanel( 
        reactableOutput("table")
        )
      
    )
)
)
