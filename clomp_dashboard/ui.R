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
           checkboxGroupInput('Heatmap_rank', 'Rank', c('D','P','C','O','F','G','S','-'), selected = c('G','S','-'), inline = TRUE),
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
        checkboxGroupInput('phyloRank', 'Rank', c('D','P','C','O','F','G','S','-'), selected = c('G','S','-'), inline = TRUE),
        radioButtons("normalizeWater", h3(""), choices = list("Include water samples" = 1, "Normalize to water samples" = 2), selected = 1)
        ),
      mainPanel( 
        reactableOutput("table")
        )
      
    )
)
)
