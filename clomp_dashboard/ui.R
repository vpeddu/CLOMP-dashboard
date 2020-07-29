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

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # # Application title
    titlePanel("ClompViz"),
    # 
    # # Sidebar with a slider input for number of bins
     sidebarLayout(
         sidebarPanel(
           numericInput("RPM_threshold", "Minmum RPMr threshold for filtering", 10, min = 1),
           checkboxGroupInput('Heatmap_rank', 'Rank', c('D','P','C','O','F','G','S','-'), selected = c('G','S','-'), inline = TRUE),
           uiOutput('selected_samples'),
         ),

        # Show a plot of the generated distribution
        mainPanel(
            #plotOutput("distPlot"),
            plotlyOutput(outputId = "heatmap")
        )
    )
)
)
