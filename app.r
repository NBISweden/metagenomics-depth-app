library(shiny)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(shinythemes)
library(plotly)


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Metagenomics Coverage Calculator"),
    theme = shinytheme("flatly"),
  sidebarLayout(
    sidebarPanel(
        radioButtons("source",
            label="Metagenome Source",
            choices=c("Human/Host" = "hos",
                    "Environmental" = "env"),
            inline=T),
        uiOutput("source_ui"), # to include contamination %
        radioButtons("target",
            label="Is there a target organism?",
            choices=c("Yes","No"),
            inline=T),
        uiOutput("target_ui"), # to include target species genome size
        sliderInput("coverage",
            "X coverage expected for target species",
            min = 0,
            max = 100,
            value = 10),
        radioButtons("read",
            label="Sequencing machine type",
            choices=c("Illumina/Short-read" = "short",
                    "PacBio/ONT/long-read" = "long"),
            inline=T),
        uiOutput("read_ui"), # to include SE and PE options
        numericInput("length",
            label="Read length",value=150,min=1,max=25000)
    ),
    
    mainPanel(
        plotOutput("rarePlot"),
        plotOutput("compPlot")
    )
  )
)
# Define server logic required to draw a histogram
server=function(input,output,session) {
  observe({
    if(something) {
      updateSelectInput(session,"select-input",label="selectInput",choices=c("D","E","F"))
      updateNumericInput(session,"numeric-input",label="numericInput",value=10,min=1,max=10)
      updateSliderInput(session,"slider-input",label="sliderInput",value=8,min=1,max=10)
    }
  })
}
# Run the application
shinyApp(ui = ui, server = server)