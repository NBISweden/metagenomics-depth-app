library(shiny)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(shinythemes)
library(plotly)

# r function to calculate total basepairs
total_basepairs <- function(read_length, num_reads, paired=FALSE) {
  if (paired) {
    read_length <- read_length * 2
  }
  return(read_length * num_reads)
}


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
        uiOutput("source_ui"),
        radioButtons("target",
            label="Is there a target organism?",
            choices=c("Yes","No"),
            inline=T),
        uiOutput("target_ui"), 
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
            label="Average Read length",value=150,min=1,max=25000),
        fluidRow(
          column(6,style=list("padding-right: 5px;"),
                 actionButton("click", "Update")
          ),
          column(6,style=list("padding-left: 5px;"),
                 downloadButton('download', 'Download')
          ),
        )
    ),
    
    mainPanel(
        plotOutput("rarePlot"),
        plotOutput("compPlot")
    )
  )
)
# Define server logic required to draw a histogram
server=function(input,output,session) {

  output$source_ui <- renderUI({
    if(input$source=="hos") {
      sliderInput("contam", 
                  "Contamination %", 
                  value=50, min=0, max=100)
    } 
  })
  output$target_ui <- renderUI({
    if(input$target=="Yes") {
      numericInput("gen_size", 
                   "Target Genome Size in MB", 
                   value=3, min=1, max=1000)
    } 
  })
  output$read_ui <- renderUI({
    if(input$read=="short") {
      radioButtons("read_type", 
                   "Read Type", 
                   choices=c("Single-end" = "SE", "Paired-end" = "PE"), 
                   inline=T)
    }
  })

  in_values <- reactiveValues(contam = 0, gen_size = 3, read_type = "SE")

  observeEvent(input$contam, {
    in_values$contam <- input$contam
  })
  observeEvent(input$gen_size, {
    in_values$gen_size <- input$gen_size
  })
  observeEvent(input$read_type, {
    in_values$read_type <- input$read_type
  })

  observeEvent(input$click, {

  })
}
# Run the application
shinyApp(ui = ui, server = server)