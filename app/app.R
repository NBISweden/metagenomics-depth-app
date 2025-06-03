library(shiny)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(shinythemes)
library(plotly)
library(ggpubr)

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
    theme = shinytheme("darkly"),
  sidebarLayout(
    sidebarPanel(
        radioButtons("source",
            label="Metagenome Source",
            choices=c("Human/Host" = "host",
                    "Environmental" = "env"),
            inline=T),
        uiOutput("source_ui"),
        radioButtons("target",
            label="Is there a target organism?",
            choices=c("Yes","No"),
            inline=T),
        uiOutput("target_ui"), 
        uiOutput("target_per_ui"),
        radioButtons("read",
            label="Sequencing machine type",
            choices=c("Illumina/Short-read" = "short",
                    "PacBio/ONT/long-read" = "long"),
            inline=T),
        uiOutput("read_ui"),
        numericInput("length",
            label="Average Read length",value=150,min=1,max=50000),
        conditionalPanel(condition = "input.target == 'No'",
                                 textOutput('coverage_text')),
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
        plotlyOutput("rarePlot"),
        tableOutput('table')
    )
  )
)
# Define server logic required to draw a histogram
server=function(input,output,session) {

  output$coverage_text <- renderText({
    if(input$target == "No") {
      paste("The assumed target is precent at 1% of metagenome in the sample 
      and its's genome size is 3 Mbp.")
    }
  })

  output$source_ui <- renderUI({
    if(input$source=="host") {
      sliderInput("contam", 
                  "Host DNA Contamination %", 
                  value=50, min=0, max=100)
    } 
  })
  output$target_ui <- renderUI({
    if(input$target=="Yes") {
      numericInput("gen_size", 
                   "Target Genome Size in MB", 
                   value=3, min=0, max=670000)
    } 
  })
  output$target_per_ui <- renderUI({
    if(input$target=="Yes") {
      sliderInput("gen_perc", 
                  "Target Organism Percentage in Sample %",
                  value=1, min=0, max=100)
    } 
  })
  output$read_ui <- renderUI({
    if(input$read=="short") {
      radioButtons("read_type", 
                   "Read Type", 
                   choices=c("Paired-end" = "PE", "Single-end" = "SE"), 
                   inline=T)
    }
  })

  in_values <- reactiveValues(contam = 0, 
                              gen_size = 3, 
                              gen_perc = 1, 
                              read_type = "SE")

  observeEvent(input$contam, {
    in_values$contam <- input$contam
  })
  observeEvent(input$gen_size, {
    in_values$gen_size <- input$gen_size
  })
  observeEvent(input$gen_perc, {
    in_values$gen_perc <- input$gen_perc
  })
  observeEvent(input$read_type, {
    in_values$read_type <- input$read_type
  })
  observeEvent(input$source, {
    if(input$source == "host") {
      in_values$contam <- input$contam
    } else {
      in_values$contam <- 0
    }
  })
    observeEvent(input$read, {
    if(input$read == "short") {
      in_values$read_type <- input$read_type
    } else {
      in_values$read_type <- "SE"
    }
  })


  observeEvent(input$click, {
    
    if(in_values$contam > 0) {
      fin_gen_perc <- in_values$gen_perc * ((100 - in_values$contam) / 100)
    } else {
      fin_gen_perc <- in_values$gen_perc
    }

    if(in_values$gen_perc < 10) {
      rel_abund_range <- seq(0, 10, by=0.25)
    } else {
      rel_abund_range <- seq(10, 100, by=1)
    }

    if(input$read == "short") {
      seq_range <- seq(10, 100, by=10)
    } else {
      seq_range <- seq(1, 10, by=1)
    }

    plot_data <- data.frame(rel_abund = numeric(),
                            reads = numeric(),
                            coverage = numeric())

    for(seq in seq_range) {
      seq_m <- seq * 1000000
      if(input$read == "short") {
        if(in_values$read_type == "PE") {
          total_bp <- total_basepairs(input$length, seq_m, paired=TRUE)
        } else {
          total_bp <- total_basepairs(input$length, seq_m, paired=FALSE)
        }
      } else {
        total_bp <- total_basepairs(input$length, seq_m, paired=FALSE)
      }
      for(rel in rel_abund_range){
        rel_contam <- rel * ((100 - in_values$contam) / 100)
        coverage <- (total_bp / (in_values$gen_size * 1e6)) * (rel_contam / 100)
        plot_data[nrow(plot_data)+1,] <- c(rel, seq_m, coverage)
      }
    }
    plot_data <- plot_data %>%
      mutate(reads = paste0(reads / 1000000, "M"))
    plot_data$reads <- factor(plot_data$reads, 
                                 levels = unique(plot_data$reads))
    
      high_contrast_12 <- c(
      "#E69F00",  # Orange
      "#56B4E9",  # Sky Blue
      "#009E73",  # Bluish Green
      "#F0E442",  # Yellow
      "#0072B2",  # Blue
      "#D55E00",  # Vermillion
      "#CC79A7",  # Reddish Purple
      "#882255",  # Dark Red
      "#44AA99",  # Teal
      "#117733",  # Green
      "#332288",  # Indigo
      "#AA4499"   # Pink
    )

    coverage_plot <- ggplot(plot_data, aes(x=rel_abund, y=coverage, color=reads)) +
      geom_point() +
      scale_y_log10() +
      scale_color_manual(values = high_contrast_12, name = "Reads in millions") +
      labs(x="Species relative abundance in metagenome (%)", y="Genome coverage (X)") +
      theme_minimal()  

    output$rarePlot <- renderPlotly({
      ggplotly(coverage_plot)
    })

    Settings <- c(
      "Metagenome Source" = input$source,
      "Host DNA Contamination %" = in_values$contam,
      "Target Organism Genome Size (Mb)" = in_values$gen_size,
      "Target Organism Percentage in Sample (%)" = in_values$gen_perc,
      "Sequencing Machine Type" = input$read,
      "Read Type" = in_values$read_type,
      "Average Read Length (bp)" = input$length
    )
    settings_df <- Settings %>% as.data.frame() %>%
        rownames_to_column("Setting") 
    colnames(settings_df) <- c("Setting", "Value")
    tab1 <- ggtexttable(settings_df, rows = NULL, 
                        theme = ttheme("mGreen"))

    output$table <- renderTable({
      settings_df
    }, rownames = FALSE, digits = 2)

    # Download handler for the plot
    output$download <- downloadHandler(
            filename ="coverage_sequencing.pdf",
            content = function(file){
                pdf(NULL)
                #plot_grid(p1, p2, p3, p4, nrow = 2)
                ggsave(file, plot=ggarrange(coverage_plot, tab1, nrow = 2), dpi = 300, width = 10, height = 8, units = "in")
            }
        )
    
  })
}
# Run the application
shinyApp(ui = ui, server = server)