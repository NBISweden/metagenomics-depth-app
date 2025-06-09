library(shiny)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(shinythemes)
library(plotly)
library(ggpubr)
library(shinyhelper)

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
    theme = shinytheme("cosmo"),
  sidebarLayout(
    sidebarPanel(
        radioButtons("source",
            label="Metagenome Source",
            choices=c("Human/Host" = "host",
                    "Environmental" = "env"),
            inline=T) %>%
            helper(type = "inline", 
              content = c("Select the source of the metagenome.", 
                          "If the source is a host, like human skin metagenome or ostrich gut metagenome",
                          "If the source is environmental like soil, water or even bioreactor"),
              title = "Source", size = "s"),
        uiOutput("source_ui"),
        radioButtons("target",
            label="Is there a target organism?",
            choices=c("Yes","No"),
            inline=T) %>%
            helper(type = "inline", 
               content = c("Is there a particular organism of interest in the metagenome?",
                          "If yes, please provide the genome size and expected relative abundance in the sample.",
                          "If no, the calculator will assume a target organism with 3 Mbp genome size and 1% relative abundance."),
              title = "Target Organism", size = "s"),
        uiOutput("target_ui"), 
        uiOutput("target_per_ui"),
        radioButtons("read",
            label="Sequencing machine type",
            choices=c("Illumina/Short-read" = "short",
                    "PacBio/ONT/long-read" = "long"),
            inline=T),
        uiOutput("read_ui"),
        numericInput("length",
            label="Average Read length in bp",value=150,min=1,max=50000) %>%
            helper(type = "inline", 
               content = c("Illumina/short-read sequencing is usually 150 bp,",
                          "PacBio/ONT/long-read sequencing is usually 10,000 bp or more."),
              title = "Average Read Length", size = "s"),
        sliderInput("coverage",
            "X coverage aim for target species",
            min = 0,
            max = 100,
            value = 10) %>%
            helper(type = "inline", 
              content = c("This is the number of times the target genome is expected to be sequenced in the sample.",
                          "For example, if you expect 10X coverage, it means that each base in the target genome will be sequenced on average 10 times.",
                          "Naturally, any organism that is more abundant in the metagenome will receive a higher coverage on average"),
              title = "Expected Coverage", size = "s"),
        conditionalPanel(condition = "input.target == 'No'",
                                 textOutput('coverage_text')),
        fluidRow(
          column(6,style=list("padding-right: 5px;"),
                 actionButton("click", "Run test")
          ),
          column(6,style=list("padding-left: 5px;"),
                 downloadButton('download', 'Download')
          ),
        )
    ),
    
    mainPanel(
        plotlyOutput("rarePlot"),
        fluidRow(
                    column(width = 6, plotlyOutput("spePlot",  width="100%")),
                    column(width = 6, tableOutput('table'))
                )
    )
  )
)
# Define server logic required to draw a histogram
server=function(input,output,session) {

  observe_helpers()

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
                  value=50, min=0, max=100) %>%
            helper(type = "inline", 
              content = c("The amount of host DNA contamination like from human can vary widely", 
                          "Anything between <10% (like faecal) to >90% (like skin) depending on the location in human.",
                          'Related research on host contamination in human metagenome can be found <a href="https://www.nature.com/articles/s41564-023-01381-3">here</a>'),
              title = "Host DNA Contamination", size = "s")
    } 
  })
  output$target_ui <- renderUI({
    if(input$target=="Yes") {
      numericInput("gen_size", 
                   "Target Genome Size in MB", 
                   value=3.5, min=0, max=670000, step=0.5) %>%
            helper(type = "inline", 
              content = "For example an average size of a bacterial genome is 3.5 MB",
              title = "Target Genome Size", size = "s")
    } 
  })
  output$target_per_ui <- renderUI({
    if(input$target=="Yes") {
      sliderInput("gen_perc", 
                  "Target Organism Percentage in Sample %",
                  value=1, min=0, max=100) %>%
            helper(type = "inline", 
               content = "If there is an expected % of the target organism in the metagenome",
               title = "Target Organism Percentage", size = "s")
    } 
  })
  output$read_ui <- renderUI({
    if(input$read=="short") {
      radioButtons("read_type", 
                   "Read Type", 
                   choices=c("Paired-end" = "PE", "Single-end" = "SE"), 
                   inline=T) %>%
            helper(type = "inline", 
              content = c("Illumina sequencing is usually Paired-end (PE)", 
                        "while Long-read sequencing technologies are Single-end (SE)"),
              title = "Read Type", size = "s")
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
      x_seq <- seq(10, 100, by=1)
    } else {
      seq_range <- seq(1, 10, by=1)
      x_seq <- seq(1, 10, by=0.1)
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
    plot1_data <- plot_data %>%
      mutate(reads = paste0(reads / 1000000, "M"))
    plot1_data$reads <- factor(plot1_data$reads, 
                                 levels = unique(plot1_data$reads))
    
    plot2_data <- plot_data %>%
      filter(rel_abund == in_values$gen_perc) 
    
    if (min(plot2_data$coverage) > input$coverage) {
      intersections <- plot2_data %>%
        mutate(estimated_sequencing_depth = reads, target_coverage = min(plot2_data$coverage)) %>%
        slice_head()
      plot2_text <- paste0("Target already achieved at \n", paste0(intersections$estimated_sequencing_depth / 1000000, "M"), " reads.")
      
    } else if (max(plot2_data$coverage) < input$coverage) {
      intersections <-plot2_data %>%
        mutate(estimated_sequencing_depth = reads, target_coverage = max(plot2_data$coverage)) %>%
        slice_tail()
      plot2_text <- paste0("Target coverage \n NOT achieved \n even at ", paste0(intersections$estimated_sequencing_depth / 1000000, "M"), " reads.")
      
    } else {
      x_seq = seq(from=min(plot2_data$reads), to=max(plot2_data$reads), length.out=101)
      intersections <-plot2_data %>%
        group_by(rel_abund) %>%
        summarise(estimated_sequencing_depth = approx(x = reads, y = coverage, xout = x_seq)$y) %>%
        mutate(x_seq = x_seq) %>%
        slice_min(abs(estimated_sequencing_depth - input$coverage))
      
      names(intersections) <- c("rel_abund", "target_coverage", "estimated_sequencing_depth")
      plot2_text <- paste0("Target Achieved at \n", paste0(intersections$estimated_sequencing_depth / 1000000, "M"), " reads.")
    }
      
    
    high_contrast_12 <- c(
      "#FED976",  # Yellow
      "#FEB24C",  # Yellow-orange
      "#FD8D3C",  # Orange
      "#FC4E2A",  # Orange-red
      "#E31A1C",  # Red
      "#BD0026",  # Dark red
      "#800026",  # Very dark red
      "#67000D",  # Darkest red
      "#4D0009",  # Extra dark red
      "#330006"   # Nearly black red
    ) 
    coverage_plot <- ggplot(plot1_data, aes(x=rel_abund, y=coverage, color=reads)) +
      geom_point() +
      scale_y_log10() +
      scale_color_manual(values = high_contrast_12, name = "Sequencing depth\n   (million reads)") +
      labs(x="Species relative abundance in metagenome (%)", y="Genome coverage (X)") +
      theme_minimal()  

    target_plot <- ggplot(plot2_data, aes(x=reads, y=coverage)) +
      geom_line(aes(group=1), color="#CC79A7") +
      geom_segment(data=intersections, aes(x=estimated_sequencing_depth, y=target_coverage, xend=0, yend=input$coverage), 
                     linetype="dashed", color="darkgrey") +
      geom_segment(data=intersections, aes(x=estimated_sequencing_depth, y=target_coverage, xend=estimated_sequencing_depth, yend=1), 
                     linetype="dashed", color="darkgrey") +
      geom_point(data=intersections, aes(x=estimated_sequencing_depth, y=target_coverage), 
                 color="#CC79A7", size=3) +
      geom_text(data=intersections, aes(x=mean(plot2_data$reads)*1.5, y=mean(plot2_data$coverage)/2, label=plot2_text),
      vjust = 4, color="black") +
      theme_minimal() +
      scale_y_log10() +
      scale_x_continuous(breaks = plot2_data$reads, 
        label = paste0(plot2_data$reads / 1000000, "M")) +
      labs(x="Sequences per sample (million reads)", 
             y="Genome coverage (X)")

    output$rarePlot <- renderPlotly({
      ggplotly(coverage_plot)
    })

    output$spePlot <- renderPlotly({
      ggplotly(target_plot)
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
        rownames_to_column("Settings") 
    colnames(settings_df) <- c("Settings", "Value")
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
                ggsave(file, plot=ggarrange(coverage_plot, ggarrange(target_plot, tab1, nrow = 1), nrow = 2), dpi = 300, width = 10, height = 8, units = "in")
            }
        )
    
  })
}
# Run the application
shinyApp(ui = ui, server = server)