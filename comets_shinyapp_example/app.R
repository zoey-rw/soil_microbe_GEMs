library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(ggplot2)
library(tidyverse)
library(conflicted)
library(shinyhelper) # For adding help tooltips

# Use conflicted to set preferences for conflicting functions
conflict_prefer("between", "data.table")
conflict_prefer("filter", "dplyr")

# Load Data
abundance_data_filt <- fread("./species_abundance_filt.csv", nThread = 8, drop = 1, header = TRUE)
df_to_subset <- fread("./organism_data_to_subset.csv", drop = 1, header = TRUE)
df_to_print <- fread("./organism_data_to_print.csv", drop = 1, header = TRUE)
biome_info <- fread("./nlcd_key.csv", drop = 1, header = TRUE)
taxonomy <- fread("./organism_taxonomy.csv", drop = 1, header = TRUE)

# Prepare Biome Choices
biome_choices <- biome_info$nlcdClass
names(biome_choices) <- biome_info$prettyNlcd

# Define UI
ui <- fluidPage(
    titlePanel("SoilMicrobeDB: An Interactive Database of Soil Microbial Genomes"),
    tags$p("The SoilMicrobeDB is a collection of over 30,000 soil microbial genomes..."),

    tags$h4("Filters"),
    fluidRow(
        column(4,
               selectInput("biome", "Select Biome", choices = unique(df_to_subset$biome), multiple = TRUE) %>%
                   shinyhelper::helper(
                       type = "inline",
                       title = "Biome Preference",
                       content = "Biome preference is assigned if a taxon is present in at least 2%...",
                       icon = "question-circle"
                   )
        ),
        column(4,
               sliderInput("pH_range", "pH Preference Range", min = 3, max = 9, value = c(3, 9)) %>%
                   shinyhelper::helper(
                       type = "inline",
                       title = "pH Preference",
                       content = "pH preference of each taxon is assigned as the peak...",
                       icon = "question-circle"
                   )
        ),
        column(4,
               sliderInput("temperature_range", "Temperature Preference Range", min = 0, max = 40, value = c(0, 40)) %>%
                   shinyhelper::helper(
                       type = "inline",
                       title = "Temperature Preference",
                       content = "Temperature preference of each taxon...",
                       icon = "question-circle"
                   )
        )
    ),

    tags$h4("Organism Data Table"),
    DT::dataTableOutput("organism_table"),
    downloadButton("download_organism", "Download Taxon List"),

    uiOutput("modal_abundance_plot")
)

# Define Server
server <- function(input, output, session) {

    # Initialize shinyhelper
    shinyhelper::observe_helpers(withMathJax = TRUE)

    # Reactive Filtered Organism DataFrame
    filtered_organism_df <- reactive({
        df_to_subset %>%
            filter(
                (is.null(input$biome) || biome %in% input$biome),
                between(pH_preference, input$pH_range[1], input$pH_range[2]),
                between(temperature_preference, input$temperature_range[1], input$temperature_range[2])
            )
    })

    # Display Organism Data Table
    output$organism_table <- DT::renderDataTable({
        filtered_organism_df() %>% select(Kingdom, Genus, `Species of interest`, `Genome source`, `Functional in COMETS?`)
    }, selection = 'single')

    # Download Filtered Organism Data
    output$download_organism <- downloadHandler(
        filename = function() { "filtered_organism_data.csv" },
        content = function(file) { write.csv(filtered_organism_df(), file, row.names = FALSE) }
    )

    # Modal Abundance Plot
    output$modal_abundance_plot <- renderUI({
        req(input$organism_table_rows_selected)
        selected_row <- filtered_organism_df()[input$organism_table_rows_selected, ]
        selected_taxon <- selected_row$taxon

        # Filter abundance data for selected taxon
        abundance_filtered <- abundance_data_filt %>% filter(taxon == selected_taxon)

        if (nrow(abundance_filtered) == 0) {
            modalDialog(
                title = paste("Abundance Analysis for", selected_taxon),
                tags$p("Error: No abundance data available for this taxon."),
                footer = modalButton("Close")
            )
        } else {
            modalDialog(
                size = "l",
                title = paste("Abundance of", selected_taxon, "in NEON soil samples"),
                plotOutput("pH_plot"),
                plotOutput("temperature_plot"),
                footer = tagList(
                    downloadButton("download_filtered_abundance", "Download Taxon Abundance Data"),
                    modalButton("Close")
                )
            )
        }
    })

    # pH Plot
    output$pH_plot <- renderPlot({
        req(input$organism_table_rows_selected)
        selected_taxon <- filtered_organism_df()[input$organism_table_rows_selected, "taxon"]
        abundance_data <- abundance_data_filt %>% filter(taxon == selected_taxon)

        ggplot(abundance_data, aes(x = pH, y = abundance, color = biome)) +
            geom_point(alpha = .5, position = position_jitter(width = .01, height = 0), size = 2) +
            geom_smooth(method = "gam", show.legend = FALSE, se = FALSE) +
            theme_bw(base_size = 18) +
            xlab("Soil pH") + ylab("Microbial abundance") +
            labs(title = paste("Abundance vs. pH for", selected_taxon))
    })

    # Temperature Plot
    output$temperature_plot <- renderPlot({
        req(input$organism_table_rows_selected)
        selected_taxon <- filtered_organism_df()[input$organism_table_rows_selected, "taxon"]
        abundance_data <- abundance_data_filt %>% filter(taxon == selected_taxon)

        ggplot(abundance_data, aes(x = temperature, y = abundance, color = biome)) +
            geom_point(alpha = .5, position = position_jitter(width = .01, height = 0), size = 2) +
            geom_smooth(method = "gam", show.legend = FALSE, se = FALSE) +
            theme_bw(base_size = 18) +
            xlab("Soil temperature") + ylab("Microbial abundance") +
            labs(title = paste("Abundance vs. temperature for", selected_taxon))
    })

    # Download Filtered Abundance Data
    output$download_filtered_abundance <- downloadHandler(
        filename = function() { "taxon_abundance_data.csv" },
        content = function(file) {
            selected_taxon <- filtered_organism_df()[input$organism_table_rows_selected, "taxon"]
            abundance_data <- abundance_data_filt %>% filter(taxon == selected_taxon)
            write.csv(abundance_data, file, row.names = FALSE)
        }
    )
}

# Run the Application
shinyApp(ui, server)
