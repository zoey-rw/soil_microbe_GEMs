
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(tidyverse)
library(shinyhelper) # For adding help tooltips


options(rsconnect.max.bundle.size=8589934592)



# Read in biome selection table 
biome_info = fread("./reference_data/nlcd_key.csv", drop = 1, header = T) 
biome_choices = biome_info$nlcdClass
names(biome_choices) = biome_info$prettyNlcd

# Read in data file for species selections - info on environmental preferences
organism_df_to_subset <- fread("./intermediate_data/organism_data_to_subset.csv", drop = 1)  %>% 
    mutate(taxon=`Species of interest`)

# Read in taxonomy and add to organism info
taxonomy_info = fread("./intermediate_data/organism_taxonomy.csv")
organism_df_to_subset <- left_join(organism_df_to_subset, taxonomy_info %>% 
                                       select(taxonomy_id, Kingdom, Phylum, Genus))


# Read in data file for species abundances, then merge in other sample information
sample_abundance_lean <- fread("./intermediate_data/sample_abundance_data_lean.csv") %>% 
    rename("abundance"="percentage")

# Sample information
soil_metadata <- fread("./reference_data/environmental_metadata_NEON.csv") %>% 
    dplyr::select(-V1) %>% 
    mutate(sample_id = paste0(genomicsSampleID, "_soil_microbe_db_filtered")) %>% 
    filter(sample_id %in% sample_abundance_lean$sample_id) %>% 
    rename("temperature"="soilTemp","pH" = "soilInCaClpH")

# Convert to dtplyr so that this join doesn't overwhelm RAM
sample_abundance_dt <- lazy_dt(sample_abundance_lean)
abundance_dt <- left_join(sample_abundance_dt, soil_metadata)
abundance_dt <- left_join(abundance_dt, 
                          organism_df_to_subset %>% select(taxonomy_id, taxon), 
                          relationship = "many-to-many")
 
# Back to data frame for the rest of the code
abundance_df <- abundance_dt %>% as.data.frame()
abundance_df$biome = biome_info[match(abundance_df$nlcdClass, biome_info$nlcdClass),]$prettyNlcd

#print(dim(organism_df_to_subset))
#print(dim(abundance_df))



# UI
ui <- fluidPage(
    # Accessibility: declare page language and skip-link target.
    tags$html(lang = "en"),
    tags$a(href = "#main-content", class = "sr-only sr-only-focusable",
           "Skip to main content"),

    titlePanel("SoilMicrobeDB: An Interactive Database of Soil Microbial Genomes"),

    tags$main(id = "main-content",
        tags$p("The SoilMicrobeDB is a collection of over 30,000 soil microbial genomes, each annotated with ecological preferences for environmental conditions such as pH, temperature, and biome type. This tool allows you to filter, analyze, and visualize data on the most widespread microbial species across different soil environments, using the sample collection from the National Ecological Observatory Network (NEON)."),

        tags$h2("Filters"),
        fluidRow(
            column(4,
                   selectInput("biome", "Select Biome", choices = biome_choices, multiple = TRUE, selected = biome_choices) %>%
                       shinyhelper::helper(
                           type = "inline",
                           title = "Biome Preference",
                           content = "Biome preference is assigned if a taxon is present in at least 2% of samples within the biome. Biomes are assigned to each NEON sample using the National Land Cover Database.",
                           icon = "question-circle",
                           buttonLabel = "Help: Biome Preference"
                       )
            ),
            column(4,
                   sliderInput("pH_range", "pH Preference Range", min = 3, max = 9, value = c(3, 9)) %>%
                       shinyhelper::helper(
                           type = "inline",
                           title = "pH Preference",
                           content = "pH preference of each taxon is assigned as the peak of a LOESS curve fit to abundance data across the range of pH values. Click on a taxon to visualize or download this data.",
                           icon = "question-circle",
                           buttonLabel = "Help: pH Preference"
                       )
            ),
            column(4,
                   sliderInput("temperature_range", "Temperature Preference Range", min = 0, max = 40, value = c(0, 40)) %>%
                       shinyhelper::helper(
                           type = "inline",
                           title = "Temperature Preference",
                           content = "Temperature preference of each taxon is assigned as the peak of a LOESS curve fit to abundance data across the range of temperature values. Click on a taxon to visualize or download this data.",
                           icon = "question-circle",
                           buttonLabel = "Help: Temperature Preference"
                       )
            )
        ),

        tags$h2("Organism Data Table"),
        tags$p(class = "sr-only",
               "Use the table below to explore filtered organisms. ",
               "Use the keyboard arrow keys to navigate rows, ",
               "and press Enter or Space to open the abundance modal."),
        DT::dataTableOutput("organism_table"),
        downloadButton("download_organism", "Download Taxon List"),

        uiOutput("modal_abundance_plot")
    )  # end tags$main
)

# Server
server <- function(input, output, session) {
    
    # Initialize shinyhelper
    shinyhelper::observe_helpers(withMathJax = TRUE)
    
    # Reactive Filtered Organism DataFrame
    filtered_organism_df <- reactive({
        # Get selected biomes
        selected_biomes <- input$biome
        
        # If no biomes are selected, return an empty data frame
        if (is.null(selected_biomes) || length(selected_biomes) == 0) {
            return(organism_df_to_subset[0, ])  # Return an empty data frame with the same structure
        }
        
        # Filter organism_df_to_subset based on selected biomes
        filtered_df <- organism_df_to_subset %>%
            filter(
                # Check if the selected biome column has a value of 1
                rowSums(dplyr::select(., all_of(selected_biomes)) == 1, na.rm = TRUE) > 0,
                between(pH_preference, input$pH_range[1], input$pH_range[2]),
                between(temperature_preference, input$temperature_range[1], input$temperature_range[2])
            )
        
        return(filtered_df)
    })
    
    # Display Organism Data Table
    output$organism_table <- DT::renderDataTable({
        filtered_organism_df() %>% dplyr::select(taxonomy_id,Kingdom, Genus, "Species of interest", "Genome source",accession) #, "Functional in COMETS?")
    }, selection = 'single')
    
    # Download Filtered Organism Data
    output$download_organism <- downloadHandler(
        filename = function() { paste("filtered_organism_data.csv") },
        content = function(file) { write.csv(filtered_organism_df(), file, row.names = FALSE) }
    )
    
    # Modal Abundance Plot
    output$modal_abundance_plot <- renderUI({
        req(input$organism_table_rows_selected)
        selected_row <- filtered_organism_df()[input$organism_table_rows_selected, ]
        selected_taxon <- selected_row$taxonomy_id
        selected_taxon_name <- selected_row$`Species of interest`
        
        
        # Filter abundance_df for selected taxon and check for data
        abundance_filtered <- abundance_df %>%
            filter(taxonomy_id == selected_taxon) #%>%
        #filter(taxon == selected_taxon)
        
        if (nrow(abundance_filtered) == 0) {
            # Show error message if no data is available for the selected taxon
            modalDialog(
                title = paste("Abundance Analysis for", selected_taxon_name),
                tags$p("Error: No abundance data available for this taxon."),
                footer = modalButton("Close")
            )
        } else {
            # Build the GEM-match summary string OUTSIDE the modalDialog argument
            # list. The previous version did `gem_match <- ifelse(...)` as one
            # of the modalDialog() args, which silently leaked the assignment
            # value into the dialog body (the "yet it does print?" comment).
            gem_match <- if (is.na(selected_row$GEM_ID)) {
                "No curated GEM at species or genus level."
            } else {
                paste0("Genome-scale model (GEM) available: ", selected_row$GEM_ID,
                       ". Match criteria: ", selected_row$`GEM match criteria`, ".")
            }

            # Render plots if data is available
            modalDialog(
                size = "l",
                title = paste("Abundance of", selected_taxon_name, "in NEON soil samples"),
                plotOutput("pH_plot"),
                plotOutput("temperature_plot"),
                tags$p(gem_match),
                tags$p("NCBI genome accession: ",
                       ifelse(!is.na(selected_row$accession), selected_row$accession, "N/A")),
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
        
        selected_row <- filtered_organism_df()[input$organism_table_rows_selected, ]
        selected_taxon <- selected_row$taxonomy_id
        selected_taxon_name <- selected_row$`Species of interest`
        
        abundance_data <- abundance_df %>% filter(taxonomy_id == selected_taxon)
        
        if (nrow(abundance_data) == 0) {
            return(NULL) # Avoid rendering if no data
        }
        
        # Encode biome with both color AND shape so colorblind users can
        # distinguish points without relying on hue alone.
        ggplot(abundance_data,
               aes(x = pH, y = abundance, color = biome, shape = biome)) +
            geom_point(alpha = .5,
                       position = position_jitter(width = .01, height = 0), size = 2) +
            geom_smooth(method = "gam", show.legend = F, se = F) +
            scale_shape_manual(values = rep(c(16, 17, 15, 3, 7, 8, 18, 4),
                                            length.out = length(unique(abundance_data$biome)))) +
            theme_bw(base_size = 18) +
            xlab("Soil pH") +
            labs(title = paste("Abundance vs. pH for", selected_taxon_name)) +
            ylab("Microbial abundance")
    }, alt = function() {
        # Plot description for screen readers (Shiny renderPlot's alt arg).
        sel <- input$organism_table_rows_selected
        if (is.null(sel)) return("Scatter plot of microbial abundance versus soil pH (no taxon selected).")
        name <- filtered_organism_df()[sel, ]$`Species of interest`
        paste0("Scatter plot of ", name, " abundance across soil pH values, ",
               "with one point per NEON sample, colored and shaped by biome.")
    })
    
    # Temperature Plot
    output$temperature_plot <- renderPlot({
        req(input$organism_table_rows_selected)
        
        selected_row <- filtered_organism_df()[input$organism_table_rows_selected, ]
        selected_taxon <- selected_row$taxonomy_id
        selected_taxon_name <- selected_row$`Species of interest`
        
        abundance_data <- abundance_df %>% filter(taxonomy_id == selected_taxon)
        
        if (nrow(abundance_data) == 0) {
            return(NULL) # Avoid rendering if no data
        }
        
        ggplot(abundance_data,
               aes(x = temperature, y = abundance, color = biome, shape = biome)) +
            geom_point(alpha = .5,
                       position = position_jitter(width = .01, height = 0), size = 2) +
            geom_smooth(method = "gam", show.legend = F, se = F) +
            scale_shape_manual(values = rep(c(16, 17, 15, 3, 7, 8, 18, 4),
                                            length.out = length(unique(abundance_data$biome)))) +
            theme_bw(base_size = 18) +
            xlab("Soil temperature") +
            ylab("Microbial abundance") +
            labs(title = paste("Abundance vs. temperature for", selected_taxon_name))

    }, alt = function() {
        sel <- input$organism_table_rows_selected
        if (is.null(sel)) return("Scatter plot of microbial abundance versus soil temperature (no taxon selected).")
        name <- filtered_organism_df()[sel, ]$`Species of interest`
        paste0("Scatter plot of ", name, " abundance across soil temperatures, ",
               "with one point per NEON sample, colored and shaped by biome.")
    })
    
    # Download Filtered Abundance Data
    output$download_filtered_abundance <- downloadHandler(
        filename = function() { paste("taxon_abundance_data.csv") },
        content = function(file) {
            selected_taxon <- filtered_organism_df()[input$organism_table_rows_selected, "taxonomy_id"]
            abundance_data <- abundance_df %>% filter(taxonomy_id == selected_taxon)
            
            if (nrow(abundance_data) > 0) {
                write.csv(abundance_data, file, row.names = FALSE)
            }
        }
    )
}

shinyApp(ui, server)

