library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(ggplot2)
library(tidyverse)
library(shinyhelper) # For adding help tooltips

options(rsconnect.max.bundle.size=8589934592)

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

abundance_df <- left_join(sample_abundance_lean, soil_metadata)
abundance_df <- left_join(abundance_df, 
                          organism_df_to_subset %>% select(taxonomy_id, taxon), 
                          relationship = "many-to-many") 

# Read in biome data file
biome_info = fread("./reference_data/nlcd_key.csv", drop = 1, header = T) 
biome_choices = biome_info$nlcdClass
names(biome_choices) = biome_info$prettyNlcd

abundance_df$biome = biome_info[match(abundance_df$nlcdClass, biome_info$nlcdClass),]$prettyNlcd



# UI
ui <- fluidPage(
    titlePanel("SoilMicrobeDB: An Interactive Database of Soil Microbial Genomes"),
    
    tags$p("The SoilMicrobeDB is a collection of over 30,000 soil microbial genomes, each annotated with ecological preferences for environmental conditions such as pH, temperature, and biome type. This tool allows you to filter, analyze, and visualize data on microbial species across different soil environments, using the sample collection from the National Ecological Observatory Network (NEON)."),
    
    tags$h4("Filters"),
    fluidRow(
        column(4, 
               selectInput("biome", "Select Biome", choices = unique(organism_df_to_subset$biome), multiple = TRUE) %>%
                   shinyhelper::helper(
                       type = "inline", 
                       title = "Biome Preference", 
                       content = "Biome preference is assigned if a taxon is present in at least 2% of samples within the biome. Biomes are assigned to each NEON sample using the National Land Cover Database.",
                       icon = "question-circle"
                   )
               
               
        ),
        column(4, 
               sliderInput("pH_range", "pH Preference Range", min = 3, max = 9, value = c(3, 9)) %>% 
                   shinyhelper::helper(
                       type = "inline",
                       title = "pH Preference",
                       content = "pH preference of each taxon is assigned as the peak of a LOESS curve fit to abundance data across the range of pH values. Click on a taxon to visualize or download this data.",
                       icon = "question-circle"
                   )
        ),
        column(4, 
               sliderInput("temperature_range", "Temperature Preference Range", min = 0, max = 40, value = c(0, 40)) %>%
                   shinyhelper::helper(
                       type = "inline",
                       title = "Temperature Preference",
                       content = "Temperature preference of each tacon is assigned as the peak of a LOESS curve fit to abundance data across the range of temperature values. Click on a taxon to visualize or download this data.",
                       icon = "question-circle"
                   )
        )
    ),
    
    tags$h4("Organism Data Table"),
    DT::dataTableOutput("organism_table"),
    downloadButton("download_organism", "Download Taxon List"),
    
    uiOutput("modal_abundance_plot")
)

# Server
server <- function(input, output, session) {
    
    # Initialize shinyhelper
    shinyhelper::observe_helpers(withMathJax = TRUE)
    
    # Reactive Filtered Organism DataFrame
    filtered_organism_df <- reactive({
        organism_df_to_subset %>%
            filter(
                (is.null(input$biome) || biome %in% input$biome),
                between(pH_preference, input$pH_range[1], input$pH_range[2]),
                between(temperature_preference, input$temperature_range[1], input$temperature_range[2])
            )
    })
    
    # Display Organism Data Table
    output$organism_table <- DT::renderDataTable({
        filtered_organism_df() %>% select(taxonomy_id,Kingdom, Genus, "Species of interest", "Genome source",accession) #, "Functional in COMETS?")
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
            # Render plots if data is available
            modalDialog(
                size = "l",
                title = paste("Abundance of", selected_taxon_name, "in NEON soil samples"),
                plotOutput("pH_plot"),
                plotOutput("temperature_plot"),
                tags$p(paste("Match Criteria:", selected_row$`Match criteria`)),
                tags$p("NCBI genome accession:", ifelse(!is.na(selected_row$accession), selected_row$accession, "N/A")),
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
        
        ggplot(abundance_data, 
               aes(x = pH, y = abundance)) +#, color=nlcdClass)) +
            geom_point(aes(color=biome),
                       alpha=.5, 
                       position=position_jitter(width = .01, height=0), size=2) + 
            #geom_smooth(method = "loess",show.legend = F, span=.7) +
            geom_smooth(method="gam",,
                        #method = "loess",
                        show.legend = F, se=F) +
            theme_bw(base_size = 18) + 
            #scale_y_sqrt() + 
            xlab("Soil pH") + 
            labs(title = paste("Abundance vs. pH for", selected_taxon_name)) +
            ylab("Microbial abundance")
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
               aes(x = temperature, y = abundance#, color=nlcdClass
                   )) +
            geom_point(aes(color=biome),
                       alpha=.5, 
                       position=position_jitter(width = .01, height=0), size=2) + 
           # geom_smooth(method = "loess", show.legend = F, span=.7) +
            geom_smooth(method="gam",
                #method = "loess", 
                show.legend = F, se=F) +
            theme_bw(base_size = 18) + 
            #scale_y_sqrt()  + 
            xlab("Soil temperature") +
            ylab("Microbial abundance") +
            labs(title = paste("Abundance vs. temperature for", selected_taxon_name))
        
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

