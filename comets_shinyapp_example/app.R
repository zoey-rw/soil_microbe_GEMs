library(shiny)
library(data.table)
library(tidyverse)
library(DT)
library(renv)

# Read in data file for plotting
# The column "taxon" matches with "Species of interest"
abundance_data_filt = fread("./species_abundance_filt.csv",
														nThread = 8, drop = 1, header = T)

# Read in data file for species selections - info on environmental preferences
df_to_subset = fread("./organism_data_to_subset.csv", drop = 1, header = T)

# Read in data file to just show species names/links
set.seed(1)
df_to_print = fread("./organism_data_to_print.csv", drop = 1, header = T)

# Read in biome data file
biome_info = fread("./nlcd_key.csv", drop = 1, header = T)
biome_choices = biome_info$nlcdClass
names(biome_choices) = biome_info$prettyNlcd

# Read in taxonomy data file, add to organism info
taxonomy = fread("./organism_taxonomy.csv", drop = 1, header = T)
df_to_print = left_join(df_to_print, taxonomy %>%
                            select(`Species of interest` = taxon,`Genome accession` = accession) %>%
                            unique()) %>%
    select(`Species of interest`,GEM_ID, `Genome source`) %>% slice_sample(n = 1000)
df_to_subset = left_join(df_to_subset, taxonomy %>%
                             select(`Species of interest` = taxon,Kingdom,Phylum,accession) %>%
                             unique())  %>% slice_sample(n = 1000)

# Potential taxa to select from
taxon_names = df_to_subset$`Species of interest`

# Create column of radio buttons
mat_to_print= as.matrix(df_to_print[,1])
rownames(mat_to_print) = mat_to_print[,1]

for (i in seq_len(nrow(mat_to_print))) {
	mat_to_print[i, ] = sprintf(
		'<input type="radio" name="%s" value="%s"/>',
		taxon_names[i], mat_to_print[i, ]
	)
}
mat_to_print = cbind(mat_to_print, df_to_print[,1:3])
colnames(mat_to_print)[1] = "Visualize?"


# Filter abundance data to create example plots for selected taxon
single_species_obs = abundance_data_filt %>%
	filter(taxon == "Rhodotorula toruloides")

# Visualize abundances that correlate with pH
p1 = ggplot(single_species_obs,
						aes(x = soilInCaClpH, y = percentage, color=taxon)) +
	geom_point(alpha=.5,
						 position=position_jitter(width = .01, height=0), size=2,
						 show.legend = F) +
	geom_smooth(show.legend = F, span=.7) +
	facet_wrap(~taxon, scales = "free") + theme_bw(base_size = 20) +
	scale_y_sqrt() + xlab("Soil pH") +
	ylab("Microbial abundance")


# Visualize abundances that correlate with temperature
p2 = ggplot(single_species_obs,
						aes(x = soilTemp, y = percentage, color=taxon)) +
	geom_point(alpha=.5,
						 position=position_jitter(width = .01, height=0), size=2,
						 show.legend = F) +
	geom_smooth(show.legend = F, span=.7) +
	facet_wrap(~taxon, scales = "free") + theme_bw(base_size = 20) +
	scale_y_sqrt()  + xlab("Soil temperature") +
	ylab("Microbial abundance")


# Actual interface setup
ui <- fluidPage(
	titlePanel("Explore the Soil Microbe Database"),
	mainPanel(

		fluidRow(p("The Soil Microbe Database (SMDB) is a collection of over 30,000 soil microbial genomes, some of which have cultured representatives. Use this portal to explore how the abundance of each genome varies across soil samples measured via shotgun metagenomics. This abundance dataset, and the SMDB, can be downloaded for further analysis using the links at the bottom of this page."),
			p("Use the filters below to identify soil microbes that are observed to peak in abundance at specific pH or temperature values. Note that these reflect trends in soils derived from sequencing, not laboratory experiments on actual pH or temperature tolerances. For more information on environmental abundances of microbes: https://doi.org/10.1111/nph.17240"),
			column(width = 4,
	sliderInput("pHrange", "Realized soil pH preference:",min = 3, max = 9, value = c(3,9))),
	column(width = 4,
	sliderInput("temperatureRange",
							"Realized soil temperature preference:",min = 0, max = 100, value = c(0,100))),
	column(width = 4,
				 textInput("taxonName", "Filter by species taxonomy instead"))),
	fluidRow(column(width=8,
	                checkboxGroupInput("biomeSelect", "Biome", biome_choices, selected = biome_choices, inline = TRUE))),
	fluidRow(
		p("All species within filters are listed below. Visualize one species at a time using the by selecting a species.")),

	fluidRow(	column(width = 4,
									 plotOutput("pH_plot")),
						column(width = 4,
									 plotOutput("temp_plot")),
						column(width = 4,
						       uiOutput("GEMtext1"),
						       uiOutput("GEMtext2"))

						#textOutput("GEMtext"))
						),
fluidRow(column(width=12,
	DT::dataTableOutput('print_table'))
)
)
)

# Server side
server <- shinyServer(function(input, output, session){


	#closest_GEM <- reactive({ # This value is not currently reactive!
	#     #		  "No curated GEM at species or genus level",
	# 	#organism_data[input$taxon,]
	# 	("Selected species: Chitinophaga pinensis", tags$br(),
	# 	 "Culture status: Cultured, with strain and media information in Bacdive",
	# 	 tags$br(),
	# 	 "The closest available species with a COMETS simulation-ready model is iRhto1108, matched by species name")
	#
	# })
	    url <- a("Download model here", href="https://github.com/zoey-rw/soil_microbe_GEMs/tree/master/iRhto1880")

	  output$GEMtext1 <- renderUI({
	      HTML(paste0(	#organism_data[input$taxon,]
	          "Selected species: Chitinophaga pinensis",
	           tags$br(),tags$br(),
	           "Culture status: Cultured, with strain and media information in Bacdive.",
	           tags$br(),tags$br(),
	           "The closest available species with a COMETS simulation-ready model is iRhto1108, matched by species name."), collapse = "<br>")
	      })

	      #closest_GEM())
	output$GEMtext2 <- renderUI({url})


#	output$GEMtext <- renderText({#closest_GEM()})

	output$pH_plot <- renderPlot({p1})
	output$temp_plot <- renderPlot({p2})
	output$print_table <- renderTable({df_to_print})


	output$print_table = DT::renderDataTable(
		mat_to_print, escape = FALSE, selection = 'none', server = FALSE,
		options = list(dom = 't', paging = FALSE, ordering = FALSE),
		callback = JS("table.rows().every(function(i, tab, row) {
          var $this = $(this.node());
          $this.attr('id', this.data()[0]);
          $this.addClass('shiny-input-radiogroup');
        });
        Shiny.unbindAll(table.table().node());
        Shiny.bindAll(table.table().node());")
	)
	output$sel = renderPrint({
		str(sapply(taxon_names, function(i) input[[i]]))
	})
})

options(rsconnect.max.bundle.size=3145728000)
shinyApp(ui = ui, server = server)

#rsconnect::deployApp(appName = "soil_microbe_db", appDir = "/projectnb/frpmars/soil_microbe_db/shiny_app")