# Prep dataframe for printing in web application

library(data.table)
library(tidyverse)
library(DT)

soil_GEM_info = read_csv("../Soil GEMs - manually curated.csv") %>% 
	select(GEM_ID, Kingdom, Genus, Species.strain = `Species/strain`, Citation, 
				 Source=`Source/link`, Filepath = `SCC Filepath`) %>% 
	mutate(Species = paste(Genus, word(Species.strain)))

organism_data = fread("./species_list_for_smartCOMETS.csv", nThread = 8, drop = 1, header = T)

# # Filter to just ~100 organisms with pH relationship, to make data smaller for now
organism_pH_signif = organism_data %>% filter(`pH sensitive` == "Yes")

# abundance_data = fread("../species_abundance_for_smartCOMETS_visualizations.csv", nThread = 8, drop = 1, header = T)
# abundance_data_filt = abundance_data %>% 
# 	filter(taxon %in% organism_pH_signif$`Species of interest`)
# 
# write_csv(abundance_data_filt, "./species_abundance_for_smartCOMETS_visualizations_filt.csv")

abundance_data_filt = fread("./species_abundance_for_smartCOMETS_visualizations_filt.csv", nThread = 8, drop = 1, header = T)

# Function to find the peak pH/temperature value
loess_maximum <- function(x, y, tol = .Machine$double.eps^0.5){
	model <- loess(y ~ x,span = 0.7)
	yfit <- model$fitted
	inx <- which(abs(yfit - max(yfit)) < tol)[[1]]
	list(x = x[inx], y.fitted = yfit[inx], ix = inx)
}

# Get pH of maximum abundance trends per taxon
pH_max = abundance_data_filt %>% 
	group_by(taxon) %>% 
	filter(!is.na(soilInCaClpH)) %>% 
	summarize(pH_preference = loess_maximum(soilInCaClpH, percentage)[[1]]) %>% 
	rename("Species of interest" = taxon)

# Get temp of maximum abundance trends per taxon
temp_max = abundance_data_filt %>% 
	group_by(taxon) %>% 
	filter(!is.na(soilTemp)) %>% 
	summarize(temperature_preference = loess_maximum(soilTemp, percentage)[[1]])  %>% 
	rename("Species of interest" = taxon)

df_to_print <- left_join(organism_pH_signif, pH_max)
df_to_print <- left_join(df_to_print, temp_max) %>% 
	select(`Species of interest`,pH_preference, temperature_preference) %>% mutate(genome_link = NA)



taxon_names = df_to_print$`Species of interest`


gem_out = list()
for (i in 1:length(taxon_names)){
	print(taxon_names[[i]])
	spec = taxon_names[[i]]
	genus = word(spec)
	if (spec %in% soil_GEM_info$Species) {
		match_out=soil_GEM_info[which(soil_GEM_info$Species==spec),]
		match_out$match_by="Species"
	} else if (genus %in% soil_GEM_info$Genus) {
		match_out=soil_GEM_info[which(soil_GEM_info$Genus==genus),]
		match_out$match_by="Genus"
	} else {
		match_out = data.frame(GEM_ID = NA,Kingdom=NA,Genus=NA,Species.strain=NA,
													 Citation=NA,Source=NA,Filepath=NA,
													 Species = NA,
													 match_by="No curated GEM at species or genus level")
		#match_out$DB_Species = species_of_interest
		#match_out$Genus = genus_of_interest
	}
	print(match_out)
	match_out <- match_out %>% select(-Species)
	match_out$`Species of interest` = spec
	gem_out[[i]] = match_out
}

gem_out_table = data.table::rbindlist(gem_out)



df_to_print <- left_join(df_to_print, gem_out_table) %>% 
	select(`Species of interest`,pH_preference, temperature_preference, match_by, GEM_ID)

write.csv(df_to_print, "./organism_data_to_print.csv")

