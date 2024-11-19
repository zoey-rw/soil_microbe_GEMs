library(tidyverse)
library(broom)

setwd("/projectnb/frpmars/soil_microbe_db/")
source("/projectnb/frpmars/soil_microbe_db/scripts/helper_functions.r")

set.seed(1)
# Read in microbe abundances, soil metadata, and GEM information
species_abundances <- fread("./soil_microbe_db_abundances.csv")  
# length(unique(species_abundances$taxonomy_id)) # 26204 unique taxa
#%>% 
	#select(-`...1`) %>% 
	# Remove some taxon categories that we probably don't want (no viruses)
#	filter(!taxon %in% c("Classified at a higher level","Unclassified")) %>% 
#	filter(!grepl("virus", taxon))

soil_metadata <- fread("../ref_data/environmental_metadata_NEON.csv") %>% 
	dplyr::select(-V1) %>% 
    mutate(sample_id = paste0(genomicsSampleID, "_soil_microbe_db_filtered")) %>% 
    filter(sample_id %in% species_abundances$sample_id)

# Combine into one dataframe 
sample_abundance = merge(species_abundances %>% filter(!is.na(source)), 
									soil_metadata, by.x = "sample_id", by.y = "sample_id") 
sample_abundance$db_name = "soil_microbe_db"
sample_abundance$taxon=sample_abundance$name


organism_info = sample_abundance %>% distinct(source, taxonomy_id, is_MAG, taxid_lineage, lineage, name, taxon)
write_csv(organism_info, "./intermediate_data/organism_info.csv")


# Filter species by prevalence in at least 50 samples
n_sample_prevalence = sample_abundance %>% group_by(name) %>% tally(name = "n_samples")
sample_abundance <- left_join(sample_abundance, n_sample_prevalence)
sample_abundance <- sample_abundance %>% filter(n_samples > 50)
length(unique(sample_abundance$taxonomy_id)) # 18589 unique taxa
write_csv(sample_abundance, "./sample_abundance_data.csv")


sample_abundance_lean = sample_abundance %>% dplyr::select(taxonomy_id, percentage, sample_id)
write_csv(sample_abundance_lean, "./intermediate_data/sample_abundance_data_lean.csv")

env_data_lean = sample_abundance %>% dplyr::select(sample_id)
write_csv(sample_abundance_lean, "./intermediate_data/sample_abundance_data_lean.csv")


# Establish 
merged_df_nest = merged_df %>% 
	#select(-c(sampleID)) %>% 
	group_by(db_name, taxon, taxonomy_id, lineage) %>% nest()

filter_by_correlation = F

if (filter_by_correlation==T){

# Filter by correlation strength
data_nest_pH <- merged_df_nest %>% 
	mutate(model = map(data, cor_fun_pH)) %>% 
	select(-data) %>% unnest(cols = c(model)) %>% 
	filter(p.value < 0.0001 & abs(estimate) > .8 )
species_of_interest_pH <- data_nest_pH$taxon

data_nest_temperature <- merged_df_nest %>% 
	mutate(model = map(data, cor_fun_temperature)) %>% 
	select(-data) %>% unnest(cols = c(model)) %>% 
	filter(p.value < 0.0001 & abs(estimate) > .3)
species_of_interest_temperature <- data_nest_temperature$taxon

# species_of_interest <- species_of_interest_carbon
# species_of_interest <- species_of_interest_pH
# species_of_interest <- species_of_interest_nitr

species_of_interest = c(species_of_interest_pH, species_of_interest_temperature) %>% 
	unique()

} else species_of_interest = unique(merged_df_nest$taxon) %>% sample(5000)

# Assign biome preferences for species with pH or temperature correlation
data_nest_biome <- merged_df_nest %>% filter(taxon %in% species_of_interest) 

# Do this in chunks otherwise the map() function freezes
df_biome_pt1 <- data_nest_biome[1:1000,] %>% 
	mutate(biome_presence = map(data, assign_biome_presence)) 

df_biome_pt2 <- data_nest_biome[1001:3000,] %>% 
	mutate(biome_presence = map(data, assign_biome_presence)) 

df_biome_pt3 <- data_nest_biome[3001:5000,] %>% 
	mutate(biome_presence = map(data, assign_biome_presence)) 

# df_biome_pt4 <- data_nest_biome[5001:nrow(data_nest_biome),] %>% 
# 	mutate(biome_presence = map(data, assign_biome_presence)) 

#df_biome = df_biome_pt1
df_biome = data.table::rbindlist(list(df_biome_pt1, df_biome_pt2, df_biome_pt3)) %>% 
	select(-data) %>% unnest(cols = c(biome_presence)) %>% ungroup 

df_biome_wide = df_biome %>% 
	select(-c(biome_count, prevalence)) %>% 
	pivot_wider(names_from = nlcdClass, values_from = present_in_biome)



# Loop through GEM data and match by species or genus
gem_out = list()
for (i in 1:length(species_of_interest)){
	print(species_of_interest[[i]])
	spec = species_of_interest[[i]]
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
	}
	match_out <- match_out %>% select(-Species)
	match_out$Species_of_interest = spec
gem_out[[i]] = match_out
}

gem_out_table = data.table::rbindlist(gem_out)

# gem_out_table$pH_significant = ifelse(gem_out_table$Species_of_interest %in% species_of_interest_pH, "Yes", "No")
# gem_out_table$temp_significant = ifelse(gem_out_table$Species_of_interest %in% species_of_interest_temperature, "Yes", "No")


functional_gems <- c("hanpo","iAA1300","iAF692","iFC579","iHN637","iJN7462","iRi1574","iSB619","iYO844","iGC535","iBM3063","iBB1018","iMG746","iMM904","iNmo686","iYY1101","iGD1348","iCY1106","iRhto1108")
gem_out_table$Functional = ifelse(gem_out_table$GEM_ID %in% functional_gems, "Yes", "No")

gem_out_table <- left_join(gem_out_table, df_biome_wide, join_by(Species_of_interest == taxon))

gem_out_table = gem_out_table %>% 
    filter(db_name == "soil_microbe_db") %>% 
	#arrange(Functional) %>% 
	#	arrange(GEM_ID)
	rename("Species of interest" = Species_of_interest, 
				 "Match criteria" = match_by, 
				 "GEM ID" = GEM_ID, 
				 #Filepath, 
				 "Functional in COMETS?" = Functional#, 
				 #"Temperature sensitive" = temp_significant,
				 #"pH sensitive" = pH_significant
				 ) %>% 
	arrange(desc(`Functional in COMETS?`)) 
#formattable::formattable(gem_out_print)

gem_out_table <- left_join(gem_out_table, merged_df %>% select(`Species of interest`=name,taxonomy_id, source) %>% distinct())

write.csv(gem_out_table, "./species_list_for_smartCOMETS.csv")

# write.csv(merged_df %>% filter(db_name == "soil_microbe_db"), 
# 					"/projectnb/frpmars/soil_microbe_db/shiny_app/species_abundance.csv")

# # Filter to just 5000 random organisms, to make data smaller for now
abundance_data_filt = merged_df %>%
	filter(taxon %in% species_of_interest)

write.csv(abundance_data_filt, "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/species_abundance_filt.csv")


# cp '/projectnb/frpmars/soil_microbe_db/Soil GEMs - manually curated.csv' /projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/
# cp /projectnb/frpmars/soil_microbe_db/species_list_for_smartCOMETS.csv /projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/

