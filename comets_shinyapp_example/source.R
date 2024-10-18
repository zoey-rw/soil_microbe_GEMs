# Prep dataframe for printing in web application
# firstrun "/projectnb/frpmars/soil_microbe_db/scripts/determine_species_of_interest.r"

library(data.table)
library(tidyverse)
library(DT)

set.seed(1)
# categories from https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend

soilData <- readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds")
soilCores <- soilData$sls_soilCoreCollection %>% 
    mutate(isForest = ifelse(grepl("Forest", nlcdClass), 
                             "forest habitat", "non forest habitat")) %>% 
    mutate(biome = recode(nlcdClass, "mixedForest" = "Forest",
                          "evergreenForest" = "Forest",
                          "deciduousForest" = "Forest",
                          "emergentHerbaceousWetlands" = "Wetlands",
                          "woodyWetlands" = "Wetlands",
                          "dwarfScrub" = "Shrubland",
                          "shrubScrub" = "Shrubland",
                          "sedgeHerbaceous" = "Herbaceous",
                          "grasslandHerbaceous" = "Herbaceous",
                          "pastureHay" = "Agricultural",
                          "cultivatedCrops" = "Agricultural"
    ))

genomicSamples <- soilData$sls_metagenomicsPooling %>%
	tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
	dplyr::select(genomicsSampleID,first,second,third)
genSampleExample <- genomicSamples %>%
	tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
	dplyr::select(sampleID,genomicsSampleID) %>%
	drop_na()
soilCores$compositeSampleID = genSampleExample[match(soilCores$sampleID,																					genSampleExample$sampleID),]$genomicsSampleID

soilCore_subset = soilCores %>% select(sampleID = compositeSampleID,plotID,biomassID,
                                       nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, 
                                       litterDepth, standingWaterDepth, soilTemp, 
                                       sampleTiming, elevation) %>% 
    distinct(sampleID, .keep_all = T) 

soil_GEM_info = read_csv("../Soil GEMs - manually curated.csv") %>% 
	select(GEM_ID, Kingdom, Genus, Species.strain = `Species/strain`, Citation, 
				 Source=`Source/link`, Filepath = `SCC Filepath`) %>% 
	mutate(Species = paste(Genus, word(Species.strain)))

# # Filter to just ~100 random organisms
organism_data = fread("./species_list_for_smartCOMETS.csv", nThread = 8, drop = 1, header = T)  %>% slice_sample(n = 100)

# # Filter to just ~100 organisms with pH relationship, to make data smaller for now
#organism_pH_signif = organism_data %>% filter(`pH sensitive` == "Yes")

abundance_data = fread("../species_abundance_filt.csv", nThread = 8, drop = 1, header = T)

abundance_data_filt = abundance_data %>% filter(taxon %in% organism_data$`Species of interest`)

#write_csv(abundance_data_filt, "./species_abundance_filt.csv")

abundance_data_filt = fread("./species_abundance_filt.csv", nThread = 8, drop = 1, header = T) 

# Function to find the peak pH/temperature value
loess_maximum <- function(x, y, tol = .Machine$double.eps^0.5){
	model <- loess(y ~ x,span = 0.8)
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

df_to_subset <- left_join(organism_data, pH_max)
df_to_subset <- left_join(df_to_subset, temp_max) %>% 
	#select(`Species of interest`,pH_preference, temperature_preference) %>% 
	mutate(genome_link = NA)



taxon_names = df_to_subset$`Species of interest`


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


df_to_subset <- left_join(df_to_subset, gem_out_table) 
df_to_subset = df_to_subset %>% distinct(.keep_all = T)

df_to_print <- df_to_subset %>% 
	select(`Species of interest`, match_by,
				 GEM_ID, "Genome source" = source,
				 #"Genome (link)" = genome_link, 
				 `Functional in COMETS?`) %>% 
    distinct(.keep_all = T)

write.csv(df_to_print, "./organism_data_to_subset.csv")
write.csv(df_to_print, "./organism_data_to_print.csv")



#lineage_df =  split_taxonomy_ncbi(bracken_with_lineage$lineage)
lineage_df =  lapply(organism_data$lineage, 
                     function(x) {
                         phyloseq::parse_taxonomy_qiime(x) %>% #as.data.frame() %>% 
                             t()  %>% as.data.frame()
                     }) %>% 
    data.table::rbindlist()

taxonomy_to_save = lineage_df %>% select(taxon=Species, #accession, 
                                         Kingdom, Phylum, Class, Order, Family, Genus, Species)

taxonomy_to_save$accession = NA
write.csv(taxonomy_to_save, "./organism_taxonomy.csv")

