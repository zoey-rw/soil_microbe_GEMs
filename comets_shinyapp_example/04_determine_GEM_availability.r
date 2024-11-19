# First run 01_merge_env_data.r
# First run 02_determine_biome_preference.r
# First run 03_determine_pH_temp_preference.r

library(data.table)
library(tidyverse)
library(DT)

set.seed(1)

# Read in list of manually-curated models
soil_GEM_info = read_csv("./Soil GEMs - manually curated.csv") %>% 
    dplyr::select(GEM_ID, Kingdom, Genus, Species.strain = `Species/strain`, Citation, 
           Source=`Source/link`, Filepath = `SCC Filepath`) %>% 
    mutate(Species = paste(Genus, word(Species.strain)))

# Read in and merge other associated information
organism_info = fread("./intermediate_data/organism_info.csv")
organism_pH_temp_preference = fread("./intermediate_data/organism_pH_temp_preference.csv")
organism_biome_preference = fread("./intermediate_data/organism_biome_preference.csv") %>% mutate(`Species of interest` = taxon)

organism_preference = full_join(organism_pH_temp_preference, organism_biome_preference)
organism_preference = left_join(organism_preference, organism_info)

# Create list to loop through
taxon_names = organism_preference$`Species of interest`

# Match by names
gem_out = list()
for (i in 1:length(taxon_names)){
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
    }
    match_out <- match_out %>% dplyr::select(-Species)
    match_out$`Species of interest` = spec
    gem_out[[i]] = match_out
}

gem_out_table = data.table::rbindlist(gem_out) %>% dplyr::select(`Species of interest`, match_by,
                                                                 GEM_ID#, #"Genome source" = source,
                                                                 #"Genome (link)" = genome_link
                                                                 ) %>% 
    distinct(.keep_all = T)


organism_df <- left_join(organism_preference, gem_out_table) 

df_to_subset <- organism_df %>% 
    dplyr::select(taxonomy_id, temperature_preference, pH_preference,
           `Species of interest`, taxon, 
           cultivatedCrops, deciduousForest, dwarfScrub, 
           emergentHerbaceousWetlands,
           evergreenForest, grasslandHerbaceous, mixedForest, pastureHay, 
           sedgeHerbaceous, shrubScrub, woodyWetlands,
           "Metagenome-assembled genome?" = is_MAG,
           "GEM match criteria" = match_by,
           GEM_ID, 
           "Genome source" = source) %>% 
    distinct(.keep_all = T)

# Add in the genome accessions from database creation!
soil_microbe_db_struo = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/soil_genome_db_struo.tsv")
refsoil_struo = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/refsoil_struo.tsv")
struo_organism_data = rbindlist(list(soil_microbe_db_struo, refsoil_struo), fill = T) %>% rename(taxonomy_id=ncbi_species_taxid) %>% dplyr::select(taxonomy_id, accession) %>% distinct()

df_to_subset <- left_join(df_to_subset, struo_organism_data, relationship = "many-to-many") 


write.csv(df_to_subset, "./intermediate_data/organism_data_to_subset.csv")




#lineage_df =  split_taxonomy_ncbi(bracken_with_lineage$lineage)
lineage_df =  lapply(organism_info$lineage, 
                     function(x) {
                         phyloseq::parse_taxonomy_qiime(x) %>% #as.data.frame() %>% 
                             t()  %>% as.data.frame()
                     }) %>% 
    data.table::rbindlist()
lineage_df <- cbind.data.frame(lineage_df, "taxonomy_id" = organism_info$taxonomy_id)

taxonomy_to_save = lineage_df %>% select(taxonomy_id, taxon=Species, #accession, 
                                         Kingdom, Phylum, Class, Order, Family, Genus, Species)

write.csv(taxonomy_to_save, "./intermediate_data/organism_taxonomy.csv")


# 
# organism_info[grepl("heiligendammensis", organism_info$ncbi_organism_name),]
# sample_abundance[sample_abundance$taxonomy_id=="92490",]
# 
# pigz -cd taxid-changelog.csv.gz \
# | csvtk grep -f taxid -p 92490 \
# | csvtk cut -f -lineage,-lineage-taxids \
# | csvtk pretty 
# 
# # Other struo files
# spire_struo = read_tsv("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/spire_MAGs_struo.tsv") %>% 
#     mutate(is_MAG=T)
# SMAG_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/SMAG_struo.tsv") %>% 
#     mutate(is_MAG=T)
# nayfach_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/nayfach_MAGs_struo.tsv") %>% 
#     mutate(source="GEM catalog") %>% mutate(is_MAG=T)
# JGI_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/ncbi_struo.tsv")
# myco_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/mycocosm_published_struo.tsv")
