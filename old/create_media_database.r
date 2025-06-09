
# sybilSBML had to be installed with the help of SCC IT
pacman::p_load(stringr, tidyr, minval, sybilSBML, tidyverse)
source("/projectnb2/talbot-lab-data/zrwerbin/interactions/source.R")
source("/projectnb/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/02_check_equation_stoichiometry.r")
options(scipen=999)


if (!exists("ref_data")) ref_data <- get_reference_data(reac_prop_path = "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/reac_prop.tsv",
																												chem_prop_path = "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/chem_prop.tsv",
																												#comp_xref_path = ,
																												chem_xref_path = "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/chem_xref.tsv",
																												reac_xref_path = "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/reac_xref.tsv"
																													)



media_tsv = read_tsv("https://raw.githubusercontent.com/franciscozorrilla/SymbNET/main/data/media_db.tsv")


head(media_tsv)




recode_df = ref_data$chem_xref %>% 
	filter(source_id %in% c(media_tsv$compound,media_tsv$name)) %>% 
	distinct(source_id, .keep_all = T) %>% filter(!is.na(source_id))

media_tsv$MNX_ID = recode_df[match(media_tsv$compound, 
																			recode_df$source_id),]$ID
media_tsv$MNX_desc = recode_df[match(media_tsv$compound, 
																	 recode_df$source_id),]$description

# Replace any deprecated identifiers
deprecated_recode <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/deprecated_recode_mets.rds")
media_tsv$updated_MNX_ID = recode(media_tsv$MNX_ID, !!!deprecated_recode)

media_tsv[,c("MNX_ID","updated_MNX_ID")] %>% unique

write_tsv(media_tsv, "/projectnb2/talbot-lab-data/zrwerbin/soil_microbe_GEMs/reference_data/media_database.tsv")



major_carbon_sources =c("algac__M_e",
	"arab__D_e",
	"cellb_e",
	"fru_e",
	"gal_e",
	"glc__D_e",
	"lcts_e",
	"lyx__L_e",
	"malt_e",
	"man_e",
	"rib__D_e",
	"sucr_e",
	"tre_e",
	"xyl__D_e",
	"acgam_e",
	"chitin_e",
	"chitos_e",
	"pectin_e",
	"starch_e",
	"raffin_e",
	"xylan4_e",
	"succ_e",
	"pyr_e")

minor_carbon_sources = c("algac__M_e",
	"arab__D_e",
	"cellb_e",
	"fru_e",
	"gal_e",
	"glc__D_e",
	"lcts_e",
	"lyx__L_e",
	"malt_e",
	"man_e",
	"rib__D_e",
	"sucr_e",
	"tre_e",
	"xyl__D_e",
	"acgam_e",
	"chitin_e",
	"chitos_e",
	"pectin_e",
	"starch_e",
	"raffin_e",
	"xylan4_e",
	"succ_e",
	"pyr_e",
	"arg__L_e",
	"his__L_e",
	"lys__L_e",
	"asp__L_e",
	"glu__L_e",
	"ser__L_e",
	"thr__L_e",
	"asn__L_e",
	"gln__L_e",
	"cys__L_e",
	"gly_e",
	"pro__L_e",
	"ala__L_e",
	"val__L_e",
	"ile__L_e",
	"leu__L_e",
	"met__L_e",
	"phe__L_e",
	"tyr__L_e",
	"trp__L_e")


major_carbon_mets = media_tsv %>% 
	mutate(bigg_metabolite = paste0(compound, "_e"),
				 MNX_metabolite = paste0(updated_MNX_ID, "_e")) %>% 
	filter(bigg_metabolite %in% major_carbon_sources) %>% 
	select(MNX_metabolite) %>% unique() %>% unlist


minor_carbon_mets = media_tsv %>% 
	mutate(bigg_metabolite = paste0(compound, "_e"),
				 MNX_metabolite = paste0(updated_MNX_ID, "_e")) %>% 
	filter(bigg_metabolite %in% minor_carbon_sources) %>% 
	select(MNX_metabolite) %>% unique() %>% unlist




m8_media = media_tsv %>% 
	mutate(MNX_metabolite = paste0(updated_MNX_ID, "_e")) %>% 
	filter(medium == "M8") %>% 
	select(MNX_metabolite) %>% unique() %>% unlist

cat(paste0('"',major_carbon_mets,'"', ":0,"))
cat(paste0('"',minor_carbon_mets,'"', ":0,"))
cat(paste0('"',m8_media,'"', ":1000,"))




