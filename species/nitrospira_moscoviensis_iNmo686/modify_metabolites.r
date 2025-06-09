

source("/projectnb2/talbot-lab-data/zrwerbin/interactions/source.R")
source("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/source.r")
source("/projectnb/dietzelab/zrwerbin/N-cycle/helper_functions.r")

options(scipen=999)

# sybilSBML had to be installed with the help of SCC IT
pacman::p_load(stringr, tidyr, minval, sybilSBML, tidyverse)#, sybilSBML)


in_fp <- paste0("/projectnb2/talbot-lab-data/zrwerbin/interactions/models/misc/lawson_nitrospira.xml")

out_path <- paste0("/projectnb2/talbot-lab-data/metabolic_models/curated_models/iNmo686/iNmo686_modified.xml")


sbml_in <- readSBMLmod(in_fp)

sbml_copy = sbml_in

if (!exists("ref_data")){
	ref_data <- get_reference_data()
}

met_df = cbind.data.frame(orig_met = sbml_copy@met_id) %>% 
	mutate(without_compartment = removeCompartment(orig_met),
				 compartment_no = sbml_copy@met_comp)


compart_key = sbml_copy@mod_compart
names(compart_key) = 1:length(sbml_copy@mod_compart)
met_df$compartment =recode(met_df$compartment_no, !!!compart_key)


recode_df = ref_data$chem_xref %>% 
	filter(source_id %in% met_df$without_compartment) %>% 
	distinct(source_id, .keep_all = T) 

met_df$new_met = recode_df[match(met_df$without_compartment, 
																 recode_df$source_id),]$ID



# Replace any deprecated identifiers
deprecated_recode <- readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/data/deprecated_recode_mets.rds")
met_df$new_met = recode(met_df$new_met, !!!deprecated_recode)


met_df$new_met_name = ifelse(!is.na(met_df$new_met), 
														 met_df$new_met, met_df$without_compartment)
met_df$new_met_out = paste0(met_df$new_met_name, "[", met_df$compartment,"]")

sbml_copy@met_id = met_df$new_met_out

sbml_copy@mod_id = "iNmo686"


# Rename exchange reactions according to new met IDs
# Only for uptake reactions; don't want to mess with anything unnecessarily!
exchReactDF<- findExchReact(sbml_copy)
exchReactDF =	exchReactDF[exchReactDF@uptake]

for (i in 1:length(exchReactDF@react_id)){
	#if (sbml_copy@react_id %in% exchReact@react_id) {
	react_id_index = which(sbml_copy@react_id==exchReactDF@react_id[[i]])
	met_name = removeCompartment(exchReactDF@met_id[[i]])
	new_react_id = paste0("EX_",met_name,"_e")
	print(exchReactDF[i])
	print(new_react_id)
	sbml_copy@react_id[[react_id_index]] = new_react_id
}
sbml_copy@react_id <- gsub("(e)","_e",	sbml_copy@react_id,fixed = T)


writeSBML(sbml_copy, level = 3, filename = out_path, validation = T)


#sbml_copy2 <- readSBMLmod(out_path)
#sybilSBML::writeSBML(sbml_copy2, level = 3, filename = out_path, validation = T)

