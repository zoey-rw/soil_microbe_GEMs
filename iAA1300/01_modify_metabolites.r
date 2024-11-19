# Tried to use the annotations this model (iAA1300) already had, but many seemed wrong, 
# e.g. same identifiers for O2 and SO2, or H and NADPH

# sybilSBML had to be installed with the help of SCC IT
pacman::p_load(here, stringr, tidyr, minval, sybilSBML, tidyverse)

source(here("source_functions.R"))

if (!exists("ref_data")) ref_data <- get_reference_data(chem_xref_path = "/projectnb2/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/chem_xref.tsv", 
																												chem_prop_path = "/projectnb2/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/chem_prop.tsv", 
																												reac_prop_path = "/projectnb2/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/reac_prop.tsv", 
																												reac_xref_path = "/projectnb2/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/reac_xref.tsv")

# Read in list of deprecated IDs
deprecated_recode <- readRDS(here("reference_data","deprecated_recode_mets.rds"))
																	
in_fp <- here("iAA1300","iAA1300_original.xml")
out_dir <- here("iAA1300")
spec <- "iAA1300"
out_path <- file.path(out_dir, paste0(spec,"_modified.xml"))

	
	cat("\nRepairing model", spec, "\n")
	out_path <- file.path(out_dir, paste0(spec,"_modified.xml"))
	
	if (file.exists(out_path)){
		cat("\nSkipping model ", spec, "\n")
	} else {
		# Relies on sybilSBML
		sbml_in <- readSBMLmod(in_fp)
		
		sbml_copy = sbml_in
		
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
			met_df[met_df$new_met %in% deprecated$`#deprecated_ID`,]
			met_df$new_met = recode(met_df$new_met, !!!deprecated_recode)
			met_df[met_df$new_met %in% deprecated$`#deprecated_ID`,]
			
#			met_df$new_met = recode(met_df$new_met, "MNXM735623" = "MNXM1137670")
			
			
		met_df$new_met_name = ifelse(!is.na(met_df$new_met), 
																 met_df$new_met, met_df$without_compartment)
		met_df$new_met_out = paste0(met_df$new_met_name, "[", met_df$compartment,"]")
		
		met_df$old_met_out = paste0(met_df$without_compartment, "[", met_df$compartment,"]")
		
		# In case the new ids are duplicated - revert to old ones
		duped_mets = met_df[duplicated(met_df$new_met_out),]$new_met_name
		met_df$new_met_out = ifelse(met_df$new_met_name %in% duped_mets,met_df$old_met_out,met_df$new_met_out)
		
		
		sbml_copy@met_id = met_df$new_met_out

		
		
		# Rename exchange reactions according to new met IDs
		# Only for uptake reactions; don't want to mess with anything unnecessarily!
		exchReactDF<- findExchReact(sbml_copy)
		exchReactDF =	exchReactDF[exchReactDF@uptake]
		
		for (i in 1:length(exchReactDF@react_id)){
			react_id_index = which(sbml_copy@react_id==exchReactDF@react_id[[i]])
			met_name = removeCompartment(exchReactDF@met_id[[i]])
			new_react_id = paste0("EX_",met_name,"_e")
			print(exchReactDF[i])
			print(new_react_id)
			sbml_copy@react_id[[react_id_index]] = new_react_id
		}
		
		
		sbml_copy@react_id <- gsub("(e)","_e",	sbml_copy@react_id,fixed = T)
		
		writeSBML(sbml_copy, level = 3, filename = out_path)
		
	}
}
