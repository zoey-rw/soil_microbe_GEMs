
# sybilSBML had to be installed with the help of SCC IT
pacman::p_load(stringr, tidyr, minval, sybilSBML, tidyverse)

source(here("source_functions.R"))

if (!exists("ref_data")) ref_data <- get_reference_data()

# Read in list of deprecated IDs
deprecated_recode <- readRDS(here("reference_data","deprecated_recode_mets.rds"))
																	
in_fp <- here("iFC579","iFC579_original.xml")
out_dir <- here("iFC579")
	spec <- "iFC579"
	out_path <- file.path(out_dir, paste0(spec,"_modified.xml"))
	
	cat("\nRepairing model", spec, "\n")
	
		# Relies on sybilSBML
		sbml_in <- readSBMLmod(in_fp)
		
		sbml_copy = sbml_in
		
		met_df = cbind.data.frame(orig_met = sbml_copy@met_id) %>% 
			mutate(without_compartment = removeCompartment(orig_met),
						 compartment_no = sbml_copy@met_comp)
		
		met_df[is.na(sbml_copy@met_comp),]$compartment_no = 1
		
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
			
		met_df$new_met_name = ifelse(!is.na(met_df$new_met), 
																 met_df$new_met, met_df$without_compartment)
		met_df$new_met_out = paste0(met_df$new_met_name, "[", met_df$compartment,"]")
		
		sbml_copy@met_id = met_df$new_met_out
		
		
		
		
		
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
		

		
		writeSBML(sbml_copy, level = 3, filename = out_path)
		
		sbml_copy2 <- readSBMLmod(out_path)
		sybilSBML::writeSBML(sbml_copy2, level = 3, filename = out_path, validation = T)
		
