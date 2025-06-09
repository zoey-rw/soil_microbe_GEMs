# Trying to map some of the fungal models to metanet namespace
pacman::p_load(tidyverse, minval, SBMLR, rsbml, sybilSBML)
library(BiGGR)
library(janitor)

# Downloaded from metanet using wget
crossref <- data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/interactions/data/chem_xref.tsv", fill=T, skip=351, nThread = 10)
crossref <- crossref %>% separate(col = "#source", into = c("source","source_id"), sep = ":")

in_fp = "/projectnb2/talbot-lab-data/metabolic_models/curated_models/hanpo_z/hanpo_original.xml"
out_fp = "/projectnb2/talbot-lab-data/metabolic_models/curated_models/hanpo_z/hanpo_modified.xml"
hanpo_in_xml = readLines(in_fp)


hanpo_sbml_in <- readSBMLmod(in_fp)

met_original = cbind.data.frame(hanpo_sbml_in@met_id, hanpo_sbml_in@met_name, hanpo_sbml_in@met_attr) %>%
	separate(col = annotation, sep = "metanetx.chemical/", into = c(NA,"metanetx",NA), remove = F) %>%
	separate(col = metanetx, sep = ";", into = c("metanetx", NA)) #%>%
#separate(col = annotation, sep = "chebi/", into = c(NA,"chebi",NA), remove = F) %>%
#separate(col = metanetx, sep = ";", into = c("chebi", NA)) %>%
#filter(!is.na(metanetx))
# %>%


met_to_recode = cbind.data.frame(hanpo_sbml_in@met_id, hanpo_sbml_in@met_name, hanpo_sbml_in@met_attr) %>%
	separate(col = annotation, sep = "metanetx.chemical/", into = c(NA,"metanetx",NA), remove = F) %>%
	separate(col = metanetx, sep = ";", into = c("metanetx", NA)) %>% filter(!is.na(metanetx))

met_to_recode$orig_met = gsub("m[","M_m_",met_to_recode$`hanpo_sbml_in@met_id`, fixed = T)
met_to_recode$orig_met = gsub("]","",met_to_recode$orig_met, fixed = T)
met_to_recode$orig_met = gsub("s[","M_s_",met_to_recode$orig_met, fixed = T)
hanpo_metanet_key = met_to_recode$metanetx
names(hanpo_metanet_key) = met_to_recode$orig_met
#
# below these have annotations for CHEBI, kegg... just found from scrolling through the SBML
# M_s_0421 = ammonium
# M_s_0420 = ammonium
# M_s_0419 = ammonium
# M_s_0794 = H+
# 	M_m_0551 = formaldehyde =	MNXM736417
# M_m_0554= methanol = MNXM729800
# M_s_0025=R lactate=MNXM731834
# M_s_0026=R lactate=MNXM731834
# M_s_0069 = M_s_0068 = M_s_0066 = S malate = MNXM1107192
# M_s_0021 = M_s_0022 = M_s_0023 = M_s_0024 = R carnitine = MNXM1105737

	hanpo_metanet_key <- c(hanpo_metanet_key,
												 "M_s_0419" = "M_s_MNXM15",
												 "M_s_0420" = "M_s_MNXM15",
												 "M_s_0421" = "M_s_MNXM15",
												 "M_s_0794" = "M_s_MNXM15",
												 "M_m_0551" = "M_m_MNXM736417",
												 "M_m_0554" = "M_m_MNXM729800",
												 "M_s_0025" = "M_s_MNXM731834",
												 "M_s_0026" = "M_s_MNXM731834",
												 "M_s_0027" = "M_s_MNXM731834",
												 "M_s_0066" = "M_s_MNXM1107192",
												 "M_s_0068" = "M_s_MNXM1107192",
												 "M_s_0069" = "M_s_MNXM1107192",
												 "M_s_0021" = "M_s_MNXM1105737",
												 "M_s_0022" = "M_s_MNXM1105737",
												 "M_s_0023" = "M_s_MNXM1105737",
												 "M_s_0024" = "M_s_MNXM1105737")

# met_original$metanetx = ifelse(!is.na(met_original$metanetx),
# 															 met_original$metanetx, met_original$`hanpo_sbml_in@met_id`)
# hanpo_metanet_key = met_original$metanetx
# names(hanpo_metanet_key) = met_original$`hanpo_sbml_in@met_id`

# separate(col = annotation, sep = ";", into = c("url1","url2","url3","url4","url5"))

hanpo_sbml_copy = hanpo_sbml_in
#hanpo_sbml_copy@met_id = recode(hanpo_sbml_copy@met_id, !!!hanpo_metanet_key)
#writeSBML(hanpo_sbml_copy,filename = out_fp)


hanpo_xml_copy = hanpo_in_xml
for (i in 1:length(hanpo_metanet_key)){
	hanpo_xml_copy <- gsub(names(hanpo_metanet_key)[[i]], hanpo_metanet_key[[i]], hanpo_xml_copy,fixed = T)
}

writeLines(text = hanpo_xml_copy, con = out_fp, sep = "\n")


