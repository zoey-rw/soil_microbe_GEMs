#install.packages("here")
library(here)
here()
#library(minval)

options(scipen=1000)

##### fix_metabolite
fix_metabolite <- function(x){
    fine <- ifelse(grepl("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", x), T, F)
    if(fine){
        x <- str_replace_all(x,"\\[c\\]", "_c")
        x <- str_replace_all(x,"\\[p\\]", "_p")
        x <- str_replace_all(x,"\\[e\\]", "_e")
        return(x)
    } else {
        # Keep everything before the second underscore
        #	simple <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', simple), ' ')[[1]][[1]]
        simple <- strsplit(x, "_")[[1]][[1]]
        simple <- gsub("__45_|__91_|__93_", "_", x)
        compartment_letter <- ifelse(grepl("__c__", x), "c",
                                     ifelse(grepl("__e__", x), "e",
                                            ifelse(grepl("__p__", x), "p", "")))
        simple <- strsplit(simple, "__[cep]")[[1]][[1]]
        simple_compartment <- paste(simple, compartment_letter, sep = "_")
        return(simple_compartment)
    }
}

##### fix_reaction_compartment
fix_reaction_compartment <- function(x){
    x <- str_replace_all(x,"\\(c\\)", "_c")
    x <- str_replace_all(x,"\\(p\\)", "_p")
    x <- str_replace_all(x,"\\(e\\)", "_e")
    return(x)
}


##### removeCompartment from metabolite
removeCompartment <- function (metabolite, rmCoef = FALSE) {
    metabolite <- minval:::removeSpaces(metabolite = metabolite)
    if (rmCoef == TRUE) {
        metabolite <- removeCoefficients(metabolite)
    }
    metabolite <- gsub("\\_[[:alnum:]]$|\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$",
                       "", metabolite)
    return(metabolite)
}

##### metabolites from minval
metabolites <- function (reactionList, woCompartment = FALSE, uniques = TRUE, metanet_replace = NULL)
{
    reaction <- strsplit(as.vector(reactionList), "[[:blank:]]+<?=>[[:blank:]]*")
    reaction <- lapply(reaction, function(reaction) {
        strsplit(unlist(reaction), "[[:blank:]]+\\+[[:blank:]]+")
    })
    reaction <- lapply(reaction, function(reaction) {
        minval:::removeSpaces(unlist(reaction))
    })
    reaction <- lapply(reaction, function(reaction) {
        minval:::removeCoefficients(reaction)
    })
    metabolites <- unlist(reaction)
    if (woCompartment == TRUE) {
        metabolites <- removeCompartment(metabolites)
    }
    if (uniques == TRUE) {
        metabolites <- unique(metabolites)
    }
    
    # if (!is.null(metanet_replace)) {
    #   if (metabolites %in% c(""))
    #   metanet_id = metanet_replace[grepl(metabolites, metanet_replace$Metabolite.name),]$metanet_id %>% unique()
    #   metabolites = metanet_id
    # }
    return(metabolites)
}


##### get coefficients from minval
coefficients <- function (metabolite) {
    metabolite <- regmatches(x = metabolite, m = gregexpr(pattern = "^[[:digit:]]{0,5}[[:punct:]]*[[:digit:]]*[[:blank:]]+",
                                                          text = metabolite))
    metabolite[lengths(metabolite) == 0] <- 1
    metabolite <- as.numeric(minval:::removeSpaces(metabolite))
    return(metabolite)
}


##### removeCoefficients from minval
removeCoefficients <- function (metabolite) {
    metabolite <- gsub(pattern = "^[[:digit:]]{0,5}[[:punct:]]*[[:digit:]]*[[:blank:]]+",
                       replacement = "", x = metabolite)
    return(metabolite)
}


##### compartments from minval 
compartments <- function (reactionList, uniques = TRUE)
{
    if (uniques == TRUE) {
        metabolites <- metabolites(reactionList = reactionList,
                                   uniques = TRUE)
    }
    else {
        metabolites <- metabolites(reactionList = reactionList,
                                   uniques = FALSE)
    }
    compartments <- unlist(regmatches(x = metabolites, m = gregexpr(pattern = "\\_[[:alnum:]]$|\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$",
                                                                    text = metabolites)))
    compartments <- gsub("\\_|\\[|\\]", "", compartments)
    if (length(compartments) == 0) {
        compartments <- ifelse(grepl("__c__", metabolites), "c",
                               ifelse(grepl("__e__", metabolites), "e",
                                      ifelse(grepl("__p__", metabolites), "p",
                                             ifelse(grepl("__91__c__93__",metabolites), "c",
                                                    ifelse(grepl("__91__e__93__",metabolites), "e",
                                                           ifelse(grepl("__91__p__93__",metabolites), "p", ""))))))
    }
    if (length(compartments) == 0) {
        compartments <- NA
    }
    if (uniques == TRUE) {
        compartments <- unique(compartments)
        return(compartments)
    }
    else {
        return(compartments)
    }
}

##### extractData from minval
extractData <- function (inputData, boundary = "b") {
    exchange <- minval:::reactionType(inputData[["REACTION"]]) == "Exchange reaction"
    if (any(exchange) == TRUE) {
        # inputData[["REACTION"]][exchange] <- as.vector(sapply(metabolites(inputData[["REACTION"]][exchange]),
        # 																											function(x) {
        # 																												paste0(x, " <=> ",
        # 																															 paste0(metabolites(reactionList = x, woCompartment = TRUE),
        # 																															 			 "[", boundary, "]"))
        # 																											}))
    }
    data <- list()
    data$COMPARTMENTS <- compartments(inputData[["REACTION"]])
    data$METABOLITES <- metabolites(inputData[["REACTION"]],
                                    uniques = TRUE)
    data$REACTIONS <- lapply(seq_along(inputData[["REACTION"]]),
                             function(reaction) {
                                 list(id = as.vector(inputData[["ID"]])[reaction],
                                      reversible = ifelse(test = grepl(pattern = "<=>",
                                                                       x = inputData[["REACTION"]][reaction]), yes = "true",
                                                          no = "false"), gpr = as.vector(inputData[["GPR"]])[reaction],
                                      reactants = unlist(minval:::getLeft(inputData[["REACTION"]][reaction])),
                                      products = unlist(minval:::getRight(inputData[["REACTION"]][reaction])),
                                      lowbnd = ifelse(test = is.numeric(inputData[["LOWER.BOUND"]][reaction]),
                                                      yes = inputData[["LOWER.BOUND"]][reaction],
                                                      no = -1000), upbnd = ifelse(test = inputData[["UPPER.BOUND"]][reaction] !=
                                                                                      "", yes = inputData[["UPPER.BOUND"]][reaction],
                                                                                  no = 1000), objective = ifelse(test = inputData[["OBJECTIVE"]][reaction] !=
                                                                                                                     "", yes = inputData[["OBJECTIVE"]][reaction],
                                                                                                                 no = 0))
                             })
    return(data)
}

##### rearmReactions from minval
rearmReactions <- minval:::rearmReactions

##### repair_SBML_model 
repair_SBML_model <- function (modelData,
                               modelID = "model",
                               boundary = "e") {
    if (class(modelData) == "data.frame") {
        modelData <- minval:::validateData(modelData = modelData)
        modelData <- minval:::removeComments(modelData = modelData)
        modelData <- modelData[minval:::validateSyntax(modelData[["REACTION"]]),
        ]
    } else if (class(modelData) == "modelorg") {
        if(length(modelData@gpr)==0) {
            modelData@gpr <- c(rep("", length(modelData@react_rev))) }
        modelData <- minval:::convertData(model = modelData)
    } else {
        stop("Input format not supported.")
    }
    modelData <- extractData(inputData = modelData, boundary = "e")
    header <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
                "<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">",
                paste0("\t<model id=\"", modelID,
                       "\" name=\"", modelID, "\">"), "\t\t<notes>", "\t\t\t<body xmlns=\"http://www.w3.org/1999/xhtml\">",
                "\t\t\t<p> Generated with MINVAL: an R package for MINimal VALidation of stoichiometric reactions </p>",
                "\t\t\t</body>", "\t\t</notes>")
    comp <- "\t\t<listOfCompartments>"
    comp <- c(comp, as.vector(sapply(modelData[["COMPARTMENTS"]],
                                     function(compartment) {
                                         paste0("\t\t\t<compartment id=\"", compartment, "\" name=\"",
                                                compartment, "\"/>")
                                     })))
    comp <- c(comp, "\t\t</listOfCompartments>")
    mets <- "\t\t<listOfSpecies>"
    mets <- c(mets, sapply(modelData[["METABOLITES"]], function(metabolite) {
        metabolite = fix_metabolite(metabolite)
        paste0("\t\t\t<species id=\"M_", metabolite, "\" name=\"",
               metabolites(metabolite, woCompartment = TRUE, metanet_replace = metanet_replace), "\" compartment=\"",
               compartments(metabolite), "\" boundaryCondition=\"",
               ifelse(test = compartments(metabolite) == boundary,
                      yes = "true", no = "false"), "\"/>")
    }))
    mets <- c(mets, "\t\t</listOfSpecies>")
    react <- "\t\t<listOfReactions>"
    react <- c(react, unlist(sapply(seq_along(modelData[["REACTIONS"]]), function(reaction) {
        modelData[["REACTIONS"]][[reaction]][["id"]] <- janitor::make_clean_names(modelData[["REACTIONS"]][[reaction]][["id"]], case= "none")
        
        c(paste0("\t\t\t<reaction id=\"R_", modelData[["REACTIONS"]][[reaction]][["id"]],
                 "\"  reversible=\"", modelData[["REACTIONS"]][[reaction]][["reversible"]],
                 "\">"),
          paste0("\t\t\t\t<notes>"),
          paste0("\t\t\t\t\t<html xmlns=\"http://www.w3.org/1999/xhtml\">",
                 if (isTRUE(modelData[["REACTIONS"]][[reaction]][["gpr"]] !=
                            "")) {
                     paste0("<p>GENE_ASSOCIATION: ", modelData[["REACTIONS"]][[reaction]][["gpr"]],
                            "</p>")
                 }, "</html>"),
          paste0("\t\t\t\t</notes>"),
          paste0("\t\t\t\t<listOfReactants>"),
          unname(sapply(X = modelData[["REACTIONS"]][[reaction]][["reactants"]], function(x) {
              metabolite <- fix_metabolite(x)
              paste0("\t\t\t\t\t<speciesReference species=\"M_",
                     metabolites(metabolite), "\" stoichiometry=\"",
                     coefficients(metabolite), "\"/>")
          })),
          paste0("\t\t\t\t</listOfReactants>"),
          if(!is.na(modelData[["REACTIONS"]][[reaction]][["products"]][[1]])) {
              c("\t\t\t\t<listOfProducts>",
                sapply(modelData[["REACTIONS"]][[reaction]][["products"]],
                       function(x) {
                           metabolite <- fix_metabolite(x)
                           c(paste0("\t\t\t\t\t<speciesReference species=\"M_",
                                    metabolites(metabolite), "\" stoichiometry=\"",
                                    coefficients(metabolite), "\"/>"))
                       }),
                "\t\t\t\t</listOfProducts>")
              #} else c("remove"))
          },
          #),
          #																			paste0("\t\t\t\t</listOfProducts>"),
          paste0("\t\t\t\t<kineticLaw>"),
          paste0("\t\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">"),
          paste0("\t\t\t\t\t\t<ci>FLUX_VALUE</ci>"), paste0("\t\t\t\t\t</math>"),
          paste0("\t\t\t\t\t<listOfParameters>"), paste0("\t\t\t\t\t\t<parameter id=\"LOWER_BOUND\" value=\"",
                                                         modelData[["REACTIONS"]][[reaction]][["lowbnd"]],
                                                         "\"/>"), paste0("\t\t\t\t\t\t<parameter id=\"UPPER_BOUND\" value=\"",
                                                                         modelData[["REACTIONS"]][[reaction]][["upbnd"]],
                                                                         "\"/>"), paste0("\t\t\t\t\t\t<parameter id=\"OBJECTIVE_COEFFICIENT\" value=\"",
                                                                                         modelData[["REACTIONS"]][[reaction]][["objective"]],
                                                                                         "\"/>"), paste0("\t\t\t\t\t\t<parameter id=\"FLUX_VALUE\" value=\"0\"/>"),
          paste0("\t\t\t\t\t</listOfParameters>"), paste0("\t\t\t\t</kineticLaw>"),
          paste0("\t\t\t</reaction>"))
    })))
    #react <- gsub(pattern = "\t\t</listOfProducts>\t\t</listOfProducts>", "", react)
    react <- c(react, "\t\t</listOfReactions>")
    end <- c("\t</model>", "</sbml>")
    model <- as.vector(c(header, comp, mets, react, end))
    return(model)
}


#' Functions to aid in balancing equations within genome-scale metabolic models
#' To be used in conjunction with the Python package memote "find_charge_unbalanced_reactions()"
#' 
#' TO DO: output from "get_bigg_charges" could be presented in a cleaner/more useful way
#' 
#' possibly useful: https://github.com/SBRG/bigg_models/blob/cf16b53ab77ea699a1e132ce10a3eca1690e0aee/bin/load_metanetx


#' DOWNLOAD REFERENCE DATA
#####
#' 
#' Get reference data on names and charges for chemical reactions, metabolites, and cellular compartments from MetanetX and BiGG
#'  Files can be downloaded from Metanetx/BiGG (e.g. using wget) or called directly from site via URL
#' @param reac_xref_path 
#' @param chem_xref_path
#' @param comp_xref_path
#' @param reac_prop_path
#' @param chem_prop_path
#' @param bigg_met_path 
#' @param bigg_rxn_path 
#'
#' @return returns a list with 7 named dataframes for reference
#' @export
#'
#' @examples
#' ref_data <- get_reference_data()
get_reference_data <- function(threads =10, 
                               reac_xref_path = here("reference_data", "reac_xref.tsv"), 
                               chem_xref_path = here("reference_data", "chem_xref.tsv"),
                               comp_xref_path = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv",
                               reac_prop_path = here("reference_data", "reac_prop.tsv"), 
                               chem_prop_path = here("reference_data", "chem_prop.tsv"),
                               bigg_met_path = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
                               bigg_rxn_path = "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt") {
    
    # Cross-referencing tables
    reac_xref <- data.table::fread(reac_xref_path, fill=T, skip=351, nThread = threads) %>% 
        separate(col = "#source", into = c("source","source_id"), sep = ":")
    chem_xref <- data.table::fread(chem_xref_path, fill=T, skip=351, nThread = threads) %>% 
        separate(col = "#source", into = c("source","source_id"), sep = ":")
    comp_xref <- data.table::fread(comp_xref_path, fill=T, skip=351, nThread = threads) %>% 
        separate(col = "#source", into = c("source","source_id"), sep = ":")
    
    # Property tables
    reac_prop <- data.table::fread(reac_prop_path, fill=T, skip=351, nThread = threads)
    chem_prop <- data.table::fread(chem_prop_path, fill=T, skip=351, nThread = threads)
    
    # BiGG data from url
    bigg_met <- data.table::fread(bigg_met_path, fill=T)
    bigg_rxn <- data.table::fread(bigg_rxn_path, fill=T) %>% 
        separate(reaction_string, sep = " <-> ", into = c("left","right"), remove = F,fill = "right")
    
    out <- list("chem_xref" = chem_xref,
                "comp_xref" = comp_xref,
                "reac_xref" = reac_xref,
                "reac_prop" = reac_prop,
                "chem_prop" = chem_prop,
                "bigg_met" = bigg_met,
                "bigg_rxn" = bigg_rxn) 
    return(out)
}


##### get BiGG metabolite ID from metanet ID
# TO DO: could be vectorized to be faster
get_bigg_met <- function(metanetxID = "MNXM3453",
                         reference_data = NULL) {
    chem_xref = reference_data$chem_xref
    bigg_met = reference_data$bigg_met
    
    chem_met_info <- chem_xref %>% filter(ID == !!metanetxID & source == "biggM" & !grepl("M_", source_id)) %>% select(source_id)
    bigg_met_ID <- unique(chem_met_info$source_id)
    return(bigg_met_ID)
    
}

#get_bigg_met(metanetxID = "MNXM3453", reference_data = ref_data)

#####
#' INVESTIGATE SPECIFIC REACTIONS
#' 
#' Get data on names and charges for specific chemical reaction from MetanetX and BiGG
#' @param metanetxID Optional - ID of problematic/unbalanced chemical reaction
#' @param biggID Optional - ID of problematic/unbalanced chemical reaction
#' @param reference_data
#'
#' @return returns series of messages describing the metabolite(s) missing charges in a reaction, and the typical charges for that metabolite in BiGG models
#' @export
#'
#' @examples
#' get_bigg_charges(metanetxID = "MNXR94775",
#' 								 reference_data = ref_data) 
#' 
#' get_bigg_charges(metanetxID = NULL, 
#' 								 biggID = "2DDARAA",
#' 								 reference_data = ref_data)
get_bigg_charges <- function(metanetxID = NULL, 
                             biggID = NULL,
                             reference_data = NULL) {
    
    pacman::p_load(tidyverse, minval, curl) 
    
    if (is.null(reference_data)) {
        message("First run get_reference_data() and pass the output as the argument to reference_data")
    } else {
        chem_xref = reference_data$chem_xref
        comp_xref = reference_data$comp_xref
        reac_xref = reference_data$reac_xref
        reac_prop = reference_data$reac_prop
        chem_prop = reference_data$chem_prop
        bigg_met = reference_data$bigg_met
        bigg_rxn = reference_data$bigg_rxn
    }
    
    # Get metanetID from BiGG ID, if necessary
    if (!is.null(biggID)) {
        rxn_info <- reac_xref %>% filter(source_id == !!biggID) 
        metanetxID <- unique(rxn_info$ID)
    }
    
    rxn <- reac_prop %>% filter(`#ID` == !!metanetxID)
    
    rxn_split <- strsplit(rxn$mnx_equation, " = ") %>% unlist()
    
    #side <- "left"
    for (side in c("left","right")) {
        rxn_split <- strsplit(rxn$mnx_equation, " = ") %>% unlist()
        rxn_side <- switch(side,
                           "left" = rxn_split[[1]],
                           "right" = rxn_split[[2]])
        
        met <- minval::metabolites(rxn_side) %>% 
            lapply(., function(x) gsub("@MNX[CD][0-9]","", x)) %>% 
            unlist()
        
        compartments <- minval::metabolites(rxn_side) %>% 
            sapply(., function(x) {strsplit(x, "@") %>% 
                    sapply(tail, 1 )})
        names(compartments) <- met
        
        init_info <- list()
        
        for (met in met) {
            db_entry <-  chem_prop %>% filter(`#ID` == met)
            init_info[[met]] <- list("charge" = db_entry$charge,
                                     "mass" = db_entry$mass)
        }
        init_charge <- lapply(init_info, "[[", 1)
        init_mass <- lapply(init_info, "[[",2)
        
        missing_charges <- init_charge[is.na(init_charge)]
        bigg_charges <- init_charge
        names(bigg_charges) <- lapply(names(bigg_charges), function(x) get_bigg_met(metanetxID = x, reference_data = reference_data)) %>% unlist()
        
        if (length(missing_charges)==0) {
            message("MetanetX has the following charges on ", side, " of equation: ")
            print(unlist(bigg_charges))
        }
        
        #met <- "MNXM3453" # for testing
        for (met in names(missing_charges)){
            bigg_met_id <- get_bigg_met(metanetxID = met, reference_data = reference_data)
            
            message("MetanetX is missing charge on ", side, " side of equation for: \nMetanetX metabolite: ", met, "\nBiGG metabolite: ", bigg_met_id)
            
            
            # Get list of BiGG reactions that involve this metabolite on the same side
            bigg_ids <- chem_xref %>% filter(ID == met) %>% select(source_id) %>% unlist()
            bigg_universal_id <- bigg_met %>% filter(universal_bigg_id %in% bigg_ids) %>% select(universal_bigg_id) %>% unique() %>% unlist()
            if (side == "left") {
                rxn_list <- bigg_rxn %>% filter(grepl(bigg_universal_id, left)) %>% select(bigg_id) %>% unlist()
            } else {
                rxn_list <- bigg_rxn %>% filter(grepl(bigg_universal_id, right)) %>% select(bigg_id) %>% unlist()
            }
            
            #rxn_of_interest <- "2AGPGAT161"
            rxn_of_interest <- "2AGPG161tipp"
            for (rxn_of_interest in rxn_list) {
                
                message("Metabolite used in reaction: ", rxn_of_interest)
                rxn_info <- bigg_rxn %>% filter(bigg_id == !!rxn_of_interest) 
                
                model_list <- rxn_info$model_list %>% strsplit("; ") %>% unlist()
                api_url <- paste0("http://bigg.ucsd.edu/api/v2/models/", model_list[[1]], "/reactions/", rxn_of_interest)
                
                #if (RCurl::url.exists(api_url)) {
                json_out <- jsonlite::fromJSON(api_url, flatten = T)
                charges <- json_out$metabolites %>% filter(bigg_id == !!bigg_universal_id) 
                
                # Get metabolite compartment for decision-making
                comp <- compartments[met]
                comp_desc <- comp_xref[comp_xref$ID==comp,]$description
                
                
                # Report output
                for (k in 1:nrow(charges)){
                    message("Observed charge for: ", met, " is ", charges[k,]$stoichiometry, " in ", "compartment: ", charges[k,]$compartment_bigg_id)
                }
            }
            # After the final loop, report the actual compartment to decide which of the charges is most appropriate
            message("Metabolite ", met, " is in compartment: ", comp, ", which has the following description from MetanetX: ", comp_desc)
        }
    }
}

#ref_data <- get_reference_data()

##### Parsing NEON soil sample IDs to pull out timing/location metadata
parseNEONsampleIDs <- function(sampleID){
    df <- data.frame(siteID = substr(sampleID, 1, 4), sampleID = sampleID, stringsAsFactors = F) %>% 
        mutate(sample = sapply(strsplit(sampleID, "-GEN|-gen"),  "[[" , 1)) %>% 
        mutate(geneticSampleID = sapply(strsplit(sampleID, "-DNA"),  "[[" , 1)) %>% 
        mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>% 
        mutate(dates = sapply(strsplit(sample, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>% 
        mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>% 
        mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>% 
        mutate(dateID = substr(as.character(dates), 1, 6)) %>% 
        mutate(plotID = substr(sample, 1, 8)) %>% 
        mutate(site_date = paste0(siteID, "-", dateID)) %>% 
        mutate(horizon = ifelse(grepl("-M-", sample), "M", "O")) %>% 
        mutate(without_horizon = gsub("-[M|O]-", "-", sample)) %>% 
        mutate(plot_date = paste0(plotID, "-", dateID)) %>% 
        as.data.frame()
    rownames(df) <- make.unique(sampleID)
    return(df)
}
# 


# if (!exists("ref_data")){
# 	ref_data <- get_reference_data()
# }
# mets_recode = ref_data$chem_xref %>%
# 	filter(ID %in% gsub("\\_e","",media_in$metabolite) &
# 				 	source %in% c("biggM")) %>%
# 	distinct(source_id, .keep_all = T)
# full_recode_mets = paste0(mets_recode$source_id,"_e") %>% gsub("^M_","",.)
# names(full_recode_mets) = paste0(mets_recode$ID,"_e")
# 
# saveRDS(full_recode_mets, here("reference_data","recode_mets.rds"))

#full_recode_mets <- readRDS(here("reference_data","recode_mets.rds"))


# # # Load list of MetanetX-deprecated identifiers
# deprecated = data.table::fread(here("reference_data","chem_depr.tsv"), fill=T, skip=351, nThread = 10)
# deprecated = deprecated %>% 
# 	filter(!ID %in% `#deprecated_ID`) %>% 
# 	distinct(`#deprecated_ID`, .keep_all=T) %>% 
# 	distinct(ID, .keep_all=T)
# 
# # recode with or without the _e suffix
# deprecated_recode = c(deprecated$ID, paste0(deprecated$ID,"_e"))
# names(deprecated_recode) = c(deprecated$`#deprecated_ID`, paste0(deprecated$`#deprecated_ID`,"_e"))
# deprecated_suffix = deprecated %>% 
# 	select(-version) %>% 
# 	mutate(`#deprecated_ID` = paste0(`#deprecated_ID`,"_e"),
# 																				 ID =	paste0(ID,"_e"))
# 
# deprecated_key = rbind(deprecated %>% select(-version), deprecated_suffix)
# 
# saveRDS(deprecated_recode, here("reference_data", "deprecated_recode_mets.rds"))
# saveRDS(deprecated_key, here("reference_data", "deprecated_recode_mets_key.rds"))

#deprecated_recode <- readRDS(here("reference_data", "deprecated_recode_mets.rds"))
#deprecated_key <- readRDS(here("reference_data", "deprecated_recode_mets_key.rds"))


recode_mets = c(MNXM1103926_e = "M02447_e", MNXM1104527_e = "malttr_e", MNXM1104679_e = "alaala_e", 
                MNXM1105029_e = "glc__D_e", MNXM1105731_e = "ala__D_e", MNXM1105765_e = "glcur_e", 
                MNXM1105810_e = "fru_e", MNXM1105842_e = "man_e", MNXM1106000_e = "malt_e", 
                MNXM1106164_e = "lys__L_e", MNXM1106762_e = "leu__L_e", MNXM1107136_e = "udcpdp_e", 
                MNXM1107769_e = "his__L_e", MNXM1107821_e = "asn__L_e", MNXM1107894_e = "murein4px4px4p_e", 
                MNXM1107902_e = "no_e", MNXM1108018_e = "h2_e", MNXM1108051_e = "rib__D_e", 
                MNXM1108092_e = "etoh_e", MNXM1108175_e = "gal_e", MNXM1108206_e = "asp__L_e", 
                MNXM114_e = "pro__L_e", MNXM118_e = "ptrc_e", MNXM124_e = "spmd_e", 
                MNXM128_e = "ca2_e", MNXM158_e = "ura_e", MNXM1734_e = "galctn__D_e", 
                MNXM174_e = "xan_e", MNXM199_e = "val__L_e", MNXM2255_e = "mn2_e", 
                MNXM23_e = "pyr_e", MNXM26_e = "ac_e", MNXM27_e = "na1_e", MNXM270_e = "ribflv_e", 
                MNXM2757_e = "btd_RR_e", MNXM282_e = "taur_e", MNXM289_e = "glyb_e", 
                MNXM29_e = "gly_e", MNXM304_e = "btn_e", MNXM320_e = "phenol_e", 
                MNXM3230_e = "1btol_e", MNXM341_e = "glcn_e", MNXM37_e = "gln__L_e", 
                MNXM371_e = "bzal_e", MNXM376_e = "hqn_e", MNXM39_e = "for_e", 
                MNXM419_e = "pydxn_e", MNXM458_e = "but_e", MNXM58_e = "so4_e", 
                MNXM5970_e = "murein5px3p_e", MNXM60_e = "h2co3_e", MNXM617_e = "fol_e", 
                MNXM653_e = "mg2_e", MNXM694_e = "ser__D_e", MNXM726092_e = "mobd_e", 
                MNXM726637_e = "co_e", MNXM726711_e = "fe2_e", MNXM728337_e = "ile__L_e", 
                MNXM729215_e = "HC02172_e", MNXM729302_e = "nh3_e", MNXM729800_e = "meoh_e", 
                MNXM730135_e = "thm_e", MNXM731166_e = "cu2_e", MNXM731834_e = "lac__D_e", 
                MNXM731835_e = "lac__L_e", MNXM732398_e = "no3_e", MNXM732471_e = "citr__L_e", 
                MNXM733031_e = "udcpp_e", MNXM733618_e = "rmn_e", MNXM734574_e = "xyl__D_e", 
                MNXM734859_e = "arab__L_e", MNXM735438_e = "o2_e", MNXM735978_e = "cl_e", 
                MNXM736226_e = "galur_e", MNXM737787_e = "ser__L_e", MNXM737824_e = "murein4px4p_e", 
                MNXM738068_e = "cys__L_e", MNXM738657_e = "murein5p4p_e", MNXM738663_e = "murein5p3p_e", 
                MNXM738804_e = "met__L_e", MNXM738834_e = "alltn_e", MNXM739527_e = "arg__L_e", 
                MNXM739590_e = "ch4_e", MNXM739668_e = "uaagmda_e", MNXM740231_e = "murein5p5p5p_e", 
                MNXM741014_e = "murein5p5p_e", MNXM741016_e = "murein4p3p_e", 
                MNXM741019_e = "murein4p4p_e", MNXM741173_e = "glu__L_e", MNXM741553_e = "trp__L_e", 
                MNXM761_e = "csn_e", MNXM87121_e = "murein5px4px4p_e", MNXM88338_e = "murein5px4p_e", 
                MNXM889_e = "tol_e", MNXM89612_e = "glyc_e", MNXM9_e = "pi_e", 
                MNXM90960_e = "cobalt2_e", MNXM95_e = "k_e", WATER_e = "h2o_e",
                "MNXM2_e" = "h2o_e", 
                "MNXM4_e" = "o2_e",
                "MNXM13_e" = "co2_e",
                "MNXM207_e" = "no3_e",
                "MNXM107_e" = "no2_e",
                "MNXM15_e" = "nh4_e",
                "MNXM714_e" = "ch4_e",
                "MNXM579_e" = "n2o_e",
                "MNXM1_e" = "h_e","MNXM25_e" = "succ_e",
                "MNXM117_e" = "urea_e",
                "MNXM738430_e" = "n2o_e")



parseSimID <- function(df) {
    
    df <- df %>%
        separate(col = "id", sep = "_", into = c(NA,"diffusion_type",NA,"diffusion",NA,"ammonium",NA,'nitrate',NA,"mois",NA,"liq_diffusion"), remove=F)
    df$scenario_label = paste0("H20: ", df$mois, ", NH4: ", 
                               df$ammonium, ", NO3: ", 
                               df$nitrate, ", D(gas): ",df$diffusion, ", D(liq): ", df$liq_diffusion)
    
    return(df)
}



library(stringr)
library(xml2)

#' SBML preprocessing with logging of changes
#' @param file_path Path to original SBML file
#' @param output_path Optional path for preprocessed file (defaults to temp file)
#' @param preserve_groups Whether to attempt to preserve group data (default: FALSE for safety)
#' @return List with preprocessed file path and detailed change log
preprocess_sbml_file <- function(file_path, output_path = NULL, preserve_groups = FALSE) {
    
    cat("=== SBML File Preprocessing ===\n")
    cat("Input file:", basename(file_path), "\n")
    cat("File size:", round(file.size(file_path) / 1024, 1), "KB\n")
    
    # Initialize change tracking
    changes_log <- list(
        original_file = file_path,
        file_size_before = file.size(file_path),
        lines_removed = 0,
        elements_removed = list(),
        attributes_removed = list(),
        warnings = list()
    )
    
    # Read the file as text
    sbml_content <- readLines(file_path, warn = FALSE)
    original_line_count <- length(sbml_content)
    changes_log$lines_before = original_line_count
    
    cat("Original file has", original_line_count, "lines\n")
    
    # Track what we're looking for and removing
    elements_to_check <- c(
        "groups:group", "listOfGroups", "groups:listOfMembers", 
        "groups:member", "fbc:geneProductAssociation", "fbc:listOfGeneProducts"
    )
    
    # Count occurrences before removal
    for (element in elements_to_check) {
        count <- sum(str_detect(sbml_content, element))
        if (count > 0) {
            changes_log$elements_found <- append(changes_log$elements_found, 
                                                 list(setNames(count, element)))
            cat("Found", count, "lines containing", element, "\n")
        }
    }
    
    # Strategy 1: Remove namespace declarations that cause issues
    namespace_removals <- 0
    
    # Remove groups namespace
    before_groups_ns <- sum(str_detect(sbml_content, 'xmlns:groups='))
    sbml_content <- str_replace_all(sbml_content, 'xmlns:groups="[^"]*"', '')
    after_groups_ns <- sum(str_detect(sbml_content, 'xmlns:groups='))
    groups_ns_removed <- before_groups_ns - after_groups_ns
    
    if (groups_ns_removed > 0) {
        changes_log$attributes_removed$groups_namespace <- groups_ns_removed
        cat("Removed", groups_ns_removed, "groups namespace declarations\n")
    }
    
    # Remove groups required attribute
    before_groups_req <- sum(str_detect(sbml_content, 'groups:required='))
    sbml_content <- str_replace_all(sbml_content, 'groups:required="[^"]*"', '')
    after_groups_req <- sum(str_detect(sbml_content, 'groups:required='))
    groups_req_removed <- before_groups_req - after_groups_req
    
    if (groups_req_removed > 0) {
        changes_log$attributes_removed$groups_required <- groups_req_removed
        cat("Removed", groups_req_removed, "groups:required attributes\n")
    }
    
    # Strategy 2: Remove problematic group elements
    if (!preserve_groups) {
        cat("\nRemoving group elements (these often cause Matrix errors):\n")
        
        # Remove entire group sections
        cleaned_content <- c()
        in_groups_section <- FALSE
        in_list_of_groups <- FALSE
        groups_depth <- 0
        removed_group_lines <- 0
        
        for (i in seq_along(sbml_content)) {
            line <- sbml_content[i]
            original_line <- line
            
            # Detect start of listOfGroups
            if (str_detect(line, "<listOfGroups")) {
                in_list_of_groups <- TRUE
                groups_depth <- str_count(line, "<[^/]") - str_count(line, "</")
                removed_group_lines <- removed_group_lines + 1
                cat("  Removing listOfGroups section starting at line", i, "\n")
                next
            }
            
            # Detect start of individual group
            if (str_detect(line, "<groups:group")) {
                in_groups_section <- TRUE
                groups_depth <- str_count(line, "<[^/]") - str_count(line, "</")
                removed_group_lines <- removed_group_lines + 1
                cat("  Removing groups:group section starting at line", i, "\n")
                next
            }
            
            # Handle lines inside group sections
            if (in_list_of_groups || in_groups_section) {
                # Update depth tracking
                groups_depth <- groups_depth + str_count(line, "<[^/]") - str_count(line, "</")
                removed_group_lines <- removed_group_lines + 1
                
                # Check if we've closed all group elements
                if (groups_depth <= 0) {
                    in_list_of_groups <- FALSE
                    in_groups_section <- FALSE
                    groups_depth <- 0
                }
                next  # Skip this line
            }
            
            # Keep non-group lines, but clean up any remaining group references
            line <- str_replace_all(line, 'groups:[^\\s>]*="[^"]*"', '')  # Remove groups: attributes
            line <- str_replace_all(line, 'groups:[^\\s>]*', '')          # Remove groups: elements
            
            cleaned_content <- c(cleaned_content, line)
        }
        
        sbml_content <- cleaned_content
        changes_log$elements_removed$group_lines <- removed_group_lines
        cat("  Total group-related lines removed:", removed_group_lines, "\n")
        
    } else {
        cat("Preserving group elements (preserve_groups = TRUE)\n")
        changes_log$warnings <- append(changes_log$warnings, 
                                       "Group elements preserved - may still cause Matrix errors")
    }
    
    # Strategy 3: Clean up gene product associations if they're problematic
    gpr_lines_before <- sum(str_detect(sbml_content, "fbc:geneProductAssociation"))
    if (gpr_lines_before > 0) {
        cat("\nFound", gpr_lines_before, "gene product association lines\n")
        
        # Only remove these if they seem to be causing issues (you can adjust this)
        # For now, we'll keep them but clean up malformed ones
        
        # Clean up empty or malformed GPR elements
        sbml_content <- str_replace_all(sbml_content, 
                                        '<fbc:geneProductAssociation[^>]*></fbc:geneProductAssociation>', 
                                        '')
        
        gpr_lines_after <- sum(str_detect(sbml_content, "fbc:geneProductAssociation"))
        gpr_removed <- gpr_lines_before - gpr_lines_after
        
        if (gpr_removed > 0) {
            changes_log$elements_removed$empty_gpr <- gpr_removed
            cat("  Removed", gpr_removed, "empty gene product associations\n")
        }
    }
    
    # Final cleanup
    final_line_count <- length(sbml_content)
    changes_log$lines_after <- final_line_count
    changes_log$lines_removed <- original_line_count - final_line_count
    
    # Write to output file
    if (is.null(output_path)) {
        output_path <- tempfile(fileext = ".xml")
    }
    
    writeLines(sbml_content, output_path)
    changes_log$output_file <- output_path
    changes_log$file_size_after <- file.size(output_path)
    changes_log$size_reduction_kb <- round((changes_log$file_size_before - changes_log$file_size_after) / 1024, 1)
    
    # Summary
    cat("\n=== Preprocessing Summary ===\n")
    cat("Lines removed:", changes_log$lines_removed, "out of", original_line_count, 
        "(", round(changes_log$lines_removed / original_line_count * 100, 1), "%)\n")
    cat("File size reduction:", changes_log$size_reduction_kb, "KB\n")
    cat("Output file:", output_path, "\n")
    
    # Validate that core elements are preserved
    validate_core_elements(sbml_content, changes_log)
    
    return(list(
        file_path = output_path,
        changes_log = changes_log
    ))
}

#' Validate that core metabolic elements are preserved after preprocessing
#' @param sbml_content Preprocessed SBML content as character vector
#' @param changes_log Changes log to update with validation results
validate_core_elements <- function(sbml_content, changes_log) {
    
    cat("\n=== Validation: Core Elements Preserved ===\n")
    
    # Check for essential SBML elements
    core_elements <- list(
        "species" = "listOfSpecies|<species",
        "reactions" = "listOfReactions|<reaction",
        "compartments" = "listOfCompartments|<compartment",
        "metabolites" = "fbc:charge|fbc:chemicalFormula",
        "bounds" = "fbc:lowerFluxBound|fbc:upperFluxBound",
        "objectives" = "fbc:objective|listOfObjectives"
    )
    
    validation_results <- list()
    
    for (element_name in names(core_elements)) {
        pattern <- core_elements[[element_name]]
        count <- sum(str_detect(sbml_content, pattern))
        validation_results[[element_name]] <- count
        
        status <- if (count > 0) "✓" else "⚠"
        cat(status, element_name, ":", count, "found\n")
        
        if (count == 0 && element_name %in% c("species", "reactions", "compartments")) {
            changes_log$warnings <- append(changes_log$warnings, 
                                           paste("WARNING: No", element_name, "found after preprocessing"))
        }
    }
    
    changes_log$validation_results <- validation_results
    
    # Check if file still looks like valid SBML
    has_sbml_root <- sum(str_detect(sbml_content, "<sbml"))
    has_model <- sum(str_detect(sbml_content, "<model"))
    
    if (has_sbml_root == 0 || has_model == 0) {
        changes_log$warnings <- append(changes_log$warnings, 
                                       "WARNING: File may not be valid SBML after preprocessing")
        cat("⚠ WARNING: File structure may be compromised\n")
    } else {
        cat("✓ SBML structure appears intact\n")
    }
}

#' Compare two SBML files to see exactly what changed
#' @param original_file Path to original SBML file
#' @param processed_file Path to processed SBML file
#' @return Detailed comparison report
compare_sbml_files <- function(original_file, processed_file) {
    
    cat("=== SBML File Comparison ===\n")
    
    # Read both files
    original_content <- readLines(original_file, warn = FALSE)
    processed_content <- readLines(processed_file, warn = FALSE)
    
    # Basic statistics
    cat("Original file:", length(original_content), "lines,", file.size(original_file), "bytes\n")
    cat("Processed file:", length(processed_content), "lines,", file.size(processed_file), "bytes\n")
    
    # Find differences
    common_lines <- intersect(original_content, processed_content)
    removed_lines <- setdiff(original_content, processed_content)
    added_lines <- setdiff(processed_content, original_content)
    
    cat("\nContent analysis:\n")
    cat("  Common lines:", length(common_lines), "\n")
    cat("  Removed lines:", length(removed_lines), "\n")
    cat("  Added lines:", length(added_lines), "\n")
    
    # Show sample of removed content
    if (length(removed_lines) > 0) {
        cat("\nSample of removed content:\n")
        sample_removed <- head(removed_lines[nchar(removed_lines) > 10], 5)
        for (i in seq_along(sample_removed)) {
            cat("  ", i, ":", str_trunc(sample_removed[i], 80), "\n")
        }
    }
    
    # Analyze what types of elements were removed
    removed_elements <- list()
    if (length(removed_lines) > 0) {
        removed_elements$groups <- sum(str_detect(removed_lines, "group"))
        removed_elements$gene_products <- sum(str_detect(removed_lines, "geneProduct"))
        removed_elements$annotations <- sum(str_detect(removed_lines, "annotation"))
        removed_elements$namespaces <- sum(str_detect(removed_lines, "xmlns:"))
    }
    
    cat("\nTypes of removed elements:\n")
    for (element_type in names(removed_elements)) {
        if (removed_elements[[element_type]] > 0) {
            cat("  ", element_type, ":", removed_elements[[element_type]], "lines\n")
        }
    }
    
    return(list(
        original_lines = length(original_content),
        processed_lines = length(processed_content),
        removed_lines = length(removed_lines),
        removed_content = removed_lines,
        removed_elements = removed_elements
    ))
}

