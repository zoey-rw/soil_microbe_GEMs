# Enhanced Generalized SBML Processing Pipeline
# Handles all annotation patterns found in the 60-species database

library(here)
library(stringr)
library(tidyr)
library(minval)
library(sybilSBML)
library(tidyverse)
library(yaml)
library(jsonlite)

#' Detect the annotation pattern for a given SBML model
detect_annotation_pattern <- function(sbml_model) {
    
    # Check if annotation field exists
    has_annotation <- !is.null(sbml_model@met_attr) && 
        "annotation" %in% names(sbml_model@met_attr) &&
        !all(is.na(sbml_model@met_attr$annotation))
    
    if (!has_annotation) {
        return(list(
            pattern = "simple_bigg",
            databases = c(),
            format = "none",
            description = "No annotation field - using direct metabolite ID matching"
        ))
    }
    
    annotations <- sbml_model@met_attr$annotation
    annotations <- annotations[!is.na(annotations)]
    sample_annotations <- head(annotations, 20)
    
    # Detect RDF format
    has_rdf <- any(str_detect(sample_annotations, "<rdf:RDF|xmlns:rdf"))
    annotation_format <- ifelse(has_rdf, "rdf", "string")
    
    databases <- c()
    
    if (has_rdf) {
        cat("Detected RDF annotation format\n")
        
        # Parse sample annotations to detect all databases
        for (annotation in sample_annotations[1:min(5, length(sample_annotations))]) {
            rdf_data <- parse_rdf_annotation_robust(annotation)
            databases <- unique(c(databases, names(rdf_data)))
        }
        
        cat("Detected databases:", paste(databases, collapse = ", "), "\n")
        
    } else {
        # Original string-based detection
        has_metanetx <- any(str_detect(sample_annotations, "metanetx\\.chemical/"))
        has_bigg <- any(str_detect(sample_annotations, "bigg\\.metabolite/"))
        has_seed <- any(str_detect(sample_annotations, "seed\\.compound/"))
        has_kegg <- any(str_detect(sample_annotations, "kegg\\.compound/"))
        has_chebi <- any(str_detect(sample_annotations, "CHEBI:|chebi\\.compound/"))
        has_pubchem <- any(str_detect(sample_annotations, "pubchem\\.compound/"))
        has_inchi <- any(str_detect(sample_annotations, "metaid_.*_inchi"))
        
        if (has_metanetx) databases <- c(databases, "metanetx")
        if (has_bigg) databases <- c(databases, "bigg") 
        if (has_seed) databases <- c(databases, "seed")
        if (has_kegg) databases <- c(databases, "kegg")
        if (has_chebi) databases <- c(databases, "chebi")
        if (has_pubchem) databases <- c(databases, "pubchem")
        if (has_inchi) databases <- c(databases, "inchi")
    }
    
    # Pattern determination with RDF support
    pattern <- if (annotation_format == "rdf" && length(databases) >= 3) {
        "rdf_multi_database"
    } else if (length(databases) >= 4) {
        paste(sort(databases), collapse = "_")
    } else if ("metanetx" %in% databases && "bigg" %in% databases && "seed" %in% databases) {
        "metanetx_bigg_seed"
    } else if ("metanetx" %in% databases && "bigg" %in% databases) {
        "metanetx_bigg"
    } else if ("metanetx" %in% databases && "seed" %in% databases) {
        "metanetx_seed"
    } else if ("bigg" %in% databases && "seed" %in% databases) {
        "bigg_seed"
    } else if ("seed" %in% databases) {
        "seed_only"
    } else {
        "simple_bigg"
    }
    
    description <- paste("Format:", annotation_format, "| Databases:", paste(databases, collapse = ", "))
    
    return(list(
        pattern = pattern,
        databases = databases,
        format = annotation_format,
        description = description
    ))
}

parse_rdf_annotation_robust <- function(rdf_annotation) {
    
    if (is.na(rdf_annotation) || rdf_annotation == "" || length(rdf_annotation) == 0) {
        return(list())
    }
    
    result <- list()
    
    tryCatch({
        # Extract all identifiers.org resources using a more robust pattern
        resources <- str_extract_all(rdf_annotation, 'http://identifiers\\.org/[^"\\s>]+')[[1]]
        
        if (length(resources) == 0) {
            # Fallback: try broader pattern
            resources <- str_extract_all(rdf_annotation, 'rdf:resource="([^"]+)"')[[1]]
            if (length(resources) > 0) {
                # Clean up the rdf:resource=" part
                resources <- str_replace(resources, 'rdf:resource="', '')
                resources <- str_replace(resources, '"$', '')
            }
        }
        
        # Process each resource
        for (resource in resources) {
            if (str_detect(resource, "identifiers\\.org/")) {
                parts <- str_split(resource, "/")[[1]]
                if (length(parts) >= 2) {
                    database_raw <- parts[length(parts) - 1]
                    id_raw <- parts[length(parts)]
                    
                    # Clean up any trailing characters
                    id_raw <- str_replace(id_raw, '[">\\s].*$', '')
                    
                    # Map database names - only keep the ones we use
                    if (database_raw == "metanetx.chemical") {
                        result$metanetx <- id_raw
                    } else if (database_raw == "bigg.metabolite") {
                        result$bigg <- id_raw
                    } else if (database_raw == "seed.compound") {
                        result$seed <- id_raw
                    } else if (database_raw == "kegg.compound") {
                        result$kegg <- id_raw
                    } else if (str_detect(database_raw, "chebi")) {
                        clean_id <- str_replace(id_raw, "CHEBI:", "")
                        result$chebi <- clean_id
                    } else if (database_raw == "hmdb") {
                        result$hmdb <- id_raw
                    }
                }
            }
        }
        
    }, error = function(e) {
        # Silently return empty list on parsing errors
        return(list())
    })
    
    return(result)
}



#' Extract metabolite annotations based on detected pattern
extract_metabolite_annotations <- function(sbml_model, pattern_info) {
    
    original_met_count <- length(sbml_model@met_id)
    cat("Processing", original_met_count, "metabolites\n")
    
    # Create base data frame - this MUST have exactly 698 rows
    met_df <- data.frame(
        orig_met = sbml_model@met_id,
        met_name = if(length(sbml_model@met_name) > 0) sbml_model@met_name else rep(NA, original_met_count),
        compartment_no = sbml_model@met_comp,
        stringsAsFactors = FALSE
    )
    
    cat("Base data frame created with", nrow(met_df), "rows\n")
    
    # Compartment processing
    compart_key <- sbml_model@mod_compart
    names(compart_key) <- 1:length(sbml_model@mod_compart)
    
    compart_key <- recode(compart_key,
                          "Cytosol" = "c", "Cytoplasm" = "c", "cytosol" = "c", 
                          "extracellular space" = "e", "extracellular" = "e",
                          "Periplasm" = "p", "periplasm" = "p")
    
    compart_key <- gsub("0$", "", compart_key)
    met_df$compartment <- recode(met_df$compartment_no, !!!compart_key)
    met_df$without_compartment <- removeCompartment(met_df$orig_met)
    
    cat("After compartment processing:", nrow(met_df), "rows\n")
    
    # Process annotations only if they exist
    if (pattern_info$pattern != "simple_bigg" && 
        !is.null(sbml_model@met_attr) && 
        "annotation" %in% names(sbml_model@met_attr)) {
        
        annotation_field <- sbml_model@met_attr$annotation
        
        if (pattern_info$format == "rdf") {
            cat("Processing RDF annotations...\n")
            
            # Pre-allocate all database columns with NAs - SAME LENGTH AS met_df
            database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi", "hmdb")
            
            for (db in database_columns) {
                met_df[[db]] <- rep(NA_character_, nrow(met_df))
            }
            
            cat("Added database columns, now have", nrow(met_df), "rows\n")
            
            # Process annotations using direct indexing (no loops that could change length)
            cat("Processing", length(annotation_field), "annotations...\n")
            
            # Use vectorized approach instead of loop to prevent row issues
            for (i in 1:min(length(annotation_field), nrow(met_df))) {
                
                annotation <- annotation_field[i]
                if (!is.na(annotation) && annotation != "") {
                    rdf_data <- parse_rdf_annotation_robust_fixed(annotation)
                    
                    # Direct assignment by index to prevent row changes
                    for (db in names(rdf_data)) {
                        if (db %in% database_columns && !is.na(rdf_data[[db]])) {
                            met_df[i, db] <- rdf_data[[db]]
                        }
                    }
                }
                
                # Progress indicator for large datasets
                if (i %% 100 == 0) {
                    cat("  Processed", i, "annotations, data frame still has", nrow(met_df), "rows\n")
                }
            }
            
            cat("RDF processing complete, final row count:", nrow(met_df), "\n")
            
        } else {
            # String format processing - use the approach that was working before
            cat("Processing string-format annotations...\n")
            
            # Create a temporary annotation column
            met_df_temp <- met_df
            met_df_temp$annotation <- annotation_field
            
            cat("Added annotation column, now have", nrow(met_df_temp), "rows\n")
            
            # Use the original separate approach but capture any row changes
            original_separate_approach <- function(df) {
                suppressWarnings({
                    if ("metanetx" %in% pattern_info$databases) {
                        df <- df %>%
                            separate(col = annotation, sep = "metanetx\\.chemical/", 
                                     into = c(NA, "metanetx", NA), remove = FALSE, fill = "right", extra = "drop") %>%
                            separate(col = metanetx, sep = ";", into = c("metanetx", NA), fill = "right", extra = "drop")
                    }
                    
                    if ("bigg" %in% pattern_info$databases) {
                        df <- df %>%
                            separate(col = annotation, sep = "bigg\\.metabolite/", 
                                     into = c(NA, "bigg", NA), remove = FALSE, fill = "right", extra = "drop") %>%
                            separate(col = bigg, sep = ";", into = c("bigg", NA), fill = "right", extra = "drop")
                    }
                    
                    if ("seed" %in% pattern_info$databases) {
                        df <- df %>%
                            separate(col = annotation, sep = "seed\\.compound/", 
                                     into = c(NA, "seed", NA), remove = FALSE, fill = "right", extra = "drop") %>%
                            separate(col = seed, sep = ";", into = c("seed", NA), fill = "right", extra = "drop")
                    }
                    
                    # Remove annotation column
                    df <- df %>% select(-annotation)
                })
                return(df)
            }
            
            met_df <- original_separate_approach(met_df_temp)
            cat("String processing complete, final row count:", nrow(met_df), "\n")
        }
    }
    
    # Final validation
    if (nrow(met_df) != original_met_count) {
        stop("CRITICAL: Row count changed from ", original_met_count, " to ", nrow(met_df), 
             " during processing. This will cause the replacement error.")
    }
    
    cat("✓ Annotation extraction completed successfully with", nrow(met_df), "rows\n")
    return(met_df)
}


#' Convert database IDs to MetanetX using reference data
convert_to_metanetx <- function(met_df, pattern_info, ref_data) {
    
    if (missing(ref_data) || is.null(ref_data)) {
        stop("ref_data is missing or NULL")
    }
    
    if (!"chem_xref" %in% names(ref_data)) {
        stop("ref_data must contain 'chem_xref' component")
    }
    
    # Original database conversions
    if ("bigg" %in% names(met_df) && !all(is.na(met_df$bigg))) {
        cat("  Converting BiGG IDs...\n")
        bigg_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$bigg[!is.na(met_df$bigg)], !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_bigg <- bigg_recode[match(met_df$bigg, bigg_recode$source_id), ]$ID
    }
    
    if ("seed" %in% names(met_df) && !all(is.na(met_df$seed))) {
        cat("  Converting SEED IDs...\n")
        seed_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$seed[!is.na(met_df$seed)], !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_seed <- seed_recode[match(met_df$seed, seed_recode$source_id), ]$ID
    }
    
    if ("kegg" %in% names(met_df) && !all(is.na(met_df$kegg))) {
        cat("  Converting KEGG IDs...\n")
        kegg_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$kegg[!is.na(met_df$kegg)], 
                   source == "kegg.compound", !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_kegg <- kegg_recode[match(met_df$kegg, kegg_recode$source_id), ]$ID
    }
    
    if ("chebi" %in% names(met_df) && !all(is.na(met_df$chebi))) {
        cat("  Converting ChEBI IDs...\n")
        # Handle ChEBI IDs with or without CHEBI: prefix
        chebi_with_prefix <- paste0("CHEBI:", met_df$chebi[!is.na(met_df$chebi)])
        
        chebi_recode <- ref_data$chem_xref %>%
            filter((source_id %in% chebi_with_prefix | source_id %in% met_df$chebi[!is.na(met_df$chebi)]), 
                   source == "chebi", !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        
        # Try both formats
        met_df$new_met_chebi <- chebi_recode[match(chebi_with_prefix, chebi_recode$source_id), ]$ID
        if (all(is.na(met_df$new_met_chebi))) {
            met_df$new_met_chebi <- chebi_recode[match(met_df$chebi, chebi_recode$source_id), ]$ID
        }
    }
    
    # Add new database conversions for RDF files
    if ("hmdb" %in% names(met_df) && !all(is.na(met_df$hmdb))) {
        cat("  Converting HMDB IDs...\n")
        hmdb_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$hmdb[!is.na(met_df$hmdb)], 
                   source == "hmdb", !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_hmdb <- hmdb_recode[match(met_df$hmdb, hmdb_recode$source_id), ]$ID
    }
    
    if ("pubchem" %in% names(met_df) && !all(is.na(met_df$pubchem))) {
        cat("  Converting PubChem IDs...\n")
        pubchem_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$pubchem[!is.na(met_df$pubchem)], 
                   source == "pubchem.compound", !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_pubchem <- pubchem_recode[match(met_df$pubchem, pubchem_recode$source_id), ]$ID
    }
    
    # Simple BiGG-style conversion
    if (pattern_info$pattern == "simple_bigg") {
        cat("  Converting simple metabolite IDs...\n")
        simple_recode <- ref_data$chem_xref %>%
            filter(source_id %in% met_df$without_compartment[!is.na(met_df$without_compartment)], 
                   !is.na(source_id)) %>%
            distinct(source_id, .keep_all = TRUE)
        met_df$new_met_simple <- simple_recode[match(met_df$without_compartment, simple_recode$source_id), ]$ID
    }
    
    return(met_df)
}


#' Apply prioritization rules based on annotation pattern
apply_prioritization <- function(met_df, pattern_info) {
    
    pattern <- pattern_info$pattern
    met_df$new_met <- met_df$without_compartment  # Default fallback
    
    # Apply priority order: MetanetX (direct) > BiGG > SEED > KEGG > ChEBI > HMDB > PubChem > Original
    
    if ("new_met_pubchem" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_pubchem, met_df$new_met)
    }
    
    if ("new_met_hmdb" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_hmdb, met_df$new_met)
    }
    
    if ("new_met_chebi" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_chebi, met_df$new_met)
    }
    
    if ("new_met_kegg" %in% names(met_df)) {
        if (pattern == "metanetx_kegg" && "metanetx" %in% names(met_df)) {
            met_df$new_met <- coalesce(met_df$new_met_kegg, met_df$metanetx, met_df$new_met)
        } else {
            met_df$new_met <- coalesce(met_df$new_met_kegg, met_df$new_met)
        }
    }
    
    if ("new_met_seed" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_seed, met_df$new_met)
    }
    
    if ("new_met_bigg" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_bigg, met_df$new_met)
    }
    
    if ("metanetx" %in% names(met_df)) {
        if (pattern != "metanetx_kegg") {
            met_df$new_met <- coalesce(met_df$metanetx, met_df$new_met)
        }
    }
    
    if ("new_met_simple" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_simple, met_df$new_met)
    }
    
    return(met_df)
}

#' Apply deprecated ID replacements
replace_deprecated_ids <- function(met_df, deprecated_recode) {
    if (length(deprecated_recode) > 0) {
        met_df$new_met <- recode(met_df$new_met, !!!deprecated_recode)
    }
    return(met_df)
}

#' Handle duplicate metabolite IDs
handle_duplicates <- function(met_df) {
    
    met_df$new_met_name <- ifelse(!is.na(met_df$new_met), met_df$new_met, met_df$without_compartment)
    met_df$new_met_out <- paste0(met_df$new_met_name, "[", met_df$compartment, "]")
    met_df$old_met_out <- paste0(met_df$without_compartment, "[", met_df$compartment, "]")
    
    duped_mets <- met_df[duplicated(met_df$new_met_out), ]$new_met_name
    met_df$new_met_out <- ifelse(met_df$new_met_name %in% duped_mets, 
                                 met_df$old_met_out, 
                                 met_df$new_met_out)
    
    problematic_metabolites <- c("E", "Cx")
    for (prob_met in problematic_metabolites) {
        if (any(met_df$without_compartment == prob_met)) {
            met_df$new_met_out[met_df$without_compartment == prob_met] <- 
                met_df$orig_met[met_df$without_compartment == prob_met]
        }
    }
    
    return(met_df)
}

#' Update exchange reactions to match new metabolite IDs
update_exchange_reactions <- function(sbml_model) {
    
    tryCatch({
        exchReactDF <- findExchReact(sbml_model)
        exchReactDF <- exchReactDF[exchReactDF@uptake]
        
        for (i in 1:length(exchReactDF@react_id)) {
            react_id_index <- which(sbml_model@react_id == exchReactDF@react_id[[i]])
            met_name <- removeCompartment(exchReactDF@met_id[[i]])
            new_react_id <- paste0("EX_", met_name, "_e")
            sbml_model@react_id[[react_id_index]] <- new_react_id
        }
    }, error = function(e) {
        cat("Warning: Could not update exchange reactions:", e$message, "\n")
    })
    
    return(sbml_model)
}

#' Standardize compartment notation in reaction IDs
standardize_compartment_notation <- function(sbml_model) {
    
    sbml_model@react_id <- gsub("(e)", "_e", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(p)", "_p", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(c)", "_c", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(m)", "_m", sbml_model@react_id, fixed = TRUE)
    
    return(sbml_model)
}

#' Clean GPR (Gene-Protein-Reaction) associations
clean_gpr_associations <- function(sbml_model) {
    
    if (length(sbml_model@gpr) > 0) {
        sbml_model@gpr <- gsub("() and ", "", sbml_model@gpr, fixed = TRUE)
        sbml_model@gpr <- gsub("() or ", "", sbml_model@gpr, fixed = TRUE)
        sbml_model@gpr <- gsub("and ()", "", sbml_model@gpr, fixed = TRUE)
        sbml_model@gpr <- gsub("or ()", "", sbml_model@gpr, fixed = TRUE)
        
        for (i in 1:length(sbml_model@gpr)) {
            gpr <- sbml_model@gpr[i]
            
            if (!is.na(gpr) && gpr != "") {
                open_parens <- str_count(gpr, "\\(")
                close_parens <- str_count(gpr, "\\)")
                
                if (open_parens > close_parens) {
                    sbml_model@gpr[i] <- paste0(gpr, paste(rep(")", open_parens - close_parens), collapse = ""))
                } else if (close_parens > open_parens) {
                    extra_close <- close_parens - open_parens
                    gpr_cleaned <- gpr
                    for (j in 1:extra_close) {
                        gpr_cleaned <- str_replace(gpr_cleaned, "\\)$", "")
                    }
                    sbml_model@gpr[i] <- gpr_cleaned
                }
                
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "\\s+\\)", ")")
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "\\(\\s+", "(")
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "\\)\\s*\\)", ")")
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "\\(\\s*\\(", "(")
                
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "^\\s*(and|or)\\s+", "")
                sbml_model@gpr[i] <- str_replace_all(sbml_model@gpr[i], "\\s+(and|or)\\s*$", "")
                
                sbml_model@gpr[i] <- str_trim(sbml_model@gpr[i])
                
                if (sbml_model@gpr[i] %in% c("", "()", "( )", " ")) {
                    sbml_model@gpr[i] <- ""
                }
            }
        }
    }
    
    return(sbml_model)
}


#' Main processing function for a single species
#' @param species_dir Path to species directory
#' @param ref_data Reference data from get_reference_data()
#' @param deprecated_recode Named vector of deprecated ID replacements
#' @param config Processing configuration list
#' @return List with processing results and metadata
process_single_species <- function(species_dir, ref_data, deprecated_recode, config = list()) {
    
    model_id <- basename(species_dir)
    actual_model_id <- str_extract(model_id, "[^_]+$")  # Extract iAA1300 from azotobacter_vinelandii_iAA1300
    cat("\n=== Processing", model_id, "===\n")
    
    # Validate inputs first
    if (missing(ref_data) || is.null(ref_data)) {
        return(list(success = FALSE, error = "ref_data is missing or NULL - please provide reference data from get_reference_data()"))
    }
    
    if (!"chem_xref" %in% names(ref_data)) {
        return(list(success = FALSE, error = "ref_data must contain 'chem_xref' component"))
    }
    
    if (missing(deprecated_recode) || is.null(deprecated_recode)) {
        cat("Warning: deprecated_recode not provided, skipping deprecated ID replacement\n")
        deprecated_recode <- c()  # Empty named vector
    }
    
    # Find input file - look for any XML file that's not already processed
    xml_files <- list.files(species_dir, pattern = "\\.xml$", full.names = TRUE)
    
    # Better filtering of processed files
    processed_patterns <- c("_processed", "_cobra_validated", "_modified_cobra", "COBRA-sbml3")
    processed_files <- xml_files[str_detect(basename(xml_files), paste(processed_patterns, collapse = "|"))]
    input_candidates <- setdiff(xml_files, processed_files)
    
    if (length(input_candidates) == 0) {
        return(list(success = FALSE, error = "No input XML file found"))
    }
    
    # Better file selection prioritizing model ID match
    input_file <- input_candidates[1]  # Default to first
    if (length(input_candidates) > 1) {
        # Prioritize files with model ID
        model_id_files <- input_candidates[str_detect(basename(input_candidates), actual_model_id)]
        if (length(model_id_files) > 0) {
            # Among model ID files, prefer shorter names (likely original)
            file_lengths <- nchar(basename(model_id_files))
            input_file <- model_id_files[which.min(file_lengths)]
        }
    }
    
    output_file <- file.path(species_dir, paste0(actual_model_id, "_processed.xml"))
    
    # Processing metadata with problem tracking
    start_time <- Sys.time()
    processing_log <- list(
        model_id = actual_model_id,
        species_directory = model_id,
        start_time = start_time,
        input_file = basename(input_file),
        input_file_original_name = basename(input_file),  # Track original filename
        output_file = basename(output_file),
        problems_detected = list(),
        warnings = list(),
        conversion_summary = list()
    )
    
    tryCatch({
        # Read SBML model with fallback options
        cat("Reading SBML model...\n")
        sbml_params <- config$sbml_params %||% list()
        
        sbml_model <- NULL
        
        read_attempts <- list(
            "standard" = sbml_params,
            "safe_mode" = list(validateSBML = FALSE, bndCond = FALSE),
            "minimal" = list(validateSBML = FALSE, bndCond = FALSE, mergeMet = FALSE, balanceReact = FALSE),
            "preprocessed" = "use_preprocessing"  # New fallback
        )
        
        for (attempt_name in names(read_attempts)) {
            tryCatch({
                if (attempt_name == "preprocessed") {
                    # Use preprocessing for problematic files
                    cat("    Trying with file preprocessing...\n")
                    preprocess_result <- preprocess_sbml_file_enhanced(input_file)
                    sbml_model <- do.call(readSBMLmod, c(list(file = preprocess_result$file_path), 
                                                         list(validateSBML = FALSE, bndCond = FALSE)))
                    unlink(preprocess_result$file_path)  # Clean up temp file
                    processing_log$used_preprocessing <- TRUE
                    processing_log$preprocessing_changes <- preprocess_result$changes_log
                } else {
                    sbml_model <- do.call(readSBMLmod, c(list(file = input_file), read_attempts[[attempt_name]]))
                }
                
                processing_log$sbml_read_mode <- attempt_name
                if (attempt_name != "standard") {
                    processing_log$warnings <- append(processing_log$warnings, 
                                                      paste("Used", attempt_name, "SBML reading mode"))
                }
                break
            }, error = function(e) {
                if (attempt_name == names(read_attempts)[length(read_attempts)]) {
                    stop("All SBML read attempts failed: ", e$message)
                }
            })
        }
     
        
        if (is.null(sbml_model)) {
            stop("Failed to read SBML model with any method")
        }
        
        # Check for basic SBML issues
        if (length(sbml_model@met_id) == 0) {
            processing_log$problems_detected <- append(processing_log$problems_detected, "No metabolites found in model")
        }
        if (length(sbml_model@react_id) == 0) {
            processing_log$problems_detected <- append(processing_log$problems_detected, "No reactions found in model")
        }
        
        # Detect annotation pattern
        cat("Detecting annotation pattern...\n")
        pattern_info <- detect_annotation_pattern(sbml_model)
        cat("Pattern detected:", pattern_info$pattern, "\n")
        cat("Description:", pattern_info$description, "\n")
        
        processing_log$pattern_detected <- pattern_info$pattern
        processing_log$databases_found <- pattern_info$databases
        
        # Extract annotations
        cat("Extracting metabolite annotations...\n")
        met_df <- extract_metabolite_annotations(sbml_model, pattern_info)
        processing_log$metabolites_total <- nrow(met_df)
        
        # Convert to MetanetX with memory management
        cat("Converting to MetanetX IDs...\n")
        
        # For large datasets, process in batches to prevent memory crashes
        if (nrow(met_df) > 3000) {
            cat("Large dataset detected (", nrow(met_df), "metabolites), using batch processing...\n")
            processing_log$used_batch_processing <- TRUE
            
            batch_size <- 1000
            batches <- split(1:nrow(met_df), ceiling(seq_along(1:nrow(met_df)) / batch_size))
            
            for (i in seq_along(batches)) {
                cat("  Processing batch", i, "of", length(batches), "\n")
                batch_indices <- batches[[i]]
                met_df_batch <- met_df[batch_indices, ]
                met_df_batch <- convert_to_metanetx(met_df_batch, pattern_info, ref_data)
                met_df[batch_indices, ] <- met_df_batch
                gc()  # Force garbage collection between batches
            }
        } else {
            met_df <- convert_to_metanetx(met_df, pattern_info, ref_data)
        }
        
        # Apply prioritization
        cat("Applying prioritization rules...\n")
        met_df <- apply_prioritization(met_df, pattern_info)
        
        # Replace deprecated IDs
        cat("Replacing deprecated IDs...\n")
        met_df_before_deprecated <- met_df
        met_df <- replace_deprecated_ids(met_df, deprecated_recode)
        
        # Track deprecated replacements
        if (length(deprecated_recode) > 0) {
            deprecated_replacements <- sum(met_df$new_met != met_df_before_deprecated$new_met, na.rm = TRUE)
            if (deprecated_replacements > 0) {
                processing_log$conversion_summary$deprecated_ids_replaced <- deprecated_replacements
            }
        }
        
        # Handle duplicates
        cat("Handling duplicate IDs...\n")
        met_df <- handle_duplicates(met_df)
        
        # Track duplicate issues
        duplicates_reverted <- sum(met_df$new_met_out == met_df$old_met_out)
        if (duplicates_reverted > 0) {
            processing_log$problems_detected <- append(processing_log$problems_detected, 
                                                       paste(duplicates_reverted, "metabolites reverted to original IDs due to duplicates"))
            processing_log$conversion_summary$duplicate_ids_reverted <- duplicates_reverted
        }
        
        # Track conversion success rates
        processing_log$metabolites_standardized <- sum(!is.na(met_df$new_met) & met_df$new_met != met_df$without_compartment)
        processing_log$conversion_summary$total_metabolites <- nrow(met_df)
        processing_log$conversion_summary$successfully_converted <- processing_log$metabolites_standardized
        processing_log$conversion_summary$conversion_rate <- round(processing_log$metabolites_standardized / nrow(met_df) * 100, 1)
        
        # Track metabolites that couldn't be converted
        unconverted <- sum(is.na(met_df$new_met) | met_df$new_met == met_df$without_compartment)
        if (unconverted > 0) {
            processing_log$conversion_summary$unconverted_metabolites <- unconverted
            processing_log$problems_detected <- append(processing_log$problems_detected,
                                                       paste(unconverted, "metabolites could not be converted to MetanetX IDs"))
        }
        
        # Update model with validation
        cat("Updating SBML model...\n")
        
        # Validate data frame before updating model
        if (nrow(met_df) != length(sbml_model@met_id)) {
            stop("Metabolite count mismatch: data frame has ", nrow(met_df), " rows but model has ", length(sbml_model@met_id), " metabolites")
        }
        
        if (any(is.na(met_df$new_met_out))) {
            na_count <- sum(is.na(met_df$new_met_out))
            processing_log$problems_detected <- append(processing_log$problems_detected,
                                                       paste(na_count, "metabolites have NA in final output"))
        }
        
        sbml_model@met_id <- met_df$new_met_out
        
        # Additional processing with problem tracking
        
        # Apply model-specific compartment updates to the SBML object
        # Extract compartment key from the metabolite data frame processing
        compart_key <- sbml_model@mod_compart
        names(compart_key) <- 1:length(sbml_model@mod_compart)
        
        compart_key <- recode(compart_key,
                              "Cytosol" = "c", "Extra_organism" = "e", "Periplasm" = "p",
                              "cytosol" = "c", "extracellular" = "e",
                              "C_p" = "p", "C_m" = "m", "C_c" = "c", "C_e" = "e")
        
        compart_key <- gsub("0$", "", compart_key)
        sbml_model@mod_compart <- compart_key        
        # Handle special compartment assignments (like iJN1462)
        # Check for specific metabolites that need compartment corrections
        if (any(str_detect(sbml_model@met_id, "acmtsoxin"))) {
            # Find and fix specific metabolite compartment assignments
            acmt_e_idx <- which(sbml_model@met_id == "MNXM1092516[e]")
            acmt_p_idx <- which(sbml_model@met_id == "MNXM1092516[p]")
            
            if (length(acmt_e_idx) > 0) sbml_model@met_comp[acmt_e_idx] <- as.integer(2)
            if (length(acmt_p_idx) > 0) sbml_model@met_comp[acmt_p_idx] <- as.integer(3)
        }
        
        # Check GPR issues before cleaning
        if (length(sbml_model@gpr) > 0) {
            gpr_issues <- 0
            for (gpr in sbml_model@gpr) {
                if (!is.na(gpr) && gpr != "") {
                    open_parens <- str_count(gpr, "\\(")
                    close_parens <- str_count(gpr, "\\)")
                    if (open_parens != close_parens) {
                        gpr_issues <- gpr_issues + 1
                    }
                }
            }
            if (gpr_issues > 0) {
                processing_log$problems_detected <- append(processing_log$problems_detected,
                                                           paste(gpr_issues, "GPR associations had unmatched parentheses (fixed)"))
                processing_log$conversion_summary$gpr_issues_fixed <- gpr_issues
            }
        }
        
        sbml_model <- update_exchange_reactions(sbml_model)
        sbml_model <- standardize_compartment_notation(sbml_model)
        sbml_model <- clean_gpr_associations(sbml_model)
        
        # Write output with validation
        cat("Writing processed SBML file...\n")
        writeSBML(sbml_model, level = 3, filename = output_file)
        
        # Verify output file was created
        if (!file.exists(output_file)) {
            stop("Output file was not created successfully")
        }
        
        # Processing summary
        end_time <- Sys.time()
        processing_log$end_time <- end_time
        processing_log$processing_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
        processing_log$success <- TRUE
        
        # Write processing metadata
        metadata_file <- file.path(species_dir, "processing_metadata.json")
        write_json(processing_log, metadata_file, pretty = TRUE, auto_unbox = TRUE)
        
        cat("Processing completed successfully!\n")
        cat("Processing time:", round(processing_log$processing_time_seconds, 2), "seconds\n")
        cat("Conversion rate:", processing_log$conversion_summary$conversion_rate, "%\n")
        
        if (length(processing_log$problems_detected) > 0) {
            cat("Problems detected:\n")
            for (problem in processing_log$problems_detected) {
                cat("  -", problem, "\n")
            }
        }
        
        return(list(
            success = TRUE,
            processing_log = processing_log,
            metabolite_data = met_df,
            pattern_info = pattern_info
        ))
        
    }, error = function(e) {
        processing_log$success <- FALSE
        processing_log$error <- e$message
        processing_log$end_time <- Sys.time()
        
        # Write failed processing metadata
        metadata_file <- file.path(species_dir, "processing_metadata.json")
        write_json(processing_log, metadata_file, pretty = TRUE, auto_unbox = TRUE)
        
        cat("Processing failed:", e$message, "\n")
        
        return(list(
            success = FALSE,
            error = e$message,
            processing_log = processing_log
        ))
    })
}


source("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/pipeline/process_sbml_species.R")
source("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/pipeline/sbml_processing_utils.R")

ref_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/metanetx_reference_data.rds")
deprecated_recode <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/deprecated_recode_mets.rds")


# result <- process_single_species("./species/nitrobacter_winogradskyi_iFC579", ref_data, deprecated_recode)
# result <- process_single_species("./species/methanosarcina_barkeri_iMG746", ref_data, deprecated_recode)
# result <- process_single_species("./species/nitrosomonas_europaea_iGC535/", ref_data, deprecated_recode)
# result <- process_single_species("./species/bacillus_subtilis_iBB1018/", ref_data, deprecated_recode)
# 
# # Validate growth rates
# validation <- validate_model_growth("species/methanosarcina_barkeri_iMG746")



validation <- validate_model_growth("/Users/zoeywerbin/soil_GEM_database/microbial_gem_database/species/methanosarcina_barkeri_iMG746")
