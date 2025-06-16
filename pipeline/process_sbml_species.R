# Generalized SBML Processing Pipeline
# Handles all annotation patterns found in the 60-species database

library(here)
library(stringr)
library(tidyr)
library(minval)
library(sybilSBML)
library(tidyverse)
library(yaml)
library(jsonlite)

# Source utilities
source(here("pipeline", "processing_utils.R"))
source(here("pipeline", "sbml_processing_utils.R"))

#' Determine annotation pattern from detected databases
#' @param databases Vector of detected database names
#' @param annotation_format Format type ("rdf" or "string")
#' @return Pattern string
determine_pattern_from_databases <- function(databases, annotation_format) {
    
    if (annotation_format == "rdf" && length(databases) >= 3) {
        return("rdf_multi_database")
    } else if (length(databases) >= 4) {
        return(paste(sort(databases), collapse = "_"))
    } else if ("metanetx" %in% databases && "bigg" %in% databases && "seed" %in% databases) {
        return("metanetx_bigg_seed")
    } else if ("metanetx" %in% databases && "bigg" %in% databases) {
        return("metanetx_bigg")
    } else if ("metanetx" %in% databases && "seed" %in% databases) {
        return("metanetx_seed")
    } else if ("bigg" %in% databases && "seed" %in% databases) {
        return("bigg_seed")
    } else if ("seed" %in% databases) {
        return("seed_only")
    } else {
        return("simple_bigg")
    }
}

#' Detect the annotation pattern for a given SBML model
#' @param sbml_model sybilSBML model object
#' @return List with pattern information
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
        # String-based detection
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
    
    # Pattern determination
    pattern <- determine_pattern_from_databases(databases, annotation_format)
    description <- paste("Format:", annotation_format, "| Databases:", paste(databases, collapse = ", "))
    
    return(list(
        pattern = pattern,
        databases = databases,
        format = annotation_format,
        description = description
    ))
}

#' Extract metabolite annotations with compartment recovery for inconsistent SBML
extract_metabolite_annotations <- function(sbml_model, pattern_info) {
    
    original_met_count <- length(sbml_model@met_id)
    cat("Processing", original_met_count, "metabolites\n")
    
    # Fix encoding issues in original metabolite IDs
    cat("  Checking for encoding issues in metabolite IDs...\n")
    original_met_ids <- sbml_model@met_id
    fixed_met_ids <- fix_metabolite_encoding(original_met_ids)
    
    encoding_fixes_applied <- sum(original_met_ids != fixed_met_ids)
    if (encoding_fixes_applied > 0) {
        cat("  🔧 Fixed encoding issues in", encoding_fixes_applied, "metabolite IDs\n")
        sbml_model@met_id <- fixed_met_ids
    }
    
    # Create base data frame - MUST maintain exact row count
    met_df <- data.frame(
        orig_met = original_met_ids,  # Keep original for logging
        fixed_met = fixed_met_ids,    # Store fixed version
        met_name = if(length(sbml_model@met_name) > 0) sbml_model@met_name else rep(NA, original_met_count),
        compartment_no = sbml_model@met_comp,
        stringsAsFactors = FALSE
    )
    
    cat("Base data frame created with", nrow(met_df), "rows\n")
    
    # Check for metabolites with missing compartment_no
    na_comp_no <- is.na(met_df$compartment_no)
    if (any(na_comp_no)) {
        cat("  Found", sum(na_comp_no), "metabolites with missing compartment_no\n")
        cat("  Attempting to recover compartment info from metabolite IDs...\n")
        
        # Try to extract compartment from metabolite ID
        for (i in which(na_comp_no)) {
            met_id <- met_df$fixed_met[i]
            
            # Extract compartment from patterns like [e], [p], [c], etc.
            comp_match <- str_extract(met_id, "\\[([a-zA-Z0-9]+)\\]$")
            if (!is.na(comp_match)) {
                # Remove brackets to get compartment letter
                comp_letter <- str_replace_all(comp_match, "\\[|\\]", "")
                
                # Map compartment letter back to compartment number
                # Find which compartment index corresponds to this letter
                compart_key <- standardize_compartments(sbml_model)
                matching_comp_no <- names(compart_key)[compart_key == comp_letter]
                
                if (length(matching_comp_no) > 0) {
                    met_df$compartment_no[i] <- as.numeric(matching_comp_no[1])
                    cat("    Recovered compartment for", met_id, ": compartment_no =", matching_comp_no[1], "\n")
                } else {
                    cat("    Could not map compartment letter '", comp_letter, "' for", met_id, "\n")
                }
            }
        }
    }
    
    # Enhanced compartment processing
    cat("  Processing compartments...\n")
    compart_key <- standardize_compartments(sbml_model)
    
    # Apply compartment mapping with validation
    met_df$compartment <- dplyr::recode(as.character(met_df$compartment_no), !!!compart_key)
    
    # Check for failed compartment mappings after recovery attempt
    na_compartments <- is.na(met_df$compartment) | met_df$compartment == ""
    if (any(na_compartments)) {
        cat("  Warning:", sum(na_compartments), "metabolites still have missing compartments after recovery\n")
        
        # For remaining failures, try to extract compartment directly from ID
        for (i in which(na_compartments)) {
            met_id <- met_df$fixed_met[i]
            
            # Try different compartment patterns
            if (grepl("\\[e\\]$|_e$", met_id)) {
                met_df$compartment[i] <- "e"
            } else if (grepl("\\[p\\]$|_p$", met_id)) {
                met_df$compartment[i] <- "p"
            } else if (grepl("\\[c\\]$|_c$", met_id)) {
                met_df$compartment[i] <- "c"
            } else if (grepl("\\[m\\]$|_m$", met_id)) {
                met_df$compartment[i] <- "m"
            } else {
                # Default to cytosol for unknown
                met_df$compartment[i] <- "c"
                cat("    No compartment info found for", met_id, ", defaulting to 'c'\n")
            }
            
            cat("    Direct extraction for", met_id, "-> compartment =", met_df$compartment[i], "\n")
        }
    }
    
    # Final validation - no compartments should be NA at this point
    final_na_compartments <- is.na(met_df$compartment) | met_df$compartment == ""
    if (any(final_na_compartments)) {
        cat("  ERROR: Still have", sum(final_na_compartments), "metabolites with missing compartments\n")
        problem_mets <- met_df[final_na_compartments, c("orig_met", "fixed_met", "compartment_no", "compartment")]
        print(head(problem_mets, 10))
        stop("CRITICAL: Could not resolve compartments for all metabolites")
    }
    
    # Use fixed metabolite IDs for compartment removal
    met_df$without_compartment <- removeCompartment(met_df$fixed_met)
    
    cat("After compartment processing:", nrow(met_df), "rows\n")
    cat("Final compartment distribution:\n")
    print(table(met_df$compartment))
    
    # Process annotations only if they exist - CRITICAL: No operations that change row count
    if (pattern_info$pattern != "simple_bigg" && 
        !is.null(sbml_model@met_attr) && 
        "annotation" %in% names(sbml_model@met_attr)) {
        
        annotation_field <- sbml_model@met_attr$annotation
        
        # Ensure annotation field has same length as metabolites
        if (length(annotation_field) != original_met_count) {
            cat("  Warning: annotation field length (", length(annotation_field), 
                ") != metabolite count (", original_met_count, "). Padding with NAs.\n")
            if (length(annotation_field) < original_met_count) {
                annotation_field <- c(annotation_field, rep(NA, original_met_count - length(annotation_field)))
            } else {
                annotation_field <- annotation_field[1:original_met_count]
            }
        }
        
        if (pattern_info$format == "rdf") {
            cat("Processing RDF annotations...\n")
            
            # Pre-allocate all database columns with NAs - EXACT same length as met_df
            database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi", "hmdb")
            
            for (db in database_columns) {
                met_df[[db]] <- rep(NA_character_, nrow(met_df))
            }
            
            cat("Added database columns, still have", nrow(met_df), "rows\n")
            
            # Process annotations by direct assignment (NO row operations)
            for (i in 1:nrow(met_df)) {
                if (i <= length(annotation_field)) {
                    annotation <- annotation_field[i]
                    if (!is.na(annotation) && annotation != "") {
                        rdf_data <- parse_rdf_annotation_robust(annotation)
                        
                        # Direct assignment by index - cannot change row count
                        for (db in names(rdf_data)) {
                            if (db %in% database_columns && !is.na(rdf_data[[db]])) {
                                met_df[i, db] <- rdf_data[[db]]
                            }
                        }
                    }
                }
                
                # Progress indicator
                if (i %% 200 == 0) {
                    cat("  Processed", i, "/", nrow(met_df), "annotations\n")
                }
            }
            
            cat("RDF processing complete, final row count:", nrow(met_df), "\n")
            
        } else {
            # String format processing - USE SAFER APPROACH
            cat("Processing string-format annotations...\n")
            
            # Initialize database columns manually to avoid separate() issues
            database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi")
            for (db in database_columns) {
                met_df[[db]] <- rep(NA_character_, nrow(met_df))
            }
            
            cat("Pre-allocated database columns, still have", nrow(met_df), "rows\n")
            
            # Manual parsing instead of separate() to avoid row changes
            for (i in 1:nrow(met_df)) {
                if (i <= length(annotation_field)) {
                    annotation <- annotation_field[i]
                    if (!is.na(annotation) && annotation != "") {
                        
                        # Extract each database ID manually
                        if ("metanetx" %in% pattern_info$databases) {
                            metanetx_match <- str_extract(annotation, "metanetx\\.chemical/([^;\\s]+)")
                            if (!is.na(metanetx_match)) {
                                met_df[i, "metanetx"] <- str_replace(metanetx_match, "metanetx\\.chemical/", "")
                            }
                        }
                        
                        if ("bigg" %in% pattern_info$databases) {
                            bigg_match <- str_extract(annotation, "bigg\\.metabolite/([^;\\s]+)")
                            if (!is.na(bigg_match)) {
                                met_df[i, "bigg"] <- str_replace(bigg_match, "bigg\\.metabolite/", "")
                            }
                        }
                        
                        if ("seed" %in% pattern_info$databases) {
                            seed_match <- str_extract(annotation, "seed\\.compound/([^;\\s]+)")
                            if (!is.na(seed_match)) {
                                met_df[i, "seed"] <- str_replace(seed_match, "seed\\.compound/", "")
                            }
                        }
                        
                        if ("kegg" %in% pattern_info$databases) {
                            kegg_match <- str_extract(annotation, "kegg\\.compound/([^;\\s]+)")
                            if (!is.na(kegg_match)) {
                                met_df[i, "kegg"] <- str_replace(kegg_match, "kegg\\.compound/", "")
                            }
                        }
                        
                        if ("chebi" %in% pattern_info$databases) {
                            chebi_match <- str_extract(annotation, "CHEBI:([^;\\s]+)")
                            if (!is.na(chebi_match)) {
                                met_df[i, "chebi"] <- str_replace(chebi_match, "CHEBI:", "")
                            }
                        }
                    }
                }
                
                # Progress indicator
                if (i %% 200 == 0) {
                    cat("  Processed", i, "/", nrow(met_df), "annotations\n")
                }
            }
            
            cat("String processing complete, final row count:", nrow(met_df), "\n")
        }
    }
    
    # Check for and apply CarveFungi mappings
    if (!is.null(carvefungi_mapping_dir)) {
        cat("Checking for CarveFungi mappings...\n")
        carvefungi_mapping <- load_carvefungi_mapping(carvefungi_mapping_dir)
        if (!is.null(carvefungi_mapping)) {
            met_df <- apply_carvefungi_mappings(met_df, carvefungi_mapping)
            
            # Update pattern info to reflect additional databases
            if (!"carvefungi_enhanced" %in% pattern_info$databases) {
                pattern_info$databases <- c(pattern_info$databases, "carvefungi_enhanced")
                pattern_info$description <- paste(pattern_info$description, "| Enhanced with CarveFungi mappings")
            }
        }
    }
    
    
    # Final validation - CRITICAL CHECK
    if (nrow(met_df) != original_met_count) {
        stop("CRITICAL: Row count changed from ", original_met_count, " to ", nrow(met_df), 
             " during processing. This will cause the replacement error.")
    }
    
    cat("✓ Annotation extraction completed successfully with", nrow(met_df), "rows\n")
    return(met_df)
}
#' Process RDF format annotations
#' @param met_df Metabolite data frame
#' @param annotation_field Annotation field vector
#' @return Updated metabolite data frame
process_rdf_annotations <- function(met_df, annotation_field) {
    
    cat("Processing RDF annotations...\n")
    
    # Pre-allocate database columns
    database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi", "hmdb")
    for (db in database_columns) {
        met_df[[db]] <- rep(NA_character_, nrow(met_df))
    }
    
    # Process annotations by direct assignment
    for (i in 1:nrow(met_df)) {
        if (i <= length(annotation_field)) {
            annotation <- annotation_field[i]
            if (!is.na(annotation) && annotation != "") {
                rdf_data <- parse_rdf_annotation_robust(annotation)
                
                # Direct assignment by index
                for (db in names(rdf_data)) {
                    if (db %in% database_columns && !is.na(rdf_data[[db]])) {
                        met_df[i, db] <- rdf_data[[db]]
                    }
                }
            }
        }
        
        if (i %% 200 == 0) {
            cat("  Processed", i, "/", nrow(met_df), "annotations\n")
        }
    }
    
    cat("RDF processing complete, final row count:", nrow(met_df), "\n")
    return(met_df)
}

#' Process string format annotations
#' @param met_df Metabolite data frame
#' @param annotation_field Annotation field vector
#' @param databases Vector of databases to extract
#' @return Updated metabolite data frame
process_string_annotations <- function(met_df, annotation_field, databases) {
    
    cat("Processing string-format annotations...\n")
    
    # Initialize database columns
    database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi")
    for (db in database_columns) {
        met_df[[db]] <- rep(NA_character_, nrow(met_df))
    }
    
    # Manual parsing to avoid row changes
    for (i in 1:nrow(met_df)) {
        if (i <= length(annotation_field)) {
            annotation <- annotation_field[i]
            if (!is.na(annotation) && annotation != "") {
                extracted_ids <- extract_string_annotation_ids(annotation, databases)
                
                # Assign extracted IDs
                for (db in names(extracted_ids)) {
                    if (db %in% database_columns) {
                        met_df[i, db] <- extracted_ids[[db]]
                    }
                }
            }
        }
        
        if (i %% 200 == 0) {
            cat("  Processed", i, "/", nrow(met_df), "annotations\n")
        }
    }
    
    cat("String processing complete, final row count:", nrow(met_df), "\n")
    return(met_df)
}

#' Convert database IDs to MetanetX
#' @param met_df Metabolite data frame
#' @param pattern_info Pattern information
#' @param ref_data Reference data
#' @return Updated metabolite data frame with MetanetX conversions
convert_to_metanetx <- function(met_df, pattern_info, ref_data) {
        
        original_row_count <- nrow(met_df)
        cat("  Starting conversion with", original_row_count, "rows\n")
        
        if (missing(ref_data) || is.null(ref_data) || !"chem_xref" %in% names(ref_data)) {
            stop("ref_data is missing, NULL, or doesn't contain 'chem_xref' component")
        }
        
        # Convert each database type using the safe utility function
        database_conversions <- list(
            list(col = "bigg", output = "new_met_bigg", source = NULL),
            list(col = "seed", output = "new_met_seed", source = NULL),
            list(col = "kegg", output = "new_met_kegg", source = "kegg.compound"),
            list(col = "chebi", output = "new_met_chebi", source = "chebi"),
            list(col = "hmdb", output = "new_met_hmdb", source = "hmdb"),
            list(col = "pubchem", output = "new_met_pubchem", source = "pubchem.compound")
        )
        
        for (conversion in database_conversions) {
            if (conversion$col %in% names(met_df)) {
                cat("  Converting", toupper(conversion$col), "IDs...\n")
                met_df <- convert_single_database_safe(
                    met_df, conversion$col, conversion$output, ref_data, conversion$source
                )
            }
        }
        
        # Simple BiGG-style conversion
        if (pattern_info$pattern == "simple_bigg") {
            cat("  Converting simple metabolite IDs...\n")
            met_df <- convert_single_database_safe(
                met_df, "without_compartment", "new_met_simple", ref_data, NULL
            )
        }
        
        validate_row_count(met_df, original_row_count, "database conversion")
        cat("  ✓ Database conversion completed with", nrow(met_df), "rows preserved\n")
        return(met_df)
}

#' Apply prioritization rules based on annotation pattern
#' @param met_df Metabolite data frame with conversions
#' @param pattern_info Pattern information
#' @return Data frame with final MetanetX IDs
apply_prioritization <- function(met_df, pattern_info) {
    
    pattern <- pattern_info$pattern
    met_df$new_met <- met_df$without_compartment  # Default fallback
    
    # Apply priority order: MetanetX (direct) > BiGG > SEED > KEGG > ChEBI > HMDB > PubChem > Original
    priority_columns <- c("new_met_pubchem", "new_met_hmdb", "new_met_chebi", 
                          "new_met_kegg", "new_met_seed", "new_met_bigg")
    
    for (col in priority_columns) {
        if (col %in% names(met_df)) {
            met_df$new_met <- coalesce(met_df[[col]], met_df$new_met)
        }
    }
    
    # Handle special case for metanetx_kegg pattern
    if (pattern == "metanetx_kegg" && "metanetx" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_kegg, met_df$metanetx, met_df$new_met)
    } else if ("metanetx" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$metanetx, met_df$new_met)
    }
    
    if ("new_met_simple" %in% names(met_df)) {
        met_df$new_met <- coalesce(met_df$new_met_simple, met_df$new_met)
    }
    
    return(met_df)
}

#' Apply deprecated ID replacements
#' @param met_df Metabolite data frame
#' @param deprecated_recode Named vector of deprecated ID replacements
#' @return Data frame with deprecated IDs replaced
replace_deprecated_ids <- function(met_df, deprecated_recode) {
    if (length(deprecated_recode) > 0) {
        met_df$new_met <- recode(met_df$new_met, !!!deprecated_recode)
    }
    return(met_df)
}

#' Handle duplicate metabolite IDs with improved encoding awareness
handle_duplicates <- function(met_df) {
    
    met_df$new_met_name <- ifelse(!is.na(met_df$new_met), met_df$new_met, met_df$without_compartment)
    met_df$new_met_out <- paste0(met_df$new_met_name, "[", met_df$compartment, "]")
    met_df$old_met_out <- paste0(met_df$without_compartment, "[", met_df$compartment, "]")
    
    # Revert duplicates to original IDs
    duped_mets <- met_df[duplicated(met_df$new_met_out), ]$new_met_name
    met_df$new_met_out <- ifelse(met_df$new_met_name %in% duped_mets, 
                                 met_df$old_met_out, 
                                 met_df$new_met_out)
    
    # Handle problematic metabolites
    problematic_metabolites <- c("E", "Cx")
    for (prob_met in problematic_metabolites) {
        if (any(met_df$without_compartment == prob_met)) {
            met_df$new_met_out[met_df$without_compartment == prob_met] <- 
                met_df$fixed_met[met_df$without_compartment == prob_met]
        }
    }
    
    return(met_df)
}


#' Main processing function for a single species with enhanced logging
#' Main processing function for a single species
#' @param species_dir Path to species directory
#' @param ref_data Reference data from get_reference_data()
#' @param deprecated_recode Named vector of deprecated ID replacements
#' @param config Processing configuration list
#' @return List with processing results and metadata
process_single_species <- function(species_dir, ref_data, deprecated_recode, config = list()) {
    
    model_id <- basename(species_dir)
    actual_model_id <- str_extract(model_id, "[^_]+$")
    cat("\n=== Processing", model_id, "===\n")
    
    # Validate inputs
    if (missing(ref_data) || is.null(ref_data)) {
        return(list(success = FALSE, error = "ref_data is missing or NULL"))
    }
    
    if (!"chem_xref" %in% names(ref_data)) {
        return(list(success = FALSE, error = "ref_data must contain 'chem_xref' component"))
    }
    
    if (missing(deprecated_recode) || is.null(deprecated_recode)) {
        cat("Warning: deprecated_recode not provided, skipping deprecated ID replacement\n")
        deprecated_recode <- c()
    }
    
    # Setup directories and find files
    log_dir <- file.path(species_dir, "conversion_logs")
    if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
    }
    
    files <- discover_sbml_files(species_dir, actual_model_id)
    input_file <- files$input_file
    output_file <- files$output_file
    
    # Initialize processing log
    start_time <- Sys.time()
    processing_log <- list(
        model_id = actual_model_id,
        species_directory = model_id,
        start_time = start_time,
        input_file = basename(input_file),
        output_file = basename(output_file),
        problems_detected = list(),
        warnings = list(),
        conversion_summary = list()
    )
    
    tryCatch({
        # Read SBML model
        cat("Reading SBML model...\n")
        sbml_result <- read_sbml_simple(input_file)
        sbml_model <- sbml_result$model
        processing_log <- c(processing_log, sbml_result$metadata)
        
        # Validate basic SBML structure
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
        carvefungi_dir <- config$carvefungi_mapping_dir %||% dirname(species_dir)
        met_df <- extract_metabolite_annotations(sbml_model, pattern_info, carvefungi_dir)
        processing_log$metabolites_total <- nrow(met_df)
        
        # Convert to MetanetX
        cat("Converting to MetanetX IDs...\n")
        met_df <- convert_to_metanetx(met_df, pattern_info, ref_data)
        
        # Apply prioritization
        cat("Applying prioritization rules...\n")
        met_df <- apply_prioritization(met_df, pattern_info)
        
        # Replace deprecated IDs
        cat("Replacing deprecated IDs...\n")
        met_df <- replace_deprecated_ids(met_df, deprecated_recode)
        
        # Handle duplicates
        cat("Handling duplicate IDs...\n")
        #met_df <- handle_duplicates(met_df)
        met_df <- handle_duplicates(met_df)
        
        # Log conversions
        cat("Generating conversion logs...\n")
        conversion_log_result <- log_metabolite_conversions(met_df, actual_model_id, log_dir)
        processing_log$conversion_log_files <- conversion_log_result
        processing_log$conversion_summary <- conversion_log_result$statistics
        
        # Update SBML model
        cat("Updating SBML model...\n")
        
        # Validate before updating
        if (nrow(met_df) != length(sbml_model@met_id)) {
            stop("Metabolite count mismatch: data frame has ", nrow(met_df), 
                 " rows but model has ", length(sbml_model@met_id), " metabolites")
        }
        
        if (any(is.na(met_df$new_met_out))) {
            processing_log$problems_detected <- append(processing_log$problems_detected,
                                                       paste(sum(is.na(met_df$new_met_out)), "metabolites have NA in final output"))
        }
        
        # Update metabolite IDs
        sbml_model@met_id <- met_df$new_met_out
        
        # Apply additional model updates
        sbml_model <- update_exchange_reactions(sbml_model)
        sbml_model <- standardize_compartment_notation(sbml_model)
        sbml_model <- clean_gpr_associations(sbml_model)
        
        # Write output
        cat("Writing processed SBML file...\n")
        writeSBML(sbml_model, level = 3, filename = output_file)
        
        if (!file.exists(output_file)) {
            stop("Output file was not created successfully")
        }
        
        # Finalize processing log
        end_time <- Sys.time()
        processing_log$end_time <- end_time
        processing_log$processing_time_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
        processing_log$success <- TRUE
        
        # Write processing metadata
        metadata_file <- file.path(species_dir, "processing_metadata.json")
        write_json(processing_log, metadata_file, pretty = TRUE, auto_unbox = TRUE)
        
        cat("✅ Processing completed successfully!\n")
        cat("Processing time:", round(processing_log$processing_time_seconds, 2), "seconds\n")
        cat("Conversion rate:", processing_log$conversion_summary$conversion_rate, "%\n")
        cat("MetanetX rate:", processing_log$conversion_summary$metanetx_rate, "%\n")
        
        if (length(processing_log$problems_detected) > 0) {
            cat("⚠️  Problems detected:\n")
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
        
        cat("❌ Processing failed:", e$message, "\n")
        
        return(list(
            success = FALSE,
            error = e$message,
            processing_log = processing_log
        ))
    })
}

# For testing with individual species 

# ref_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/metanetx_reference_data.rds")
# deprecated_recode <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/deprecated_recode_mets.rds")

# result <- process_single_species("./species/nitrobacter_winogradskyi_iFC579", ref_data, deprecated_recode)
# result <- process_single_species("./species/methanosarcina_barkeri_iMG746", ref_data, deprecated_recode)
# result <- process_single_species("./species/nitrosomonas_europaea_iGC535/", ref_data, deprecated_recode)
#result <- process_single_species("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species/bacillus_subtilis_iBB1018/", ref_data, deprecated_recode)
