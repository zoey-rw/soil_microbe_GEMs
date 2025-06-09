# SBML Processing Utilities
# Modular functions extracted from main processing pipeline

library(stringr)
library(tidyverse)
library(jsonlite)
library(sybilSBML)

# =============================================================================
# CHARACTER ENCODING AND VALIDATION UTILITIES
# =============================================================================

#' Fix character encoding issues in metabolite IDs
#' @param metabolite_ids Vector of metabolite IDs that may have encoding issues
#' @return Vector of metabolite IDs with fixed encoding
fix_metabolite_encoding <- function(metabolite_ids) {
    
    # Create a comprehensive mapping of common encoding issues
    encoding_fixes <- list(
        # Bracket encoding issues
        "__91__" = "[",
        "__93__" = "]", 
        "_LSQB_" = "[",
        "_RSQB_" = "]",
        
        # Parentheses encoding issues  
        "_LPAREN_" = "(",
        "_RPAREN_" = ")",
        "__40__" = "(",
        "__41__" = ")",
        
        # Dash/hyphen encoding issues
        "_DASH_" = "-",
        "_MINUS_" = "-",
        "__45__" = "-",
        
        # Other common encoding issues
        "_DOT_" = ".",
        "_COMMA_" = ",",
        "_COLON_" = ":",
        "_SEMICOLON_" = ";",
        "_SLASH_" = "/",
        "_BACKSLASH_" = "\\",
        "_PLUS_" = "+",
        "_EQUALS_" = "=",
        "_SPACE_" = " ",
        "_UNDERSCORE_" = "_",
        
        # Numeric HTML entity codes
        "__46__" = ".",   # period
        "__44__" = ",",   # comma
        "__58__" = ":",   # colon
        "__59__" = ";",   # semicolon
        "__47__" = "/",   # forward slash
        "__92__" = "\\",  # backslash
        "__43__" = "+",   # plus
        "__61__" = "=",   # equals
        "__32__" = " ",   # space
        "__95__" = "_",   # underscore
        
        # Common compartment encoding fixes
        "__91__c0__93__" = "[c]",
        "__91__e0__93__" = "[e]", 
        "__91__p0__93__" = "[p]",
        "__91__m0__93__" = "[m]",
        "__91__x0__93__" = "[x]",
        "__91__v0__93__" = "[v]",
        "__91__n0__93__" = "[n]",
        "__91__r0__93__" = "[r]",
        
        # Alternative compartment patterns
        "_LSQB_c0_RSQB_" = "[c]",
        "_LSQB_e0_RSQB_" = "[e]",
        "_LSQB_p0_RSQB_" = "[p]",
        "_LSQB_m0_RSQB_" = "[m]"
    )
    
    # Apply fixes in order of specificity (more specific patterns first)
    fixed_ids <- metabolite_ids
    
    # Sort patterns by length (longest first) to handle overlapping patterns correctly
    patterns_ordered <- names(encoding_fixes)[order(nchar(names(encoding_fixes)), decreasing = TRUE)]
    
    for (pattern in patterns_ordered) {
        replacement <- encoding_fixes[[pattern]]
        fixed_ids <- str_replace_all(fixed_ids, fixed(pattern), replacement)
    }
    
    # Additional cleanup for any remaining numeric codes
    # Look for pattern __NUMBER__ and try to convert to ASCII character
    numeric_pattern_matches <- str_extract_all(fixed_ids, "__\\d+__")
    unique_numeric_codes <- unique(unlist(numeric_pattern_matches))
    
    for (code in unique_numeric_codes) {
        if (!is.na(code)) {
            # Extract the number
            number <- as.numeric(str_extract(code, "\\d+"))
            if (!is.na(number) && number >= 32 && number <= 126) {
                # Convert to ASCII character if it's in printable range
                char <- rawToChar(as.raw(number))
                fixed_ids <- str_replace_all(fixed_ids, fixed(code), char)
            }
        }
    }
    
    return(fixed_ids)
}

#' Validate data frame row count preservation
#' @param df Data frame to check
#' @param expected_rows Expected number of rows
#' @param operation_name Name of operation for error message
validate_row_count <- function(df, expected_rows, operation_name) {
    actual_rows <- nrow(df)
    if (actual_rows != expected_rows) {
        stop("CRITICAL: Row count changed during ", operation_name, 
             " from ", expected_rows, " to ", actual_rows)
    }
    return(TRUE)
}

# =============================================================================
# SBML READING AND PREPROCESSING UTILITIES
# =============================================================================

#' Enhanced SBML preprocessing with detailed logging
#' @param file_path Path to SBML file
#' @param output_path Optional output path (creates temp file if NULL)
#' @return List with processed file path and changes log
preprocess_sbml_file_enhanced <- function(file_path, output_path = NULL) {
    
    cat("=== SBML File Preprocessing ===\n")
    cat("Input file:", basename(file_path), "\n")
    
    changes_log <- list(
        original_file = file_path,
        lines_removed = 0,
        elements_removed = list()
    )
    
    # Read file content
    sbml_content <- readLines(file_path, warn = FALSE)
    original_line_count <- length(sbml_content)
    
    # Remove problematic group elements
    cat("Removing group elements...\n")
    cleaned_content <- c()
    in_groups_section <- FALSE
    groups_depth <- 0
    removed_lines <- 0
    
    for (line in sbml_content) {
        # Detect group sections
        if (str_detect(line, "<listOfGroups|<groups:group")) {
            in_groups_section <- TRUE
            groups_depth <- str_count(line, "<[^/]") - str_count(line, "</")
            removed_lines <- removed_lines + 1
            next
        }
        
        if (in_groups_section) {
            groups_depth <- groups_depth + str_count(line, "<[^/]") - str_count(line, "</")
            removed_lines <- removed_lines + 1
            
            if (groups_depth <= 0) {
                in_groups_section <- FALSE
                groups_depth <- 0
            }
            next
        }
        
        line <- str_replace_all(line, 'groups:required="[^"]*"', '')
        
        cleaned_content <- c(cleaned_content, line)
    }
    
    changes_log$lines_removed <- removed_lines
    changes_log$elements_removed$group_lines <- removed_lines
    
    # Write output
    if (is.null(output_path)) {
        output_path <- tempfile(fileext = ".xml")
    }
    
    writeLines(cleaned_content, output_path)
    changes_log$output_file <- output_path
    
    cat("Preprocessing complete. Removed", removed_lines, "lines\n")
    cat("Output file:", output_path, "\n")
    
    return(list(file_path = output_path, changes_log = changes_log))
}

#' Attempt to read SBML with multiple fallback strategies
#' @param input_file Path to SBML file
#' @param config Configuration list
#' @return List with SBML model and read metadata
read_sbml_with_fallbacks <- function(input_file, config = list()) {
    
    cat("Reading SBML model...\n")
    sbml_params <- config$sbml_params %||% list()
    
    read_attempts <- list(
        "standard" = sbml_params,
        "safe_mode" = list(validateSBML = FALSE, bndCond = FALSE),
        "minimal" = list(validateSBML = FALSE, bndCond = FALSE, mergeMet = FALSE, balanceReact = FALSE),
        "preprocessed" = "use_preprocessing"
    )
    
    sbml_model <- NULL
    read_metadata <- list()
    
    for (attempt_name in names(read_attempts)) {
        tryCatch({
            if (attempt_name == "preprocessed") {
                cat("    Trying with file preprocessing...\n")
                preprocess_result <- preprocess_sbml_file_enhanced(input_file)
                sbml_model <- do.call(readSBMLmod, c(list(file = preprocess_result$file_path), 
                                                     list(validateSBML = FALSE, bndCond = FALSE)))
                unlink(preprocess_result$file_path)
                read_metadata$used_preprocessing <- TRUE
                read_metadata$preprocessing_changes <- preprocess_result$changes_log
            } else {
                sbml_model <- do.call(readSBMLmod, c(list(file = input_file), read_attempts[[attempt_name]]))
            }
            
            read_metadata$sbml_read_mode <- attempt_name
            if (attempt_name != "standard") {
                read_metadata$warnings <- append(read_metadata$warnings %||% list(), 
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
    
    return(list(model = sbml_model, metadata = read_metadata))
}

# =============================================================================
# COMPARTMENT PROCESSING UTILITIES
# =============================================================================

#' Standardize compartment names and mappings
#' @param sbml_model SBML model object
#' @return List with standardized compartment key and model
standardize_compartments <- function(sbml_model) {
    
    compart_key <- sbml_model@mod_compart
    names(compart_key) <- 1:length(sbml_model@mod_compart)
    
    # Standardize compartment names
    compart_key <- recode(compart_key,
                          "Cytosol" = "c", "Cytoplasm" = "c", "cytosol" = "c", 
                          "extracellular space" = "e", "extracellular" = "e", "Extra_organism" = "e",
                          "Periplasm" = "p", "periplasm" = "p",
                          "Mitochondria" = "m", "mitochondria" = "m",
                          "Nucleus" = "n", "nucleus" = "n",
                          "Vacuole" = "v", "vacuole" = "v",
                          "Endoplasmic_reticulum" = "r", "endoplasmic_reticulum" = "r",
                          "C_p" = "p", "C_m" = "m", "C_c" = "c", "C_e" = "e")
    
    # Remove trailing zeros from compartment names
    compart_key <- gsub("0$", "", compart_key)
    
    # Update model compartments
    sbml_model@mod_compart <- compart_key
    
    return(list(
        model = sbml_model,
        compartment_key = compart_key
    ))
}

# =============================================================================
# ANNOTATION PARSING UTILITIES
# =============================================================================

#' Parse RDF annotation to extract database IDs
#' @param rdf_annotation RDF annotation string
#' @return List of database IDs
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

#' Extract database IDs from string format annotations
#' @param annotation String annotation
#' @param databases Vector of database names to extract
#' @return Named list of extracted IDs
extract_string_annotation_ids <- function(annotation, databases) {
    
    result <- list()
    
    if (is.na(annotation) || annotation == "") {
        return(result)
    }
    
    if ("metanetx" %in% databases) {
        metanetx_match <- str_extract(annotation, "metanetx\\.chemical/([^;\\s]+)")
        if (!is.na(metanetx_match)) {
            result$metanetx <- str_replace(metanetx_match, "metanetx\\.chemical/", "")
        }
    }
    
    if ("bigg" %in% databases) {
        bigg_match <- str_extract(annotation, "bigg\\.metabolite/([^;\\s]+)")
        if (!is.na(bigg_match)) {
            result$bigg <- str_replace(bigg_match, "bigg\\.metabolite/", "")
        }
    }
    
    if ("seed" %in% databases) {
        seed_match <- str_extract(annotation, "seed\\.compound/([^;\\s]+)")
        if (!is.na(seed_match)) {
            result$seed <- str_replace(seed_match, "seed\\.compound/", "")
        }
    }
    
    if ("kegg" %in% databases) {
        kegg_match <- str_extract(annotation, "kegg\\.compound/([^;\\s]+)")
        if (!is.na(kegg_match)) {
            result$kegg <- str_replace(kegg_match, "kegg\\.compound/", "")
        }
    }
    
    if ("chebi" %in% databases) {
        chebi_match <- str_extract(annotation, "CHEBI:([^;\\s]+)")
        if (!is.na(chebi_match)) {
            result$chebi <- str_replace(chebi_match, "CHEBI:", "")
        }
    }
    
    return(result)
}

# =============================================================================
# DATABASE CONVERSION UTILITIES
# =============================================================================

#' Convert single database type to MetanetX with safe row preservation
#' @param met_df Metabolite data frame
#' @param database_column Column name in met_df containing database IDs
#' @param metanetx_column Output column name for MetanetX IDs
#' @param ref_data Reference data
#' @param source_filter Optional source filter for reference data
#' @return Updated metabolite data frame
convert_single_database_safe <- function(met_df, database_column, metanetx_column, ref_data, source_filter = NULL) {
    
    original_row_count <- nrow(met_df)
    
    # Get unique non-NA IDs for lookup
    if (database_column %in% names(met_df)) {
        unique_ids <- unique(met_df[[database_column]][!is.na(met_df[[database_column]]) & met_df[[database_column]] != ""])
    } else {
        met_df[[metanetx_column]] <- rep(NA_character_, nrow(met_df))
        return(met_df)
    }
    
    if (length(unique_ids) > 0) {
        # Create lookup table
        lookup_query <- ref_data$chem_xref %>%
            filter(source_id %in% unique_ids, !is.na(source_id))
        
        if (!is.null(source_filter)) {
            lookup_query <- lookup_query %>% filter(source == source_filter)
        }
        
        lookup_table <- lookup_query %>%
            distinct(source_id, .keep_all = TRUE) %>%
            select(source_id, ID)
        
        # Create named vector for fast lookup
        lookup_vector <- setNames(lookup_table$ID, lookup_table$source_id)
        
        # Apply conversion with exact length preservation
        met_df[[metanetx_column]] <- rep(NA_character_, nrow(met_df))
        valid_indices <- which(!is.na(met_df[[database_column]]) & met_df[[database_column]] != "")
        
        for (i in valid_indices) {
            db_id <- met_df[[database_column]][i]
            if (db_id %in% names(lookup_vector)) {
                met_df[[metanetx_column]][i] <- lookup_vector[db_id]
            }
        }
    } else {
        met_df[[metanetx_column]] <- rep(NA_character_, nrow(met_df))
    }
    
    # Validate row count preservation
    validate_row_count(met_df, original_row_count, paste("converting", database_column))
    
    return(met_df)
}

# =============================================================================
# MODEL VALIDATION AND UPDATING UTILITIES
# =============================================================================

#' Update exchange reactions to match new metabolite IDs
#' @param sbml_model sybilSBML model object
#' @return Updated sybilSBML model object
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
#' @param sbml_model sybilSBML model object
#' @return Updated sybilSBML model object
standardize_compartment_notation <- function(sbml_model) {
    
    sbml_model@react_id <- gsub("(e)", "_e", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(p)", "_p", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(c)", "_c", sbml_model@react_id, fixed = TRUE)
    sbml_model@react_id <- gsub("(m)", "_m", sbml_model@react_id, fixed = TRUE)
    
    return(sbml_model)
}

#' Clean GPR (Gene-Protein-Reaction) associations
#' @param sbml_model sybilSBML model object
#' @return Updated sybilSBML model object
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


# =============================================================================
# LOGGING AND REPORTING UTILITIES
# =============================================================================

#' Log metabolite conversions for detailed tracking
#' @param met_df Metabolite data frame with conversion information
#' @param model_id Model identifier
#' @param output_dir Directory to save conversion logs
#' @return List with conversion statistics and log file path
log_metabolite_conversions <- function(met_df, model_id, output_dir) {
    
    # Create detailed conversion log
    conversion_log <- met_df %>%
        select(orig_met, without_compartment, new_met, new_met_out, compartment) %>%
        mutate(
            conversion_successful = !is.na(new_met) & new_met != without_compartment,
            conversion_type = case_when(
                is.na(new_met) ~ "No conversion",
                new_met == without_compartment ~ "Unchanged",
                str_detect(new_met, "^MNXM\\d+$") ~ "MetanetX conversion",
                TRUE ~ "Other conversion"
            ),
            encoding_issues_fixed = orig_met != fix_metabolite_encoding(orig_met)
        )
    
    # Add database source information if available
    database_columns <- c("metanetx", "bigg", "seed", "kegg", "chebi", "hmdb", "pubchem")
    available_db_cols <- intersect(database_columns, names(met_df))
    
    if (length(available_db_cols) > 0) {
        conversion_log$database_sources <- apply(met_df[available_db_cols], 1, function(row) {
            sources <- names(row)[!is.na(row) & row != ""]
            if (length(sources) == 0) return("None")
            return(paste(sources, collapse = ", "))
        })
    }
    
    # Calculate conversion statistics
    conversion_stats <- list(
        model_id = model_id,
        total_metabolites = nrow(conversion_log),
        successful_conversions = sum(conversion_log$conversion_successful, na.rm = TRUE),
        metanetx_conversions = sum(conversion_log$conversion_type == "MetanetX conversion", na.rm = TRUE),
        unchanged_metabolites = sum(conversion_log$conversion_type == "Unchanged", na.rm = TRUE),
        failed_conversions = sum(conversion_log$conversion_type == "No conversion", na.rm = TRUE),
        encoding_fixes_applied = sum(conversion_log$encoding_issues_fixed, na.rm = TRUE),
        conversion_rate = round(sum(conversion_log$conversion_successful, na.rm = TRUE) / nrow(conversion_log) * 100, 1),
        metanetx_rate = round(sum(conversion_log$conversion_type == "MetanetX conversion", na.rm = TRUE) / nrow(conversion_log) * 100, 1)
    )
    
    # Create summary by conversion type
    conversion_summary <- conversion_log %>%
        group_by(conversion_type) %>%
        summarise(
            count = n(),
            percentage = round(n() / nrow(conversion_log) * 100, 1),
            .groups = "drop"
        )
    
    # Create summary by compartment
    compartment_summary <- conversion_log %>%
        group_by(compartment, conversion_type) %>%
        summarise(count = n(), .groups = "drop") %>%
        pivot_wider(names_from = conversion_type, values_from = count, values_fill = 0)
    
    # Save detailed conversion log
    log_file_path <- file.path(output_dir, paste0(model_id, "_metabolite_conversions.csv"))
    write_csv(conversion_log, log_file_path)
    
    # Save conversion summary
    summary_file_path <- file.path(output_dir, paste0(model_id, "_conversion_summary.json"))
    summary_data <- list(
        statistics = conversion_stats,
        by_type = conversion_summary,
        by_compartment = compartment_summary
    )
    write_json(summary_data, summary_file_path, pretty = TRUE, auto_unbox = TRUE)
    
    cat("  📊 Conversion log saved:", log_file_path, "\n")
    cat("  📈 Conversion rate:", conversion_stats$conversion_rate, "% (", 
        conversion_stats$successful_conversions, "/", conversion_stats$total_metabolites, ")\n")
    cat("  🔧 MetanetX conversions:", conversion_stats$metanetx_rate, "% (", 
        conversion_stats$metanetx_conversions, ")\n")
    if (conversion_stats$encoding_fixes_applied > 0) {
        cat("  🔧 Encoding issues fixed:", conversion_stats$encoding_fixes_applied, "metabolites\n")
    }
    
    return(list(
        statistics = conversion_stats,
        log_file = log_file_path,
        summary_file = summary_file_path
    ))
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Merge two lists
#' @param list1 First list
#' @param list2 Second list
#' @return Merged list
merge_lists <- function(list1, list2) {
    for (name in names(list2)) {
        list1[[name]] <- list2[[name]]
    }
    return(list1)
}
