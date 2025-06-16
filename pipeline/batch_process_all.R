# batch_processing.R
# Batch processing functions for remote and local operations

library(here)
library(stringr)
library(jsonlite)

# Source required functions
source("process_sbml_species.R")      # Main processing pipeline
source("validate_model_growth.R")     # Validation functions

#' Batch process all species on remote machine (where sybilSBML works)
#' @param base_dir Path to microbial_gem_database/species directory
#' @param ref_data Reference data from get_reference_data()
#' @param deprecated_recode Deprecated ID mappings
#' @param species_filter Optional vector of species to process (for testing)
#' @return Summary of processing results
#' 
batch_process_remote <- function(base_dir, ref_data, deprecated_recode, species_filter = NULL) {
cat("=== BATCH PROCESSING (REMOTE) ===\n")
cat("Processing SBML files with sybilSBML...\n\n")

# Handle multiple directories
if (length(base_dir) > 1) {
    # Process multiple directories
    all_species_dirs <- c()
    all_species_names <- c()
    
    for (dir in base_dir) {
        if (dir.exists(dir)) {
            species_dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
            species_names <- basename(species_dirs)
            all_species_dirs <- c(all_species_dirs, species_dirs)
            all_species_names <- c(all_species_names, species_names)
            cat("Found", length(species_dirs), "species in", dir, "\n")
        } else {
            cat("Warning: Directory not found:", dir, "\n")
        }
    }
    species_dirs <- all_species_dirs
    species_names <- all_species_names
} else {
    # Single directory (original behavior)
    species_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    species_names <- basename(species_dirs)
}

if (!is.null(species_filter)) {
    species_dirs <- species_dirs[species_names %in% species_filter]
    species_names <- species_names[species_names %in% species_filter]
}

cat("Found", length(species_dirs), "species directories total\n")

# Initialize results tracking with proper defaults
results_summary <- data.frame(
    species = species_names,
    processing_success = FALSE,
    processing_time = NA_real_,
    pattern_detected = NA_character_,
    metabolites_total = NA_integer_,
    metabolites_converted = NA_integer_,
    conversion_rate = NA_real_,
    problems_detected = NA_integer_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
)

start_time <- Sys.time()

# Process each species
for (i in seq_along(species_dirs)) {
    species_dir <- species_dirs[i]
    species_name <- species_names[i]
    
    cat("\n--- Processing", i, "of", length(species_dirs), ":", species_name, "---\n")
    
    tryCatch({
        result <- process_single_species(species_dir, ref_data, deprecated_recode)
        
        if (result$success) {
            results_summary[i, "processing_success"] <- TRUE
            
            # Safe assignment with null checking
            if (!is.null(result$processing_log$processing_time_seconds)) {
                results_summary[i, "processing_time"] <- result$processing_log$processing_time_seconds
            }
            if (!is.null(result$processing_log$pattern_detected)) {
                results_summary[i, "pattern_detected"] <- result$processing_log$pattern_detected
            }
            if (!is.null(result$processing_log$metabolites_total)) {
                results_summary[i, "metabolites_total"] <- result$processing_log$metabolites_total
            }
            if (!is.null(result$processing_log$metabolites_standardized)) {
                results_summary[i, "metabolites_converted"] <- result$processing_log$metabolites_standardized
            }
            if (!is.null(result$processing_log$conversion_summary$conversion_rate)) {
                results_summary[i, "conversion_rate"] <- result$processing_log$conversion_summary$conversion_rate
            }
            if (!is.null(result$processing_log$problems_detected)) {
                results_summary[i, "problems_detected"] <- length(result$processing_log$problems_detected)
            }
            
            cat("✓ Success! Conversion rate:", 
                ifelse(is.null(result$processing_log$conversion_summary$conversion_rate), "Unknown", 
                       paste0(result$processing_log$conversion_summary$conversion_rate, "%")), "\n")
        } else {
            if (!is.null(result$error)) {
                results_summary[i, "error_message"] <- result$error
            }
            cat("✗ Failed:", ifelse(is.null(result$error), "Unknown error", result$error), "\n")
        }
        
    }, error = function(e) {
        results_summary[i, "error_message"] <- e$message
        cat("✗ Unexpected error:", e$message, "\n")
    })
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Summary statistics with safe calculations
cat("\n=== BATCH PROCESSING SUMMARY ===\n")
cat("Total time:", round(total_time, 2), "minutes\n")
cat("Species processed:", nrow(results_summary), "\n")
cat("Successful:", sum(results_summary$processing_success, na.rm = TRUE), "\n")
cat("Failed:", sum(!results_summary$processing_success, na.rm = TRUE), "\n")

if (any(results_summary$processing_success, na.rm = TRUE)) {
    valid_rates <- results_summary$conversion_rate[!is.na(results_summary$conversion_rate)]
    valid_times <- results_summary$processing_time[!is.na(results_summary$processing_time)]
    
    if (length(valid_rates) > 0) {
        cat("Average conversion rate:", round(mean(valid_rates), 1), "%\n")
    }
    if (length(valid_times) > 0) {
        cat("Average processing time:", round(mean(valid_times), 2), "seconds\n")
    }
}

# Show failed species
failed_species <- results_summary[!results_summary$processing_success, ]
if (nrow(failed_species) > 0) {
    cat("\nFailed species:\n")
    for (i in 1:nrow(failed_species)) {
        error_msg <- failed_species[i, "error_message"]
        if (is.na(error_msg)) error_msg <- "Unknown error"
        cat("  -", failed_species[i, "species"], ":", error_msg, "\n")
    }
}

# Write summary report with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_file <- file.path(dirname(base_dir[1]), "reports", paste0("batch_processing_summary_", timestamp, ".csv"))
dir.create(dirname(summary_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results_summary, summary_file, row.names = FALSE)
cat("\nDetailed results saved to:", summary_file, "\n")

return(results_summary)
}

# Update validation function to handle multiple directories
run_post_processing_validation <- function(species_base_dirs = "species", 
                                           python_cmd = NULL, 
                                           force_revalidate = FALSE) {
    
    cat("=== Post-Processing Validation ===\n")
    
    # Handle multiple directories
    if (length(species_base_dirs) > 1) {
        all_species_dirs <- c()
        for (dir in species_base_dirs) {
            if (dir.exists(dir)) {
                species_dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
                species_dirs <- species_dirs[!str_detect(basename(species_dirs), "^\\.|temp|backup")]
                all_species_dirs <- c(all_species_dirs, species_dirs)
                cat("Found", length(species_dirs), "species in", dir, "\n")
            }
        }
        species_dirs <- all_species_dirs
    } else {
        # Single directory
        if (!dir.exists(species_base_dirs)) {
            stop("Species directory not found: ", species_base_dirs)
        }
        species_dirs <- list.dirs(species_base_dirs, recursive = FALSE, full.names = TRUE)
        species_dirs <- species_dirs[!str_detect(basename(species_dirs), "^\\.|temp|backup")]
        cat("Base directory:", species_base_dirs, "\n")
    }
    
    # Rest of function remains the same...
    cat("Found", length(species_dirs), "species directories total\n")
    
    # Filter to species that have processed files
    processed_species <- c()
    for (species_dir in species_dirs) {
        processed_files <- list.files(species_dir, pattern = "_processed\\.xml$", full.names = TRUE)
        if (length(processed_files) > 0) {
            processed_species <- c(processed_species, species_dir)
        }
    }
    
    cat("Found", length(processed_species), "species with processed files\n")
    
    if (length(processed_species) == 0) {
        cat("No processed species found. Run batch processing first.\n")
        return(NULL)
    }
    
    # Initialize validation summary
    validation_summary <- list(
        total_species = length(processed_species),
        validation_attempted = 0,
        validation_successful = 0,
        validation_failed = 0,
        cobra_readable_input = 0,
        cobra_readable_processed = 0,
        simulation_equivalent = 0,
        results = list()
    )
    
    # Load validation function
    if (!exists("validate_model_growth")) {
        tryCatch({
            source(here("validate_model_growth.R"))
        }, error = function(e) {
            cat("Warning: Could not load validate_model_growth.R\n")
            cat("Make sure validate_model_growth.R is in your working directory\n")
            return(NULL)
        })
    }
    
    # Process each species
    for (i in seq_along(processed_species)) {
        species_dir <- processed_species[i]
        species_name <- basename(species_dir)
        
        cat("\n--- Validating", species_name, "(", i, "of", length(processed_species), ") ---\n")
        
        # Check if validation already exists
        validation_file <- file.path(species_dir, "validation_results.json")
        if (file.exists(validation_file) && !force_revalidate) {
            cat("Validation results already exist, skipping (use force_revalidate=TRUE to override)\n")
            
            # Load existing results for summary
            tryCatch({
                existing_results <- fromJSON(validation_file)
                validation_summary$results[[species_name]] <- existing_results
                
                if (existing_results$success) {
                    validation_summary$validation_successful <- validation_summary$validation_successful + 1
                    if (existing_results$input_readable_cobra) {
                        validation_summary$cobra_readable_input <- validation_summary$cobra_readable_input + 1
                    }
                    if (existing_results$processed_readable_cobra) {
                        validation_summary$cobra_readable_processed <- validation_summary$cobra_readable_processed + 1
                    }
                    if (!is.null(existing_results$simulation_equivalent) && existing_results$simulation_equivalent) {
                        validation_summary$simulation_equivalent <- validation_summary$simulation_equivalent + 1
                    }
                } else {
                    validation_summary$validation_failed <- validation_summary$validation_failed + 1
                }
            }, error = function(e) {
                cat("Warning: Could not read existing validation file\n")
            })
            
            next
        }
        
        # Run validation
        validation_summary$validation_attempted <- validation_summary$validation_attempted + 1
        
        tryCatch({
            validation_result <- validate_model_growth(species_dir, python_cmd)
            
            # Store result
            validation_summary$results[[species_name]] <- validation_result
            
            # Update counters
            if (validation_result$success) {
                validation_summary$validation_successful <- validation_summary$validation_successful + 1
                
                if (validation_result$input_readable_cobra) {
                    validation_summary$cobra_readable_input <- validation_summary$cobra_readable_input + 1
                }
                
                if (validation_result$processed_readable_cobra) {
                    validation_summary$cobra_readable_processed <- validation_summary$cobra_readable_processed + 1
                }
                
                if (!is.null(validation_result$simulation_equivalent) && validation_result$simulation_equivalent) {
                    validation_summary$simulation_equivalent <- validation_summary$simulation_equivalent + 1
                }
                
                # Print summary
                cat("✓ Validation successful\n")
                if (!is.null(validation_result$growth_rate_input)) {
                    cat("  Input growth rate:", round(validation_result$growth_rate_input, 6), "\n")
                }
                if (!is.null(validation_result$growth_rate_processed)) {
                    cat("  Processed growth rate:", round(validation_result$growth_rate_processed, 6), "\n")
                }
                if (!is.null(validation_result$simulation_equivalent)) {
                    cat("  Simulation equivalent:", validation_result$simulation_equivalent, "\n")
                }
                
            } else {
                validation_summary$validation_failed <- validation_summary$validation_failed + 1
                cat("✗ Validation failed:", validation_result$error, "\n")
            }
            
        }, error = function(e) {
            validation_summary$validation_failed <- validation_summary$validation_failed + 1
            cat("✗ Validation error:", e$message, "\n")
            
            # Store error result
            validation_summary$results[[species_name]] <- list(
                success = FALSE,
                error = e$message,
                model_id = species_name
            )
        })
        
        # Progress indicator
        if (i %% 5 == 0 || i == length(processed_species)) {
            cat("\nProgress:", i, "of", length(processed_species), "completed\n")
        }
    }
    
    # Print final summary
    cat("\n=== Validation Summary ===\n")
    cat("Total species processed:", validation_summary$total_species, "\n")
    cat("Validation attempted:", validation_summary$validation_attempted, "\n")
    cat("Validation successful:", validation_summary$validation_successful, "\n")
    cat("Validation failed:", validation_summary$validation_failed, "\n")
    
    if (validation_summary$validation_successful > 0) {
        cat("\nCOBRApy Compatibility:\n")
        cat("Input files readable:", validation_summary$cobra_readable_input, "/", validation_summary$validation_successful, "\n")
        cat("Processed files readable:", validation_summary$cobra_readable_processed, "/", validation_summary$validation_successful, "\n")
        cat("Simulation equivalent:", validation_summary$simulation_equivalent, "/", validation_summary$validation_successful, "\n")
        
        # Calculate percentages
        if (validation_summary$validation_successful > 0) {
            input_pct <- round(validation_summary$cobra_readable_input / validation_summary$validation_successful * 100, 1)
            processed_pct <- round(validation_summary$cobra_readable_processed / validation_summary$validation_successful * 100, 1)
            equiv_pct <- round(validation_summary$simulation_equivalent / validation_summary$validation_successful * 100, 1)
            
            cat("Percentages: Input", input_pct, "%, Processed", processed_pct, "%, Equivalent", equiv_pct, "%\n")
        }
    }
    
    # Save summary
    base_dir_name <- ifelse(length(species_base_dirs) > 1, "multi_dir", basename(species_base_dirs[1]))
    summary_file <- file.path(base_dir_name, "validation_summary.json")
    write_json(validation_summary, summary_file, pretty = TRUE, auto_unbox = TRUE)
    cat("\nValidation summary saved to:", summary_file, "\n")
    
    return(validation_summary)
}


#' Create a validation report from existing validation results
#' @param species_base_dir Base directory containing species folders
#' @param output_file Output file for the report (default: "validation_report.txt")
create_validation_report <- function(species_base_dir = "species", 
                                     output_file = "validation_report.txt") {
    
    cat("Creating validation report...\n")
    
    # Find all validation results
    species_dirs <- list.dirs(species_base_dir, recursive = FALSE, full.names = TRUE)
    validation_files <- file.path(species_dirs, "validation_results.json")
    validation_files <- validation_files[file.exists(validation_files)]
    
    if (length(validation_files) == 0) {
        cat("No validation results found\n")
        return(NULL)
    }
    
    # Collect all results
    all_results <- list()
    for (val_file in validation_files) {
        species_name <- basename(dirname(val_file))
        tryCatch({
            result <- fromJSON(val_file)
            all_results[[species_name]] <- result
        }, error = function(e) {
            cat("Warning: Could not read", val_file, "\n")
        })
    }
    
    # Generate report
    report_lines <- c(
        "# SBML Processing and Validation Report",
        paste("Generated:", Sys.time()),
        paste("Total species:", length(all_results)),
        "",
        "## Summary Statistics"
    )
    
    # Calculate statistics with proper logical handling
    successful <- sum(sapply(all_results, function(x) {
        success_val <- x$success
        if (is.null(success_val)) return(FALSE)
        if (is.logical(success_val)) return(success_val)
        return(FALSE)
    }))
    
    input_readable <- sum(sapply(all_results, function(x) {
        readable_val <- x$input_readable_cobra
        if (is.null(readable_val)) return(FALSE)
        if (is.logical(readable_val)) return(readable_val)
        return(FALSE)
    }))
    
    processed_readable <- sum(sapply(all_results, function(x) {
        readable_val <- x$processed_readable_cobra
        if (is.null(readable_val)) return(FALSE)
        if (is.logical(readable_val)) return(readable_val)
        return(FALSE)
    }))
    
    equivalent <- sum(sapply(all_results, function(x) {
        equiv_val <- x$simulation_equivalent
        if (is.null(equiv_val)) return(FALSE)
        if (is.logical(equiv_val)) return(equiv_val)
        return(FALSE)
    }))
    
    report_lines <- c(report_lines,
                      paste("- Validation successful:", successful, "/", length(all_results)),
                      paste("- Input files COBRA-readable:", input_readable, "/", successful),
                      paste("- Processed files COBRA-readable:", processed_readable, "/", successful),
                      paste("- Simulation equivalent:", equivalent, "/", successful),
                      ""
    )
    
    # Individual species results
    report_lines <- c(report_lines, "## Individual Results")
    
    for (species_name in names(all_results)) {
        result <- all_results[[species_name]]
        
        # Safe logical extraction
        success_val <- result$success
        if (is.null(success_val) || !is.logical(success_val)) success_val <- FALSE
        
        status <- if (success_val) "✓" else "✗"
        line <- paste(status, species_name)
        
        if (success_val) {
            details <- c()
            
            input_readable_val <- result$input_readable_cobra
            if (!is.null(input_readable_val) && is.logical(input_readable_val) && input_readable_val) {
                details <- c(details, "input-readable")
            }
            
            processed_readable_val <- result$processed_readable_cobra
            if (!is.null(processed_readable_val) && is.logical(processed_readable_val) && processed_readable_val) {
                details <- c(details, "processed-readable")
            }
            
            equiv_val <- result$simulation_equivalent
            if (!is.null(equiv_val) && is.logical(equiv_val) && equiv_val) {
                details <- c(details, "equivalent")
            }
            
            if (length(details) > 0) {
                line <- paste(line, "(", paste(details, collapse = ", "), ")")
            }
            
            # Add growth rates if available
            if (!is.null(result$growth_rate_input) && !is.null(result$growth_rate_processed) &&
                is.numeric(result$growth_rate_input) && is.numeric(result$growth_rate_processed)) {
                growth_info <- paste("Growth:", round(result$growth_rate_input, 4), "→", round(result$growth_rate_processed, 4))
                line <- paste(line, "-", growth_info)
            }
        } else {
            error_msg <- result$error
            if (is.null(error_msg)) error_msg <- "Unknown error"
            line <- paste(line, "-", str_trunc(error_msg, 60))
        }
        
        report_lines <- c(report_lines, line)
    }
    
    # Problems section
    failed_species <- names(all_results)[sapply(all_results, function(x) {
        success_val <- x$success
        if (is.null(success_val) || !is.logical(success_val)) return(TRUE)
        return(!success_val)
    })]
    
    if (length(failed_species) > 0) {
        report_lines <- c(report_lines, "", "## Failed Validations")
        for (species in failed_species) {
            result <- all_results[[species]]
            error_msg <- result$error
            if (is.null(error_msg)) error_msg <- "Unknown error"
            report_lines <- c(report_lines, paste("-", species, ":", error_msg))
        }
    }
    
    # Non-equivalent species
    non_equiv <- names(all_results)[sapply(all_results, function(x) {
        success_val <- x$success
        equiv_val <- x$simulation_equivalent
        
        if (is.null(success_val) || !is.logical(success_val) || !success_val) return(FALSE)
        if (is.null(equiv_val) || !is.logical(equiv_val)) return(TRUE)
        return(!equiv_val)
    })]
    
    if (length(non_equiv) > 0) {
        report_lines <- c(report_lines, "", "## Non-Equivalent Growth Rates")
        for (species in non_equiv) {
            result <- all_results[[species]]
            if (!is.null(result$growth_rate_input) && !is.null(result$growth_rate_processed) &&
                is.numeric(result$growth_rate_input) && is.numeric(result$growth_rate_processed)) {
                diff <- abs(result$growth_rate_input - result$growth_rate_processed)
                report_lines <- c(report_lines, 
                                  paste("-", species, ": input =", round(result$growth_rate_input, 6), 
                                        ", processed =", round(result$growth_rate_processed, 6), 
                                        ", difference =", round(diff, 8)))
            } else {
                report_lines <- c(report_lines, paste("-", species, ": growth rate comparison failed"))
            }
        }
    }
    
    # Write report
    writeLines(report_lines, output_file)
    cat("Validation report written to:", output_file, "\n")
    
    return(list(
        total_species = length(all_results),
        successful = successful,
        input_readable = input_readable,
        processed_readable = processed_readable,
        equivalent = equivalent,
        failed_species = failed_species,
        non_equivalent = non_equiv
    ))
}

#' Quick validation status check
#' @param species_base_dir Base directory containing species folders
check_validation_status <- function(species_base_dir = "species") {
  
  # Find all species directories
  species_dirs <- list.dirs(species_base_dir, recursive = FALSE, full.names = TRUE)
  species_dirs <- species_dirs[!str_detect(basename(species_dirs), "^\\.|temp|backup")]
  
  # Check processing and validation status
  status_summary <- data.frame(
    species = basename(species_dirs),
    has_processed = FALSE,
    has_validation = FALSE,
    validation_success = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(species_dirs)) {
    species_dir <- species_dirs[i]
    
    # Check for processed files
    processed_files <- list.files(species_dir, pattern = "_processed\\.xml$")
    status_summary$has_processed[i] <- length(processed_files) > 0
    
    # Check for validation results
    validation_file <- file.path(species_dir, "validation_results.json")
    status_summary$has_validation[i] <- file.exists(validation_file)
    
    if (status_summary$has_validation[i]) {
      tryCatch({
        validation_result <- fromJSON(validation_file)
        status_summary$validation_success[i] <- validation_result$success %||% FALSE
      }, error = function(e) {
        status_summary$validation_success[i] <- FALSE
      })
    }
  }
  
  # Print summary
  cat("=== Validation Status Summary ===\n")
  cat("Total species:", nrow(status_summary), "\n")
  cat("With processed files:", sum(status_summary$has_processed), "\n")
  cat("With validation results:", sum(status_summary$has_validation), "\n")
  cat("Successful validations:", sum(status_summary$validation_success), "\n")
  
  # Show species needing validation
  needs_validation <- status_summary[status_summary$has_processed & !status_summary$has_validation, ]
  if (nrow(needs_validation) > 0) {
    cat("\nSpecies needing validation:\n")
    for (species in needs_validation$species) {
      cat("  -", species, "\n")
    }
  }
  
  # Show failed validations
  failed_validation <- status_summary[status_summary$has_validation & !status_summary$validation_success, ]
  if (nrow(failed_validation) > 0) {
    cat("\nSpecies with failed validation:\n")
    for (species in failed_validation$species) {
      cat("  -", species, "\n")
    }
  }
  
  return(status_summary)
}


# Usage for remote processing to convert to MetanetX namespace:
# 
# # Load reference data
# if (!exists("ref_data")) ref_data <- get_reference_data()
# if (!exists("deprecated_recode")) deprecated_recode <- readRDS("reference_data/deprecated_recode_mets.rds")
# 
# # Process all species 
# results <- batch_process_remote("microbial_gem_database/species", ref_data, deprecated_recode)
# 
# # Or test with just a few species:
# test_species <- c("nitrosomonas_europaea_iGC535", "azotobacter_vinelandii_iAA1300")
# results <- batch_process_remote("microbial_gem_database/species", ref_data, deprecated_recode, test_species)

library(here)
library(stringr)
library(tidyr)
library(tidyverse)
library(yaml)
library(jsonlite)
library(xml2)  # For RDF parsing


source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/pipeline/process_sbml_species.R")
source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/pipeline/sbml_processing_utils.R")
source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/pipeline/processing_utils.R")

ref_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/metanetx_reference_data.rds")
deprecated_recode <- readRDS("/projectnb/talbot-lab-data/zrwerbin/microbial_gem_database/reference_data/deprecated_recode_mets.rds")

# Load reference data
if (!exists("ref_data")) ref_data <- get_reference_data()

# ran successfully
test_species <- c("azotobacter_vinelandii_iAA1300", 
                  "bacillus_subtilis_iBB1018", "bradyrhizobium_diazoefficiens_iYY1101", 
                  "clostridium_ljungdahlii_iHN637", "ensifer_meliloti_iGD1348", 
                  "lachancea_thermotolerans_iBM3063", 
                  "methanosarcina_barkeri_iMG746", "mortierella_alpina_iCY1106", 
                  "nitrobacter_winogradskyi_iFC579", "nitrosomonas_europaea_iGC535", 
                  "nitrospira_moscoviensis_iNmo686", "pseudomonas_putida_iJN1462", 
                  "rhizophagus_irregularis_iRi1574", "saccharomyces_cerevisiae_iMM904", 
                  "staphylococcus_aureus_iSB619", "streptomyces_coelicolor_iKS1317")

failed_species = c("agrobacterium_tumefaciens_iNX1344",#"aspergillus_terreus_iJL1454",
                   #"bacillus_pseudofirmus_Xu_bacill",
                   "methanosarcina_barkeri_iAF692")
                   #"rhodococcus_jostii_iMT1174","sphingopyxis_granuli_iIG743")
    
test_species = "pseudomonas_putida_iJN1462"

caused_failure = "rhodopseudomonas_palustris_iDT1294"

results <- batch_process_remote("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species", 
                                ref_data, deprecated_recode, test_species)

results <- batch_process_remote("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species", 
                               ref_data, deprecated_recode, to_finish)


results <- batch_process_remote("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species", 
                                ref_data, deprecated_recode, "rhizobium_leguminosarum_iCS1224")
# Running on all avail files
results <- batch_process_remote("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species", ref_data, deprecated_recode)

/projectnb/talbot-lab-data/metabolic_models/curated_models/CarveFungi
# USES PYTHON

# Run validation on all processed species
validation_results <- run_post_processing_validation("/Users/zoeywerbin/soil_GEM_database/microbial_gem_database/species/", force_revalidate = T)
#
# # Check current status
status <- check_validation_status("/Users/zoeywerbin/soil_GEM_database/microbial_gem_database/species/")

# Create a readable report
report_summary <- create_validation_report()

