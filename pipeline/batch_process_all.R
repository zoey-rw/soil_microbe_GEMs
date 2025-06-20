# batch_processing.R
# Batch processing functions for remote and local operations

library(here)
library(stringr)
library(jsonlite)

# Source required functions
source("process_sbml_species.R")      # Main processing pipeline

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





#' Process SBML files from a flat directory structure (like carvefungi)
#' @param input_dir Directory containing SBML files directly
#' @param output_dir Directory to save processed files
#' @param ref_data Reference data from get_reference_data()
#' @param deprecated_recode Deprecated ID mappings
#' @param file_filter Optional pattern to filter files
#' @return Summary of processing results
process_flat_directory <- function(input_dir, output_dir, ref_data, deprecated_recode, file_filter = NULL) {
    
    cat("=== PROCESSING FLAT DIRECTORY ===\n")
    cat("Input directory:", input_dir, "\n")
    cat("Output directory:", output_dir, "\n")
    
    # Create output directory if needed
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Find all SBML files
    xml_files <- list.files(input_dir, pattern = "\\.xml$", full.names = TRUE)
    sbml_files <- list.files(input_dir, pattern = "\\.sbml$", full.names = TRUE)
    all_files <- c(xml_files, sbml_files)
    
    # Apply filter if provided
    if (!is.null(file_filter)) {
        all_files <- all_files[str_detect(basename(all_files), file_filter)]
    }
    
    # Remove already processed files
    processed_patterns <- c("_processed", "_cobra_validated", "_modified_cobra")
    input_files <- all_files[!str_detect(basename(all_files), paste(processed_patterns, collapse = "|"))]
    
    cat("Found", length(input_files), "input files to process\n")
    
    if (length(input_files) == 0) {
        cat("No files to process\n")
        return(NULL)
    }
    
    # Initialize results tracking
    results_summary <- data.frame(
        file_name = basename(input_files),
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
    
    # Process each file
    for (i in seq_along(input_files)) {
        input_file <- input_files[i]
        file_name <- basename(input_file)
        model_id <- str_replace(file_name, "\\.(xml|sbml)$", "")
        
        cat("\n--- Processing", i, "of", length(input_files), ":", file_name, "---\n")
        
        # Create temporary species directory structure
        temp_species_dir <- tempfile(pattern = "temp_species_")
        dir.create(temp_species_dir)
        
        # Copy file to temp directory
        temp_input_file <- file.path(temp_species_dir, file_name)
        file.copy(input_file, temp_input_file)
        
        tryCatch({
            # Use existing process_single_species function
            result <- process_single_species(temp_species_dir, ref_data, deprecated_recode)
            
            if (result$success) {
                results_summary[i, "processing_success"] <- TRUE
                
                # Copy processed file to output directory
                processed_files <- list.files(temp_species_dir, pattern = "_processed\\.xml$", full.names = TRUE)
                if (length(processed_files) > 0) {
                    output_file <- file.path(output_dir, paste0(model_id, "_processed.xml"))
                    file.copy(processed_files[1], output_file, overwrite = TRUE)
                    cat("Saved processed file to:", output_file, "\n")
                }
                
                # Copy conversion logs if they exist
                log_dir <- file.path(temp_species_dir, "conversion_logs")
                if (dir.exists(log_dir)) {
                    output_log_dir <- file.path(output_dir, "conversion_logs")
                    if (!dir.exists(output_log_dir)) {
                        dir.create(output_log_dir, recursive = TRUE)
                    }
                    log_files <- list.files(log_dir, full.names = TRUE)
                    for (log_file in log_files) {
                        file.copy(log_file, file.path(output_log_dir, basename(log_file)), overwrite = TRUE)
                    }
                }
                
                # Store results
                if (!is.null(result$processing_log$processing_time_seconds)) {
                    results_summary[i, "processing_time"] <- result$processing_log$processing_time_seconds
                }
                if (!is.null(result$processing_log$pattern_detected)) {
                    results_summary[i, "pattern_detected"] <- result$processing_log$pattern_detected
                }
                if (!is.null(result$processing_log$metabolites_total)) {
                    results_summary[i, "metabolites_total"] <- result$processing_log$metabolites_total
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
        }, finally = {
            # Clean up temp directory
            unlink(temp_species_dir, recursive = TRUE)
        })
    }
    
    end_time <- Sys.time()
    total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    
    # Summary statistics
    cat("\n=== FLAT DIRECTORY PROCESSING SUMMARY ===\n")
    cat("Total time:", round(total_time, 2), "minutes\n")
    cat("Files processed:", nrow(results_summary), "\n")
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
    
    # Show failed files
    failed_files <- results_summary[!results_summary$processing_success, ]
    if (nrow(failed_files) > 0) {
        cat("\nFailed files:\n")
        for (i in 1:min(10, nrow(failed_files))) {
            error_msg <- failed_files[i, "error_message"]
            if (is.na(error_msg)) error_msg <- "Unknown error"
            cat("  -", failed_files[i, "file_name"], ":", error_msg, "\n")
        }
    }
    
    # Save summary report
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    summary_file <- file.path(output_dir, paste0("processing_summary_", timestamp, ".csv"))
    write.csv(results_summary, summary_file, row.names = FALSE)
    cat("\nDetailed results saved to:", summary_file, "\n")
    
    return(results_summary)
}


# Usage for remote processing to convert to MetanetX namespace:

# # Load reference data
# if (!exists("ref_data")) ref_data <- get_reference_data()
# if (!exists("deprecated_recode")) deprecated_recode <- readRDS("reference_data/deprecated_recode_mets.rds")

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


test_results <- process_flat_directory(
    input_dir = "/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/carvefungi_species/input",
    output_dir = "/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/carvefungi_species/processed",
    ref_data = ref_data,
    deprecated_recode = deprecated_recode,
    file_filter = "^[C]"  # Only files starting with A, B, or C
)
