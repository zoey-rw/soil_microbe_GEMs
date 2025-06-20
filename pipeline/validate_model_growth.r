# Functions to validate that models grow appropriately after standardizing formats and converting identifiers to the MetaNetX namespace
# Growth rate should be the same as the input model

library(stringr)
library(jsonlite)

# Source utility functions
source("sbml_processing_utils.R")

#' Validate SBML model growth in COBRApy and update validation metadata
#' @param species_dir Path to species directory
#' @param python_cmd Python command to use (e.g., "python3", "python", "/path/to/python")
#' @return List with validation results
validate_model_growth <- function(species_dir, python_cmd = NULL) {
    
    # Extract model_id more safely
    species_name <- basename(normalizePath(species_dir, mustWork = FALSE))
    if (is.null(species_name) || length(species_name) == 0 || species_name == "" || is.na(species_name)) {
        stop("Invalid species directory path: ", species_dir)
    }
    
    # More robust model_id extraction
    model_id_match <- str_extract(species_name, "[^_]+$")
    if (is.null(model_id_match) || length(model_id_match) == 0 || is.na(model_id_match) || model_id_match == "") {
        model_id <- species_name  # Fallback to full species name
    } else {
        model_id <- model_id_match
    }
    
    cat("Validating species:", species_name, "with model_id:", model_id, "\n")
    
    # Find input and processed files
    files <- tryCatch({
        discover_sbml_files(species_dir, model_id)
    }, error = function(e) {
        stop("File discovery failed for ", species_name, ": ", e$message)
    })
    input_file <- files$input_file
    
    # Check if processed files exist
    if (length(files$processed_files) == 0) {
        validation_results <- list(
            model_id = model_id,
            validation_timestamp = Sys.time(),
            input_file = basename(input_file),
            processed_file = NA,
            success = FALSE,
            error = "No processed SBML files found in species directory"
        )
        validation_file <- file.path(species_dir, "validation_results.json")
        write_json(validation_results, validation_file, pretty = TRUE, auto_unbox = TRUE)
        return(validation_results)
    }
    
    processed_file <- files$processed_files[1]  # Take first processed file
    
    validation_results <- list(
        model_id = model_id,
        validation_timestamp = Sys.time(),
        input_file = basename(input_file),
        processed_file = basename(processed_file),
        input_readable_cobra = FALSE,
        processed_readable_cobra = FALSE,
        growth_rate_input = NULL,
        growth_rate_processed = NULL,
        simulation_equivalent = NULL,
        validation_notes = list(),
        python_setup = list()
    )
    
    # Try to find working Python installation
    if (is.null(python_cmd)) {
        python_candidates <- c("python3", "python", "/usr/bin/python3", "/usr/local/bin/python3")
        
        # Check for conda/virtual environments
        if (Sys.getenv("CONDA_DEFAULT_ENV") != "") {
            conda_python <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python")
            if (file.exists(conda_python)) {
                python_candidates <- c(conda_python, python_candidates)
            }
        }
        
        # Test each candidate
        python_cmd <- NULL
        for (candidate in python_candidates) {
            test_result <- tryCatch({
                system2(candidate, args = "--version", stdout = TRUE, stderr = TRUE)
            }, error = function(e) NULL)
            
            if (!is.null(test_result) && length(test_result) > 0) {
                python_cmd <- candidate
                validation_results$python_setup$python_cmd <- candidate
                validation_results$python_setup$python_version <- test_result[1]
                break
            }
        }
        
        if (is.null(python_cmd)) {
            validation_results$success <- FALSE
            validation_results$error <- "No working Python installation found"
            validation_results$validation_notes <- append(validation_results$validation_notes,
                                                          "Tried: python3, python, /usr/bin/python3, /usr/local/bin/python3")
            
            # Write validation metadata even if failed
            validation_file <- file.path(species_dir, "validation_results.json")
            write_json(validation_results, validation_file, pretty = TRUE, auto_unbox = TRUE)
            return(validation_results)
        }
    }
    
    # Test if COBRApy is available
    cat("Testing COBRApy availability...\n")
    
    # Create a temporary script to test COBRApy
    cobra_test_script <- tempfile(fileext = ".py")
    writeLines(c(
        "try:",
        "    import cobra", 
        "    print('COBRApy version:', cobra.__version__)",
        "except ImportError:",
        "    print('COBRApy not available')",
        "    exit(1)"
    ), cobra_test_script)
    
    cobra_test <- tryCatch({
        system2(python_cmd, args = cobra_test_script, stdout = TRUE, stderr = TRUE)
    }, error = function(e) NULL)
    
    unlink(cobra_test_script)  # Clean up
    
    if (is.null(cobra_test) || any(attr(cobra_test, "status") %in% c(1, 127))) {
        validation_results$success <- FALSE
        validation_results$error <- "COBRApy not available"
        validation_results$validation_notes <- append(validation_results$validation_notes,
                                                      "COBRApy import failed - please install with: pip install cobra")
        validation_results$python_setup$cobra_available <- FALSE
        
        if (!is.null(cobra_test)) {
            validation_results$cobra_test_output <- paste(cobra_test, collapse = "\n")
        }
        
        # Write validation metadata
        validation_file <- file.path(species_dir, "validation_results.json")
        write_json(validation_results, validation_file, pretty = TRUE, auto_unbox = TRUE)
        return(validation_results)
    }
    
    validation_results$python_setup$cobra_available <- TRUE
    validation_results$python_setup$cobra_version <- paste(cobra_test, collapse = "\n")
    
    # Create Python validation script
    temp_script <- tempfile(fileext = ".py")
    
    python_code <- sprintf('
import sys
import json
import warnings
import os

# Aggressive warning suppression
warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

# Redirect stderr to devnull to suppress COBRApy warnings
import contextlib
from io import StringIO

@contextlib.contextmanager
def suppress_stderr():
    with open(os.devnull, "w") as devnull:
        old_stderr = sys.stderr
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stderr = old_stderr

try:
    with suppress_stderr():
        import cobra
    cobra_available = True
except ImportError:
    cobra_available = False
    print(json.dumps({"error": "COBRApy not available"}))
    sys.exit(1)

def test_model_growth(filepath):
    try:
        with suppress_stderr():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = cobra.io.read_sbml_model(filepath)
                growth_rate = model.slim_optimize()
        return {
            "readable": True, 
            "growth_rate": float(growth_rate),
            "metabolites": len(model.metabolites),
            "reactions": len(model.reactions)
        }
    except Exception as e:
        return {
            "readable": False, 
            "error": str(e), 
            "growth_rate": None,
            "metabolites": None,
            "reactions": None
        }

# Test input file
with suppress_stderr():
    input_result = test_model_growth("%s")

# Test processed file  
with suppress_stderr():
    processed_result = test_model_growth("%s")

# Output results
results = {
    "input": input_result,
    "processed": processed_result,
    "cobra_version": cobra.__version__
}

print(json.dumps(results))
', input_file, processed_file)
    
    writeLines(python_code, temp_script)
    
    tryCatch({
        cat("Running Python validation script...\n")
        
        # Run Python validation with timeout and stderr suppression
        result <- system2(python_cmd, args = temp_script, stdout = TRUE, stderr = FALSE, timeout = 60)
        
        if (length(result) > 0 && !any(attr(result, "status") %in% c(1, 127))) {
            # Parse results
            result_text <- paste(result, collapse = "\n")
            cat("Python output length:", nchar(result_text), "characters\n")
            
            # Look for JSON at the end of output (skip warnings/errors at the beginning)
            lines <- strsplit(result_text, "\n")[[1]]
            json_line <- NULL
            
            # Look for the last line that looks like JSON
            for (i in length(lines):1) {
                if (str_detect(lines[i], "^\\s*\\{.*\\}\\s*$")) {
                    json_line <- lines[i]
                    break
                }
            }
            
            if (!is.null(json_line) && nchar(str_trim(json_line)) > 0) {
                tryCatch({
                    json_result <- fromJSON(str_trim(json_line))
                    
                    validation_results$input_readable_cobra <- json_result$input$readable %||% FALSE
                    validation_results$processed_readable_cobra <- json_result$processed$readable %||% FALSE
                    
                    if (!is.null(json_result$input$readable) && json_result$input$readable) {
                        validation_results$growth_rate_input <- json_result$input$growth_rate
                        validation_results$input_model_size <- list(
                            metabolites = json_result$input$metabolites,
                            reactions = json_result$input$reactions
                        )
                    } else {
                        validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                      paste("Input file not readable by COBRA:", json_result$input$error %||% "unknown error"))
                    }
                    
                    if (!is.null(json_result$processed$readable) && json_result$processed$readable) {
                        validation_results$growth_rate_processed <- json_result$processed$growth_rate
                        validation_results$processed_model_size <- list(
                            metabolites = json_result$processed$metabolites,
                            reactions = json_result$processed$reactions
                        )
                    } else {
                        validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                      paste("Processed file not readable by COBRA:", json_result$processed$error %||% "unknown error"))
                    }
                    
                    # Check equivalence
                    # Check equivalence - only when both files are readable
if (!is.null(validation_results$growth_rate_input) && !is.null(validation_results$growth_rate_processed) &&
    validation_results$input_readable_cobra && validation_results$processed_readable_cobra) {
    
    growth_diff <- abs(validation_results$growth_rate_input - validation_results$growth_rate_processed)
    validation_results$simulation_equivalent <- growth_diff < 1e-6
    validation_results$growth_rate_difference <- growth_diff
    
    if (!validation_results$simulation_equivalent) {
        validation_results$validation_notes <- append(validation_results$validation_notes,
                                                      paste("Growth rates differ by", round(growth_diff, 8)))
    }
} else if (!validation_results$input_readable_cobra && !validation_results$processed_readable_cobra) {
    # Both files unreadable - cannot assess pipeline impact
    validation_results$simulation_equivalent <- NA
    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                  "Cannot assess equivalence - both files unreadable by COBRApy")
} else if (!validation_results$input_readable_cobra) {
    # Input unreadable but processed is readable - potential pipeline improvement
    validation_results$simulation_equivalent <- NA
    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                  "Input file unreadable - pipeline may have improved readability")
} else {
    # Input readable but processed is not - potential pipeline issue
    validation_results$simulation_equivalent <- FALSE
    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                  "Processing broke COBRApy readability - pipeline issue")
}
                    validation_results$success <- TRUE
                    
                }, error = function(e) {
                    validation_results$success <- FALSE
                    validation_results$error <- paste("JSON parsing error:", e$message)
                    validation_results$json_line <- json_line
                    # Save first and last few lines of Python output for debugging
                    all_lines <- strsplit(result_text, "\n")[[1]]
                    validation_results$python_output_sample <- list(
                        first_10_lines = head(all_lines, 10),
                        last_10_lines = tail(all_lines, 10),
                        total_lines = length(all_lines)
                    )
                })
            } else {
                validation_results$success <- FALSE
                validation_results$error <- "No JSON found in Python output"
                # Save first and last few lines for debugging
                all_lines <- strsplit(result_text, "\n")[[1]]
                validation_results$python_output_sample <- list(
                    first_10_lines = head(all_lines, 10),
                    last_10_lines = tail(all_lines, 10),
                    total_lines = length(all_lines)
                )
            }
            
        } else {
            validation_results$success <- FALSE
            validation_results$error <- "Python validation script failed"
            validation_results$python_exit_status <- attr(result, "status")
            validation_results$python_output <- paste(result, collapse = "\n")
            validation_results$validation_notes <- append(validation_results$validation_notes,
                                                          "Check that COBRApy is installed and files are valid SBML")
        }
        
    }, error = function(e) {
        validation_results$success <- FALSE
        validation_results$error <- paste("R error running Python:", e$message)
    }, finally = {
        unlink(temp_script)
    })
    
    # Write validation metadata
    validation_file <- file.path(species_dir, "validation_results.json")
    write_json(validation_results, validation_file, pretty = TRUE, auto_unbox = TRUE)
    
    # Print summary
    if (validation_results$success) {
        cat("Validation completed successfully!\n")
        cat("Input readable:", validation_results$input_readable_cobra, "\n")
        cat("Processed readable:", validation_results$processed_readable_cobra, "\n")
        if (!is.null(validation_results$growth_rate_input)) {
            cat("Input growth rate:", round(validation_results$growth_rate_input, 6), "\n")
        }
        if (!is.null(validation_results$growth_rate_processed)) {
            cat("Processed growth rate:", round(validation_results$growth_rate_processed, 6), "\n")
        }
        if (!is.null(validation_results$simulation_equivalent)) {
            cat("Simulation equivalent:", validation_results$simulation_equivalent, "\n")
        }
    } else {
        cat("Validation failed:", validation_results$error, "\n")
    }
    
    return(validation_results)
}


# Can handle multiple subdirectories (e.g. curated models) or batch from flat directory (e.g. carvefungi)
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
    
    # Calculate statistics with proper logical handling and safe value extraction
    successful <- sum(sapply(all_results, function(x) {
        success_val <- x$success
        if (is.null(success_val) || !is.logical(success_val)) return(FALSE)
        return(success_val)
    }))
    
    input_readable <- sum(sapply(all_results, function(x) {
        readable_val <- x$input_readable_cobra
        if (is.null(readable_val) || !is.logical(readable_val)) return(FALSE)
        return(readable_val)
    }))
    
    processed_readable <- sum(sapply(all_results, function(x) {
        readable_val <- x$processed_readable_cobra
        if (is.null(readable_val) || !is.logical(readable_val)) return(FALSE)
        return(readable_val)
    }))
    
    # Safe equivalent calculation - only count TRUE values
    equivalent <- sum(sapply(all_results, function(x) {
        equiv_val <- x$simulation_equivalent
        if (is.null(equiv_val) || !is.logical(equiv_val) || is.na(equiv_val)) return(FALSE)
        return(equiv_val)
    }))

    pipeline_improved <- sum(sapply(all_results, function(x) {
        input_readable <- x$input_readable_cobra
        processed_readable <- x$processed_readable_cobra
        if (is.null(input_readable) || is.null(processed_readable) || 
            !is.logical(input_readable) || !is.logical(processed_readable)) return(FALSE)
        return(!input_readable && processed_readable)
    }))

    pipeline_broke <- sum(sapply(all_results, function(x) {
        input_readable <- x$input_readable_cobra
        processed_readable <- x$processed_readable_cobra
        if (is.null(input_readable) || is.null(processed_readable) || 
            !is.logical(input_readable) || !is.logical(processed_readable)) return(FALSE)
        return(input_readable && !processed_readable)
    }))

    both_unreadable <- sum(sapply(all_results, function(x) {
        input_readable <- x$input_readable_cobra
        processed_readable <- x$processed_readable_cobra
        if (is.null(input_readable) || is.null(processed_readable) || 
            !is.logical(input_readable) || !is.logical(processed_readable)) return(FALSE)
        return(!input_readable && !processed_readable)
    }))

    report_lines <- c(report_lines,
                      paste("- Validation successful:", successful, "/", length(all_results)),
                      paste("- Input files COBRA-readable:", input_readable, "/", successful),
                      paste("- Processed files COBRA-readable:", processed_readable, "/", successful),
                      paste("- Simulation equivalent (when comparable):", equivalent, "/", processed_readable),
                      paste("- Pipeline improved readability:", pipeline_improved),
                      paste("- Pipeline broke readability:", pipeline_broke),
                      paste("- Both files unreadable:", both_unreadable),
                      ""
    )
    
    # Individual species results
    report_lines <- c(report_lines, "## Individual Results")
    
    for (species_name in names(all_results)) {
        result <- all_results[[species_name]]
        
        # Safe logical extraction with validation
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
            if (!is.null(equiv_val) && is.logical(equiv_val) && !is.na(equiv_val) && equiv_val) {
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
    
    # Non-equivalent species - safe processing of simulation_equivalent
    non_equiv <- names(all_results)[sapply(all_results, function(x) {
        success_val <- x$success
        equiv_val <- x$simulation_equivalent
        
        # Only consider successful validations
        if (is.null(success_val) || !is.logical(success_val) || !success_val) return(FALSE)
        
        # Check if equivalent - return TRUE if not equivalent or can't determine
        if (is.null(equiv_val) || !is.logical(equiv_val) || is.na(equiv_val)) return(TRUE)
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