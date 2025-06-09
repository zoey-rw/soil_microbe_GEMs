# validate_model_growth.R
# COBRApy validation functions for SBML models

library(stringr)
library(jsonlite)

# Source utility functions
source("./pipeline/sbml_processing_utils.R")

#' Validate SBML model growth in COBRApy and update validation metadata
#' @param species_dir Path to species directory
#' @param python_cmd Python command to use (e.g., "python3", "python", "/path/to/python")
#' @return List with validation results
validate_model_growth <- function(species_dir, python_cmd = NULL) {
    
    model_id <- str_extract(basename(species_dir), "[^_]+$")
    
    # Find input and processed files
    input_files <- list.files(species_dir, pattern = "\\.xml$", full.names = TRUE)
    processed_files <- input_files[str_detect(basename(input_files), "_processed")]
    
    # Better filtering for input files
    processed_patterns <- c("_processed", "_cobra_validated", "_modified_cobra", "COBRA-sbml3")
    processed_files_all <- input_files[str_detect(basename(input_files), paste(processed_patterns, collapse = "|"))]
    input_candidates <- setdiff(input_files, processed_files_all)
    
    if (length(input_candidates) == 0 || length(processed_files) == 0) {
        return(list(success = FALSE, error = "Missing input or processed files"))
    }
    
    # Use same file selection logic as main processor
    input_file <- input_candidates[1]
    if (length(input_candidates) > 1) {
        model_id_files <- input_candidates[str_detect(basename(input_candidates), model_id)]
        if (length(model_id_files) > 0) {
            file_lengths <- nchar(basename(model_id_files))
            input_file <- model_id_files[which.min(file_lengths)]
        }
    }
    
    processed_file <- processed_files[1]
    
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
warnings.filterwarnings("ignore")

try:
    import cobra
    cobra_available = True
except ImportError:
    cobra_available = False
    print(json.dumps({"error": "COBRApy not available"}))
    sys.exit(1)

def test_model_growth(filepath):
    try:
        # Suppress COBRApy warnings
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
input_result = test_model_growth("%s")

# Test processed file  
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
        
        # Run Python validation with timeout
        result <- system2(python_cmd, args = temp_script, stdout = TRUE, stderr = TRUE, timeout = 60)
        
        if (length(result) > 0 && !any(attr(result, "status") %in% c(1, 127))) {
            # Parse results
            result_text <- paste(result, collapse = "\n")
            
            # Extract JSON from output (in case there are warnings)
            json_match <- regexpr("\\{.*\\}", result_text)
            if (json_match > 0) {
                json_text <- regmatches(result_text, json_match)
                json_result <- fromJSON(json_text)
                
                validation_results$input_readable_cobra <- json_result$input$readable
                validation_results$processed_readable_cobra <- json_result$processed$readable
                
                if (json_result$input$readable) {
                    validation_results$growth_rate_input <- json_result$input$growth_rate
                    validation_results$input_model_size <- list(
                        metabolites = json_result$input$metabolites,
                        reactions = json_result$input$reactions
                    )
                } else {
                    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                  paste("Input file not readable by COBRA:", json_result$input$error))
                }
                
                if (json_result$processed$readable) {
                    validation_results$growth_rate_processed <- json_result$processed$growth_rate
                    validation_results$processed_model_size <- list(
                        metabolites = json_result$processed$metabolites,
                        reactions = json_result$processed$reactions
                    )
                } else {
                    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                  paste("Processed file not readable by COBRA:", json_result$processed$error))
                }
                
                # Check equivalence
                if (!is.null(validation_results$growth_rate_input) && !is.null(validation_results$growth_rate_processed)) {
                    growth_diff <- abs(validation_results$growth_rate_input - validation_results$growth_rate_processed)
                    validation_results$simulation_equivalent <- growth_diff < 1e-6
                    validation_results$growth_rate_difference <- growth_diff
                    
                    if (!validation_results$simulation_equivalent) {
                        validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                      paste("Growth rates differ by", round(growth_diff, 8)))
                    }
                } else {
                    validation_results$simulation_equivalent <- FALSE
                    validation_results$validation_notes <- append(validation_results$validation_notes,
                                                                  "Cannot compare growth rates - one or both files unreadable")
                }
                
                validation_results$success <- TRUE
                
            } else {
                validation_results$success <- FALSE
                validation_results$error <- "Could not parse Python output"
                validation_results$python_output <- result_text
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
