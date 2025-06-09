# Exchange Metabolite Analysis Script
# Analyzes processed SBML files to assess standardization success for exchange metabolites

library(here)
library(stringr)
library(tidyverse)
library(sybilSBML)
library(jsonlite)
library(yaml)

#' Analyze exchange metabolites across all processed SBML files
#' @param species_base_dir Path to species directory (default: "species")
#' @param output_dir Directory to save analysis results (default: "analysis_results")
#' @return List with comprehensive analysis results
analyze_exchange_metabolites <- function(species_base_dir = "species", 
                                         output_dir = "analysis_results") {
    
    cat("=== Exchange Metabolite Standardization Analysis ===\n")
    
    # Create output directory
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Find all processed files
    species_dirs <- list.dirs(species_base_dir, recursive = FALSE, full.names = TRUE)
    processed_files <- c()
    
    for (species_dir in species_dirs) {
        processed_candidates <- list.files(species_dir, pattern = "_processed\\.xml$", full.names = TRUE)
        if (length(processed_candidates) > 0) {
            processed_files <- c(processed_files, processed_candidates[1])  # Take first if multiple
        }
    }
    
    cat("Found", length(processed_files), "processed SBML files\n")
    
    if (length(processed_files) == 0) {
        stop("No processed SBML files found in ", species_base_dir)
    }
    
    # Initialize analysis data structures
    analysis_results <- list(
        total_models = length(processed_files),
        model_summaries = list(),
        exchange_metabolite_database = data.frame(),
        standardization_summary = list(),
        cross_model_analysis = list()
    )
    
    all_exchange_metabolites <- data.frame(
        model_id = character(),
        metabolite_id = character(),
        base_metabolite = character(),
        is_metanetx = logical(),
        appears_standardized = logical(),
        stringsAsFactors = FALSE
    )
    
    # Process each model
    cat("\nProcessing models:\n")
    
    for (i in seq_along(processed_files)) {
        file_path <- processed_files[i]
        species_name <- basename(dirname(file_path))
        model_id <- str_extract(basename(file_path), "^[^_]+")
        
        cat("  [", i, "/", length(processed_files), "]", model_id, "\n")
        
        tryCatch({
            # Read SBML model
            sbml_model <- readSBMLmod(file_path, validateSBML = FALSE)
            
            # Extract exchange metabolites (ending in "_e")
            all_metabolites <- sbml_model@met_id
            exchange_metabolites <- all_metabolites[str_detect(all_metabolites, "_e\\]?$")]
            
            if (length(exchange_metabolites) == 0) {
                # Try alternative patterns: [e], (e), _e0, etc.
                exchange_metabolites <- all_metabolites[str_detect(all_metabolites, "\\[e\\]$|\\(e\\)$|_e0$")]
            }
            
            # Analyze each exchange metabolite
            model_exchange_data <- data.frame(
                model_id = model_id,
                metabolite_id = exchange_metabolites,
                stringsAsFactors = FALSE
            )
            
            if (nrow(model_exchange_data) > 0) {
                # Extract base metabolite name (remove compartment suffix)
                model_exchange_data$base_metabolite <- str_replace(model_exchange_data$metabolite_id, 
                                                                   "\\[e\\]$|\\(e\\)$|_e$|_e0$", "")
                
                # Check if appears to be MetanetX format (MNXM followed by numbers)
                model_exchange_data$is_metanetx <- str_detect(model_exchange_data$base_metabolite, "^MNXM\\d+$")
                
                # Check if appears standardized (either MetanetX or common biochemical names)
                common_metabolites <- c("glc__D", "ac", "etoh", "succ", "pyr", "lac__L", "for", "co2", 
                                        "h2o", "pi", "h", "o2", "nh3", "so4", "h2s", "no3", "no2")
                
                model_exchange_data$appears_standardized <- model_exchange_data$is_metanetx | 
                    model_exchange_data$base_metabolite %in% common_metabolites
                
                # Add to global database
                all_exchange_metabolites <- rbind(all_exchange_metabolites, model_exchange_data)
            }
            
            # Model summary
            model_summary <- list(
                model_id = model_id,
                species_name = species_name,
                file_path = file_path,
                total_metabolites = length(all_metabolites),
                exchange_metabolites = length(exchange_metabolites),
                metanetx_exchange = sum(model_exchange_data$is_metanetx, na.rm = TRUE),
                standardized_exchange = sum(model_exchange_data$appears_standardized, na.rm = TRUE),
                conversion_rate = if(length(exchange_metabolites) > 0) {
                    round(sum(model_exchange_data$appears_standardized, na.rm = TRUE) / length(exchange_metabolites) * 100, 1)
                } else { 0 }
            )
            
            analysis_results$model_summaries[[model_id]] <- model_summary
            
        }, error = function(e) {
            cat("    ERROR reading", model_id, ":", e$message, "\n")
            
            analysis_results$model_summaries[[model_id]] <- list(
                model_id = model_id,
                species_name = species_name,
                file_path = file_path,
                error = e$message,
                total_metabolites = NA,
                exchange_metabolites = NA,
                metanetx_exchange = NA,
                standardized_exchange = NA,
                conversion_rate = NA
            )
        })
    }
    
    # Store complete exchange metabolite database
    analysis_results$exchange_metabolite_database <- all_exchange_metabolites
    
    # Cross-model analysis
    cat("\nPerforming cross-model analysis...\n")
    
    if (nrow(all_exchange_metabolites) > 0) {
        
        # Metabolite frequency analysis
        base_metabolite_counts <- all_exchange_metabolites %>%
            group_by(base_metabolite) %>%
            summarise(
                frequency = n(),
                models_present = list(unique(model_id)),
                always_metanetx = all(is_metanetx),
                sometimes_metanetx = any(is_metanetx),
                standardization_rate = mean(appears_standardized, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            arrange(desc(frequency))
        
        # Most common exchange metabolites
        analysis_results$cross_model_analysis$most_common_metabolites <- head(base_metabolite_counts, 20)
        
        # Metabolites present in all models
        universal_metabolites <- base_metabolite_counts %>%
            filter(frequency == analysis_results$total_models)
        
        analysis_results$cross_model_analysis$universal_metabolites <- universal_metabolites
        
        # Metabolites unique to single models
        unique_metabolites <- base_metabolite_counts %>%
            filter(frequency == 1)
        
        analysis_results$cross_model_analysis$unique_metabolites <- unique_metabolites
        
        # Standardization success by metabolite
        standardization_by_metabolite <- base_metabolite_counts %>%
            mutate(
                standardization_category = case_when(
                    standardization_rate == 1.0 ~ "Always standardized",
                    standardization_rate > 0.5 ~ "Usually standardized", 
                    standardization_rate > 0 ~ "Sometimes standardized",
                    TRUE ~ "Never standardized"
                )
            )
        
        analysis_results$cross_model_analysis$standardization_by_metabolite <- standardization_by_metabolite
        
        # Overall statistics
        analysis_results$standardization_summary <- list(
            total_exchange_metabolites = nrow(all_exchange_metabolites),
            unique_base_metabolites = length(unique(all_exchange_metabolites$base_metabolite)),
            metanetx_converted = sum(all_exchange_metabolites$is_metanetx, na.rm = TRUE),
            appears_standardized = sum(all_exchange_metabolites$appears_standardized, na.rm = TRUE),
            overall_conversion_rate = round(sum(all_exchange_metabolites$is_metanetx, na.rm = TRUE) / 
                                                nrow(all_exchange_metabolites) * 100, 1),
            overall_standardization_rate = round(sum(all_exchange_metabolites$appears_standardized, na.rm = TRUE) / 
                                                     nrow(all_exchange_metabolites) * 100, 1),
            models_with_good_conversion = sum(sapply(analysis_results$model_summaries, 
                                                     function(x) !is.na(x$conversion_rate) && x$conversion_rate >= 70)),
            universal_metabolite_count = nrow(universal_metabolites),
            unique_metabolite_count = nrow(unique_metabolites)
        )
        
    } else {
        cat("Warning: No exchange metabolites found in any model\n")
    }
    
    # Generate summary report
    cat("\nGenerating summary report...\n")
    
    # Model-level summary table
    model_summary_df <- do.call(rbind, lapply(analysis_results$model_summaries, function(x) {
        data.frame(
            model_id = x$model_id,
            species_name = x$species_name,
            total_metabolites = x$total_metabolites %||% NA,
            exchange_metabolites = x$exchange_metabolites %||% NA,
            metanetx_exchange = x$metanetx_exchange %||% NA,
            conversion_rate = x$conversion_rate %||% NA,
            has_error = !is.null(x$error),
            stringsAsFactors = FALSE
        )
    }))
    
    analysis_results$model_summary_table <- model_summary_df
    
    # Save results
    cat("Saving analysis results...\n")
    
    # Save main results as JSON
    write_json(analysis_results, file.path(output_dir, "exchange_metabolite_analysis.json"), 
               pretty = TRUE, auto_unbox = TRUE)
    
    # Save detailed metabolite database as CSV
    write_csv(all_exchange_metabolites, file.path(output_dir, "exchange_metabolite_database.csv"))
    
    # Save model summaries as CSV
    write_csv(model_summary_df, file.path(output_dir, "model_conversion_summary.csv"))
    
    # Save cross-model analysis
    if (nrow(all_exchange_metabolites) > 0) {
        write_csv(base_metabolite_counts, file.path(output_dir, "metabolite_frequency_analysis.csv"))
        write_csv(universal_metabolites, file.path(output_dir, "universal_metabolites.csv"))
        write_csv(unique_metabolites, file.path(output_dir, "unique_metabolites.csv"))
    }
    
    # Print summary
    print_analysis_summary(analysis_results)
    
    return(analysis_results)
}

#' Print a comprehensive summary of the analysis results
print_analysis_summary <- function(analysis_results) {
    
    cat("\n" , rep("=", 60), "\n")
    cat("EXCHANGE METABOLITE STANDARDIZATION ANALYSIS SUMMARY\n")
    cat(rep("=", 60), "\n")
    
    # Overall statistics
    if (!is.null(analysis_results$standardization_summary)) {
        ss <- analysis_results$standardization_summary
        
        cat("\nOVERALL STATISTICS:\n")
        cat("  Total models analyzed:", analysis_results$total_models, "\n")
        cat("  Total exchange metabolites:", ss$total_exchange_metabolites, "\n")
        cat("  Unique base metabolites:", ss$unique_base_metabolites, "\n")
        cat("  MetanetX converted:", ss$metanetx_converted, "(", ss$overall_conversion_rate, "%)\n")
        cat("  Appears standardized:", ss$appears_standardized, "(", ss$overall_standardization_rate, "%)\n")
        cat("  Models with >70% conversion:", ss$models_with_good_conversion, "out of", analysis_results$total_models, "\n")
        
        cat("\nCROSS-MODEL PATTERNS:\n")
        cat("  Universal metabolites (in all models):", ss$universal_metabolite_count, "\n")
        cat("  Model-specific metabolites:", ss$unique_metabolite_count, "\n")
    }
    
    # Model performance
    if (!is.null(analysis_results$model_summary_table)) {
        model_df <- analysis_results$model_summary_table
        successful_models <- model_df[!model_df$has_error, ]
        
        if (nrow(successful_models) > 0) {
            cat("\nMODEL PERFORMANCE:\n")
            cat("  Best performing models (>90% conversion):\n")
            
            high_performers <- successful_models[!is.na(successful_models$conversion_rate) & 
                                                     successful_models$conversion_rate >= 90, ]
            if (nrow(high_performers) > 0) {
                for (i in 1:min(5, nrow(high_performers))) {
                    cat("    ", high_performers$model_id[i], ":", high_performers$conversion_rate[i], "%\n")
                }
            } else {
                cat("    None found\n")
            }
            
            cat("  Models needing attention (<50% conversion):\n")
            low_performers <- successful_models[!is.na(successful_models$conversion_rate) & 
                                                    successful_models$conversion_rate < 50, ]
            if (nrow(low_performers) > 0) {
                for (i in 1:min(5, nrow(low_performers))) {
                    cat("    ", low_performers$model_id[i], ":", low_performers$conversion_rate[i], "%\n")
                }
            } else {
                cat("    None found\n")
            }
        }
    }
    
    # Most common metabolites
    if (!is.null(analysis_results$cross_model_analysis$most_common_metabolites)) {
        common_mets <- analysis_results$cross_model_analysis$most_common_metabolites
        
        cat("\nMOST COMMON EXCHANGE METABOLITES:\n")
        for (i in 1:min(10, nrow(common_mets))) {
            cat("  ", common_mets$base_metabolite[i], 
                " (", common_mets$frequency[i], " models, ",
                round(common_mets$standardization_rate[i] * 100, 1), "% standardized)\n")
        }
    }
    
    # Files saved
    cat("\nOUTPUT FILES SAVED:\n")
    cat("  - exchange_metabolite_analysis.json (complete results)\n")
    cat("  - exchange_metabolite_database.csv (all exchange metabolites)\n") 
    cat("  - model_conversion_summary.csv (per-model statistics)\n")
    cat("  - metabolite_frequency_analysis.csv (cross-model patterns)\n")
    cat("  - universal_metabolites.csv (metabolites in all models)\n")
    cat("  - unique_metabolites.csv (model-specific metabolites)\n")
    
    cat("\n", rep("=", 60), "\n")
}

#' Generate a detailed report for specific metabolites of interest
generate_metabolite_focus_report <- function(analysis_results, 
                                             metabolites_of_interest = c("glc__D", "ac", "etoh", "co2", "o2", "h2o"),
                                             output_dir = "analysis_results") {
    
    if (is.null(analysis_results$exchange_metabolite_database) || 
        nrow(analysis_results$exchange_metabolite_database) == 0) {
        cat("No exchange metabolite data available for focus report\n")
        return(NULL)
    }
    
    focus_data <- analysis_results$exchange_metabolite_database %>%
        filter(base_metabolite %in% metabolites_of_interest) %>%
        group_by(base_metabolite) %>%
        summarise(
            models_present = list(unique(model_id)),
            total_occurrences = n(),
            metanetx_conversions = sum(is_metanetx),
            standardized_instances = sum(appears_standardized),
            conversion_rate = round(sum(is_metanetx) / n() * 100, 1),
            standardization_rate = round(sum(appears_standardized) / n() * 100, 1),
            example_ids = list(head(unique(metabolite_id), 5)),
            .groups = "drop"
        )
    
    write_csv(focus_data, file.path(output_dir, "focus_metabolites_report.csv"))
    
    cat("\nFOCUS METABOLITES REPORT:\n")
    for (i in 1:nrow(focus_data)) {
        met <- focus_data[i, ]
        cat("  ", met$base_metabolite, ":\n")
        cat("    Present in", met$total_occurrences, "model instances\n")
        cat("    MetanetX conversion:", met$metanetx_conversions, "/", met$total_occurrences, 
            "(", met$conversion_rate, "%)\n")
        cat("    Example IDs:", paste(unlist(met$example_ids), collapse = ", "), "\n\n")
    }
    
    return(focus_data)
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# Main execution
if (!interactive()) {
    # Run analysis when script is executed directly
    cat("Starting exchange metabolite analysis...\n")
    results <- analyze_exchange_metabolites()
    
    # Generate focus report for common metabolites
    focus_report <- generate_metabolite_focus_report(results)
    
    cat("\nAnalysis complete! Check the 'analysis_results' directory for detailed outputs.\n")
}



results <- analyze_exchange_metabolites("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/species")
focus_report <- generate_metabolite_focus_report(results)

