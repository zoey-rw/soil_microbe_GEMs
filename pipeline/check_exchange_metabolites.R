# Exchange Metabolite Analysis Script
# Analyzes processed SBML files to assess standardization success for exchange metabolites

library(here)
library(stringr)
library(dplyr)
library(readr)
library(sybilSBML)
library(jsonlite)
library(yaml)

#' Create metabolite description mapping
#' @return Named vector mapping metabolite IDs to human-readable descriptions
create_metabolite_descriptions <- function() {
    
    # MetanetX ID mappings to common names
    metanetx_descriptions <- c(
        "MNXM2" = "Water (H2O)",
        "MNXM4" = "Hydrogen ion (H+)", 
        "MNXM5" = "Ammonia (NH3)",
        "MNXM9" = "Phosphate (Pi)",
        "MNXM11" = "Sulfate (SO4)",
        "MNXM13" = "Carbon dioxide (CO2)",
        "MNXM15" = "Oxygen (O2)",
        "MNXM32" = "Inorganic phosphate (Pi)",
        "MNXM57" = "Acetate",
        "MNXM75" = "Sulfate",
        "MNXM114" = "Ethanol",
        "MNXM210" = "Nitrate (NO3)",
        "MNXM364" = "Nitrite (NO2)",
        "MNXM89557" = "D-Glucose",
        "MNXM89558" = "L-Lactate",
        "MNXM730" = "Hydrogen sulfide (H2S)",
        "MNXM1105" = "Formate",
        "MNXM162260" = "Pyruvate",
        "MNXM162989" = "Succinate",
        "MNXM1369759" = "Fatty acid (long chain)"
    )
    
    # Common biochemical abbreviations
    common_descriptions <- c(
        "glc__D" = "D-Glucose",
        "ac" = "Acetate", 
        "etoh" = "Ethanol",
        "succ" = "Succinate",
        "pyr" = "Pyruvate",
        "lac__L" = "L-Lactate",
        "for" = "Formate",
        "co2" = "Carbon dioxide (CO2)",
        "h2o" = "Water (H2O)",
        "pi" = "Inorganic phosphate (Pi)",
        "h" = "Hydrogen ion (H+)",
        "o2" = "Oxygen (O2)",
        "nh3" = "Ammonia (NH3)",
        "so4" = "Sulfate (SO4)",
        "h2s" = "Hydrogen sulfide (H2S)",
        "no3" = "Nitrate (NO3)",
        "no2" = "Nitrite (NO2)",
        "glc" = "Glucose",
        "fru" = "Fructose",
        "xyl" = "Xylose",
        "ara" = "Arabinose",
        "gal" = "Galactose",
        "man" = "Mannose",
        "mal__L" = "L-Malate",
        "cit" = "Citrate",
        "alpha_KG" = "Alpha-ketoglutarate",
        "ala__L" = "L-Alanine",
        "gly" = "Glycine",
        "ser__L" = "L-Serine",
        "thr__L" = "L-Threonine",
        "val__L" = "L-Valine",
        "leu__L" = "L-Leucine",
        "ile__L" = "L-Isoleucine",
        "phe__L" = "L-Phenylalanine",
        "trp__L" = "L-Tryptophan",
        "tyr__L" = "L-Tyrosine",
        "pro__L" = "L-Proline",
        "his__L" = "L-Histidine",
        "lys__L" = "L-Lysine",
        "arg__L" = "L-Arginine",
        "asp__L" = "L-Aspartate",
        "asn__L" = "L-Asparagine",
        "glu__L" = "L-Glutamate",
        "gln__L" = "L-Glutamine",
        "met__L" = "L-Methionine",
        "cys__L" = "L-Cysteine"
    )
    
    # Combine all descriptions
    all_descriptions <- c(metanetx_descriptions, common_descriptions)
    
    return(all_descriptions)
}

#' Get human-readable description for a metabolite
#' @param metabolite_id The metabolite identifier
#' @param descriptions Named vector of descriptions
#' @return Human-readable description or the original ID if no mapping exists
get_metabolite_description <- function(metabolite_id, descriptions) {
    desc <- descriptions[metabolite_id]
    if (is.na(desc)) {
        # Try without compartment suffixes
        base_id <- str_replace(metabolite_id, "\\[e\\]$|\\(e\\)$|_e$|_e0$", "")
        desc <- descriptions[base_id]
        if (is.na(desc)) {
            return(metabolite_id)  # Return original if no description found
        }
    }
    return(desc)
}

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
        
        # Metabolite frequency analysis - USE BASE METABOLITE (converted names)
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
        cat("  Models with >70% conversion:", ss$models_with_good_conversion, "/", analysis_results$total_models, "\n")
        cat("  Universal metabolites:", ss$universal_metabolite_count, "\n")
        cat("  Model-specific metabolites:", ss$unique_metabolite_count, "\n")
        
        # Most common metabolites
        if (!is.null(analysis_results$cross_model_analysis$most_common_metabolites)) {
            common_mets <- analysis_results$cross_model_analysis$most_common_metabolites
            cat("\nMOST COMMON METABOLITES (top 10):\n")
            for (i in 1:min(10, nrow(common_mets))) {
                met <- common_mets[i, ]
                cat("  ", i, ".", met$base_metabolite, "(", met$frequency, "models,", 
                    round(met$standardization_rate * 100, 1), "% standardized)\n")
            }
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

#' Generate detailed report with human-readable metabolite descriptions
#' @param analysis_results Results from analyze_exchange_metabolites
#' @param metabolites_of_interest Vector of MetanetX IDs or common names to focus on
#' @param output_dir Output directory
#' @return Focus report data with descriptions
generate_metabolite_focus_report <- function(analysis_results, 
                                             metabolites_of_interest = NULL, 
                                             output_dir = "analysis_results") {
    
    if (is.null(analysis_results$exchange_metabolite_database) || 
        nrow(analysis_results$exchange_metabolite_database) == 0) {
        cat("No exchange metabolite data available for focus report\n")
        return(NULL)
    }
    
    # Create metabolite descriptions
    descriptions <- create_metabolite_descriptions()
    
    # If no specific metabolites provided, use the most common ones
    if (is.null(metabolites_of_interest)) {
        # Get the top 20 most common metabolites from the analysis
        common_mets <- analysis_results$cross_model_analysis$most_common_metabolites
        if (!is.null(common_mets) && nrow(common_mets) > 0) {
            metabolites_of_interest <- head(common_mets$base_metabolite, 20)
            cat("Using top 20 most common metabolites for focus report\n")
        } else {
            cat("No common metabolites found for focus report\n")
            return(NULL)
        }
    }
    
    # Key metabolites of biological interest
    key_metabolites <- c(
        "MNXM4", "MNXM2", "MNXM15", "MNXM13", "MNXM89557", "MNXM57", "MNXM114", 
        "MNXM89558", "MNXM32", "MNXM5", "MNXM75", "h", "h2o", "o2", "co2", 
        "glc__D", "ac", "etoh", "pi", "nh3", "so4"
    )
    
    # Combine user-specified and key metabolites
    all_focus_metabolites <- unique(c(metabolites_of_interest, key_metabolites))
    
    focus_data <- analysis_results$exchange_metabolite_database %>%
        filter(base_metabolite %in% all_focus_metabolites) %>%
        group_by(base_metabolite) %>%
        summarise(
            models_present = list(unique(model_id)),
            total_occurrences = n(),
            metanetx_conversions = sum(is_metanetx),
            standardized_instances = sum(appears_standardized),
            conversion_rate = round(sum(is_metanetx) / n() * 100, 1),
            standardization_rate = round(sum(appears_standardized) / n() * 100, 1),
            example_ids = list(head(unique(metabolite_id), 5)),
            models_list = paste(unique(model_id), collapse = ", "),
            num_models = length(unique(model_id)),
            .groups = "drop"
        ) %>%
        arrange(desc(total_occurrences))
    
    # Add human-readable descriptions
    focus_data$description <- sapply(focus_data$base_metabolite, 
                                     function(x) get_metabolite_description(x, descriptions))
    
    # Reorder columns to put description second
    focus_data <- focus_data %>%
        select(base_metabolite, description, everything())
    
    write_csv(focus_data, file.path(output_dir, "focus_metabolites_report.csv"))
    
    cat("\n", rep("=", 60), "\n")
    cat("FOCUS METABOLITES REPORT WITH DESCRIPTIONS\n")
    cat(rep("=", 60), "\n")
    cat("Analyzed", nrow(focus_data), "metabolites of interest\n\n")
    
    # Report high-priority metabolites first
    priority_metabolites <- c("MNXM15", "MNXM13", "MNXM2", "MNXM4", "MNXM89557", 
                              "o2", "co2", "h2o", "h", "glc__D")
    
    cat("KEY BIOLOGICAL METABOLITES:\n")
    priority_data <- focus_data %>% filter(base_metabolite %in% priority_metabolites)
    
    if (nrow(priority_data) > 0) {
        for (i in 1:nrow(priority_data)) {
            met <- priority_data[i, ]
            cat("• ", met$description, " (", met$base_metabolite, "):\n", sep="")
            cat("  Present in ", met$num_models, " models (", met$total_occurrences, " instances)\n", sep="")
            cat("  Standardization: ", met$standardization_rate, "% (", met$standardized_instances, "/", met$total_occurrences, ")\n", sep="")
            
            # Identify models missing this metabolite
            all_models <- unique(analysis_results$exchange_metabolite_database$model_id)
            missing_models <- setdiff(all_models, unlist(met$models_present))
            if (length(missing_models) > 0 && length(missing_models) <= 10) {
                cat("  Missing from: ", paste(head(missing_models, 5), collapse=", "), 
                    if(length(missing_models) > 5) paste(" (+", length(missing_models)-5, "more)") else "", "\n", sep="")
            } else if (length(missing_models) > 10) {
                cat("  Missing from ", length(missing_models), " models\n", sep="")
            }
            cat("\n")
        }
    }
    
    cat("OTHER COMMON METABOLITES:\n")
    other_data <- focus_data %>% filter(!base_metabolite %in% priority_metabolites) %>% head(15)
    
    for (i in 1:min(15, nrow(other_data))) {
        met <- other_data[i, ]
        cat("• ", met$description, " (", met$base_metabolite, "):\n", sep="")
        cat("  ", met$num_models, " models, ", met$standardization_rate, "% standardized\n", sep="")
    }
    
    cat("\n", rep("=", 60), "\n")
    
    return(focus_data)
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# Main execution
if (!interactive()) {
    # Run analysis when script is executed directly
    cat("Starting exchange metabolite analysis...\n")
    results <- analyze_exchange_metabolites(species_base_dir = "../species", output_dir = "../analysis_results")
    
    # Generate focus report for common metabolites
    focus_report <- generate_metabolite_focus_report(results, output_dir = "../analysis_results")
    
    cat("\nAnalysis complete! Check the 'analysis_results' directory for detailed outputs.\n")
}
