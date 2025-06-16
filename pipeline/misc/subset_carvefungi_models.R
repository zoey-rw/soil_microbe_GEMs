library(dplyr)
library(stringr)

#' Process CarveFungi directory to extract final ensemble models
#' @param input_dir Directory containing all CarveFungi SBML files
#' @param output_dir Directory to create for final ensemble models (default: "ensemble")
#' @param dry_run If TRUE, show what would be done without actually moving/deleting files
#' @return Summary of processing results
#' 
process_carvefungi_directory <- function(input_dir, output_dir = "ensemble", dry_run = FALSE) {
  
  cat("=== CarveFungi Directory Processing ===\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Dry run mode:", dry_run, "\n\n")
  
  # Get all SBML files
  all_files <- list.files(input_dir, pattern = "\\.sbml$", full.names = FALSE)
  
  if (length(all_files) == 0) {
    stop("No SBML files found in ", input_dir)
  }
  
  cat("Found", length(all_files), "SBML files\n")
  
  # Parse all filenames to identify final ensemble models
  file_data <- data.frame(
    filename = all_files,
    full_path = file.path(input_dir, all_files),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # Identify final ensemble models (double-dash pattern)
      is_final_ensemble = str_detect(filename, "--[0-9.]+\\-[0-9.]+\\.sbml$"),
      # Extract genome identifier (everything before first dash with scores)
      genome_id = str_replace(filename, "(-[OM]-[0-9.]+\\.sbml$|--[0-9.]+\\-[0-9.]+\\.sbml$)", ""),
      # Extract dual scores for ranking if multiple finals exist per genome
      dual_scores = str_extract(filename, "--([0-9.]+)\\-([0-9.]+)\\.sbml$"),
      score1 = as.numeric(str_extract(dual_scores, "(?<=--)[0-9.]+")),
      score2 = as.numeric(str_extract(dual_scores, "(?<=-)[0-9.]+(?=\\.sbml)")),
      combined_score = ifelse(is_final_ensemble, score1 + score2, NA)
    )
  
  # Summary of file types
  final_ensemble_count <- sum(file_data$is_final_ensemble)
  intermediate_count <- sum(!file_data$is_final_ensemble)
  unique_genomes <- length(unique(file_data$genome_id))
  
  cat("\nFile type summary:\n")
  cat("- Final ensemble models:", final_ensemble_count, "\n")
  cat("- Intermediate models:", intermediate_count, "\n")
  cat("- Unique genomes:", unique_genomes, "\n")
  
  # Parse individual model scores for fallback
  file_data <- file_data %>%
    mutate(
      # Extract individual model info
      individual_pattern = str_extract(filename, "-[OM]-[0-9.]+\\.sbml$"),
      bounds_type = str_extract(individual_pattern, "[OM]"),
      individual_score = as.numeric(str_extract(individual_pattern, "[0-9.]+")),
      # Prefer O (open) over M (marginal/closed) bounds
      bounds_priority = ifelse(bounds_type == "O", 2, ifelse(bounds_type == "M", 1, 0))
    )
  
  # Select best model per genome (ensemble first, then best individual)
  final_models <- file_data %>%
    group_by(genome_id) %>%
    arrange(
      desc(is_final_ensemble),           # Ensemble models first
      desc(combined_score),              # Best ensemble score
      desc(bounds_priority),             # Then O models over M models  
      desc(individual_score)             # Then highest individual score
    ) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Summary of selection types
  ensemble_selected <- sum(final_models$is_final_ensemble, na.rm = TRUE)
  individual_selected <- sum(!final_models$is_final_ensemble, na.rm = TRUE)
  
  cat("\nModel selection summary:\n")
  cat("- Final ensemble models selected:", ensemble_selected, "\n")
  cat("- Individual models selected (fallback):", individual_selected, "\n")
  
  if (individual_selected > 0) {
    # Show breakdown of individual model types
    individual_breakdown <- final_models %>%
      filter(!is_final_ensemble) %>%
      mutate(bounds_display = ifelse(is.na(bounds_type), "Unknown", bounds_type)) %>%
      count(bounds_display, name = "count")
    
    cat("  Individual model breakdown:\n")
    for (i in 1:nrow(individual_breakdown)) {
      bounds <- individual_breakdown$bounds_display[i]
      count <- individual_breakdown$count[i]
      bounds_desc <- case_when(
        bounds == "O" ~ "Open bounds",
        bounds == "M" ~ "Marginal/closed bounds", 
        TRUE ~ "Unknown pattern"
      )
      cat("    -", bounds, "(", bounds_desc, "):", count, "\n")
    }
  }
  
  # Check for unusual patterns
  unusual_files <- final_models %>%
    filter(!is_final_ensemble & is.na(bounds_type))
  
  if (nrow(unusual_files) > 0) {
    cat("\nNote:", nrow(unusual_files), "files selected with non-standard naming patterns\n")
    if (nrow(unusual_files) <= 5) {
      cat("Examples:\n")
      for (f in unusual_files$filename) {
        cat("  -", f, "\n")
      }
    }
  }
  
  # Prepare file operations
  if (!dry_run) {
    # Create output directory
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("\nCreated output directory:", output_dir, "\n")
    }
  }
  
  # Initialize processed files dataframe with model_type column
  processed_files <- data.frame(
    original_name = character(),
    new_name = character(),
    genome_id = character(),
    combined_score = numeric(),
    model_type = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(final_models)) {
    row <- final_models[i, ]
    
    # Create new filename: genome_id.sbml
    new_filename <- paste0(row$genome_id, ".sbml")
    new_path <- file.path(output_dir, new_filename)
    
    processed_files <- rbind(processed_files, data.frame(
      original_name = row$filename,
      new_name = new_filename,
      genome_id = row$genome_id,
      combined_score = ifelse(is.na(row$combined_score), row$individual_score, row$combined_score),
      model_type = ifelse(row$is_final_ensemble, "Ensemble", paste0("Individual-", row$bounds_type)),
      stringsAsFactors = FALSE
    ))
    
    if (!dry_run) {
      # Copy and rename file
      file.copy(row$full_path, new_path, overwrite = TRUE)
    }
    
    # Progress indicator
    if (i %% 100 == 0 || i == nrow(final_models)) {
      cat("Processed", i, "of", nrow(final_models), "files\n")
    }
  }
  
  # Delete intermediate files (everything not selected)
  intermediate_files <- file_data %>%
    filter(!filename %in% final_models$filename)
  
  files_to_delete <- nrow(intermediate_files)
  
  if (!dry_run && files_to_delete > 0) {
    cat("\nDeleting", files_to_delete, "intermediate files...\n")
    
    deleted_count <- 0
    for (file_path in intermediate_files$full_path) {
      if (file.exists(file_path)) {
        unlink(file_path)
        deleted_count <- deleted_count + 1
      }
    }
    
    cat("Deleted", deleted_count, "files\n")
  }
  
  # Final summary
  cat("\n=== Processing Summary ===\n")
  cat("Total input files:", length(all_files), "\n")
  cat("Final models selected:", nrow(final_models), "\n")
  cat("Intermediate/unselected files to delete:", files_to_delete, "\n")
  cat("Model type breakdown:\n")
  
  type_summary <- table(processed_files$model_type)
  for (i in 1:length(type_summary)) {
    cat("  -", names(type_summary)[i], ":", type_summary[i], "\n")
  }
  
  if (dry_run) {
    cat("\n** DRY RUN MODE - No files were actually moved or deleted **\n")
  } else {
    cat("\nFinal models saved to:", output_dir, "\n")
  }
  
  return(list(
    final_models = processed_files,
    total_input_files = length(all_files),
    intermediate_files_deleted = files_to_delete,
    output_directory = output_dir,
    selection_summary = list(
      ensemble_selected = ensemble_selected,
      individual_selected = individual_selected,
      total_genomes = length(unique(processed_files$genome_id))
    )
  ))
}
#
# # Run in dry mode first to see what would happen
# results <- process_carvefungi_directory("~/Downloads/CarveFungi", 
# "~/soil_microbe_GEMs/noncurated_species", 
#                                        dry_run = TRUE)

# results <- process_carvefungi_directory("~/Downloads/CarveFungi", 
# "~/soil_microbe_GEMs/noncurated_species", 
#                                        dry_run = FALSE)

# Move soil-related species into their own folder to limit what is being pushed to Github

# Permissive list of genera found in soil...from various sources/domain knowledge
soil_genera <- c("Agaricus", "Amanita", "Cenococcum", "Colletotrichum", "Elaphomyces", 
"Gigaspora", "Laetiporus", "Oidiodendron", "Pleurotus", "Suillus", 
"Tremella", "Xylaria", "Trichoderma", "Penicillium", "Fusarium", 
"Mortierella", "Rhizophagus", "Glomus", "Chaetomium", "Alternaria", 
"Cladosporium", "Mucor", "Rhizopus", "Acremonium", "Geosmithia", 
"Pseudallescheria", "Scedosporium", "Phialophora", "Cladophialophora", 
"Aspergillus", "Clavispora", "Bipolaris", "Aureobasidium", "Cryptococcus", 
"Exophiala", "Gymnopus", "Gymnopilus", "Leucosporidium", "Metarhizium", 
"Rhodotorula", "Melanogaster", "Tomentella", "Piloderma", "Cortinarius", 
"Entoloma", "Tricholoma", "Scleroderma", "Gymnomyces")

# Moving 518 SBML files - multiple genomes per species
sbml_list = list.files("~/soil_microbe_GEMs/noncurated_species", pattern = paste0(soil_genera,"",collapse = "|"), full.names = T)

new_names = gsub("noncurated_species", "carvefungi_species", sbml_list)
file.rename(sbml_list, new_names)

