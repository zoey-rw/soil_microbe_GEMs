# Loads validation helpers. Sourced by scripts that want
# run_post_processing_validation() / create_validation_report() /
# check_validation_status(). Requires python3 + cobrapy.
#
# This file no longer runs anything by itself when sourced; use the
# CLI entry points or call the helpers from an interactive session.
#
# Example interactive usage (run from the repo root):
#
#   source("pipeline/batch_validate_all.R")
#
#   # Validate every species under species/ that has a *_processed.xml.
#   validation_results <- run_post_processing_validation(
#       here::here("species"), force_revalidate = TRUE
#   )
#
#   # Aggregate report for a cohort (curated or carvefungi).
#   create_validation_report(
#       here::here("species"),
#       output_file = "curated_model_validation_report.txt"
#   )
#   create_validation_report(
#       here::here("carvefungi_species"),
#       output_file = "carvefungi_validation_report.txt"
#   )

library(here)
source(here("pipeline", "validate_model_growth.r"))     # Validation functions
