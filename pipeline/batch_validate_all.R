# Script to batch validate that models grow appropriately after standardizing formats and converting identifiers to the MetaNetX namespace
# Growth rate should be the same as the input model

source("validate_model_growth.r")     # Validation functions

# requires python3 with cobrapy installation

# Run validation on all processed species
validation_results <- run_post_processing_validation("../species/", force_revalidate = T)

# Check current status
# status <- check_validation_status("../species/")

# Create readable report: if models were improved, consistent, or corrupted throughout the pipeline
report_summary <- create_validation_report("../species")


report_summary <- create_validation_report("./carvefungi_species/", output_file = "carvefungi_validation_report.txt")

report_summary <- create_validation_report("./species/", output_file = "curated_model_validation_report.txt")
