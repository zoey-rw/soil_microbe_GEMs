# Installs R packages required by the pipeline and helper scripts.
#
# Usage:
#   Rscript scripts/install_r_deps.R
#
# Notes:
#  - sybilSBML and sybil were archived from CRAN in 2020. Direct CRAN install
#    will not work; the script will fall back to the `cran` GitHub mirror
#    via remotes::install_github(). The pipeline is in the process of moving
#    off sybilSBML (see issue #2 and the planned cobra-via-reticulate shim);
#    this install path remains for reference and historical reproducibility.
#  - System dependencies on Debian/Ubuntu: libxml2-dev libsbml5-dev libglpk-dev
#    (the Dockerfile in comets_shinyapp_example/ shows the equivalent set
#    for the Shiny side).

cran_packages <- c(
    # Core pipeline
    "here", "stringr", "dplyr", "tidyr", "tidyverse", "readr",
    "data.table", "jsonlite", "yaml", "xml2", "minval",
    # Shiny app
    "shiny", "shinydashboard", "shinyhelper", "DT", "ggplot2",
    "broom", "dtplyr",
    # Stage 1 cobra-via-reticulate shim (planned; safe to install now)
    "reticulate",
    # Optional (used by some preprocessing scripts)
    "segmented", "furrr",
    # Helper for archived-package install
    "remotes"
)

installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)

if (length(missing_cran) > 0) {
    cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Archived packages: sybil + sybilSBML (transitional; see issue #2).
archived <- c(sybil = "cran/sybil", sybilSBML = "cran/sybilSBML")
for (pkg in names(archived)) {
    if (!pkg %in% rownames(installed.packages())) {
        cat("Attempting archived install of", pkg, "from", archived[[pkg]], "\n")
        tryCatch(
            remotes::install_github(archived[[pkg]], upgrade = "never"),
            error = function(e) {
                message("Could not install ", pkg, ": ", e$message,
                        ". The cobra-via-reticulate shim will replace this dependency.")
            }
        )
    }
}

cat("\nDone. To verify:\n  R -e 'library(stringr); library(jsonlite)'\n")
