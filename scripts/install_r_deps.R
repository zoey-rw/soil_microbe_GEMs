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

# Archived: sybil installs cleanly from the GitHub mirror.
if (!"sybil" %in% rownames(installed.packages())) {
    cat("Installing archived sybil from cran/sybil\n")
    tryCatch(
        remotes::install_github("cran/sybil", upgrade = "never"),
        error = function(e) {
            message("Could not install sybil: ", e$message)
        }
    )
}

# Archived: sybilSBML 3.1.2 needs source patches to build against libSBML 5.19+.
# Use the vendored installer rather than direct GitHub install.
if (!"sybilSBML" %in% rownames(installed.packages())) {
    repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                               mustWork = FALSE)
    installer <- file.path(repo_root, "pipeline", "vendor", "install_sybilSBML.sh")
    if (file.exists(installer)) {
        cat("Installing sybilSBML via the vendored patched installer:\n  ", installer, "\n")
        rc <- system2("bash", installer)
        if (rc != 0) {
            message("Vendored sybilSBML install failed (exit ", rc, "). ",
                    "Set SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1 to use the cobra shim instead.")
        }
    } else {
        message("sybilSBML not installed and vendored installer not found at ", installer)
    }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("\nDone. To verify:\n  R -e 'library(stringr); library(jsonlite)'\n")
