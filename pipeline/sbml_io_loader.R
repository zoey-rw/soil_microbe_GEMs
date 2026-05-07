#' Conditional SBML I/O loader.
#'
#' Sourced by process_sbml_species.R, processing_utils.R, and
#' check_exchange_metabolites.R in place of the old `library(sybilSBML)` line.
#'
#' Behaviour:
#'   - SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1 -> source pipeline/sbml_io_cobra.R
#'     (cobra-via-reticulate shim; transitional, see issue #2 / #5).
#'   - Otherwise -> attempt `library(sybilSBML)`. If that fails (e.g. on
#'     a system where sybilSBML can't compile), print a one-line error
#'     pointing at the env-var fallback and re-throw.
#'
#' This indirection means the swap window is a single env var — no R code
#' changes are required to flip between back-ends. Once the equivalence
#' harness validates the shim across the full curated cohort (issue #5),
#' the default flips and the sybilSBML branch can be removed.

local({
    use_shim <- Sys.getenv("SOIL_MICROBE_GEMS_USE_COBRA_SHIM", unset = "")
    if (use_shim == "1" || tolower(use_shim) == "true") {
        message("[sbml_io_loader] using cobra-via-reticulate shim ",
                "(SOIL_MICROBE_GEMS_USE_COBRA_SHIM is set)")
        suppressWarnings(source(here::here("pipeline", "sbml_io_cobra.R")))
    } else {
        tryCatch(
            suppressPackageStartupMessages(library(sybilSBML)),
            error = function(e) {
                stop("sybilSBML failed to load: ", conditionMessage(e),
                     "\nFallback: set SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1 to ",
                     "use the cobra-via-reticulate shim instead. See ",
                     "pipeline/sbml_io_cobra.R and issue #2.",
                     call. = FALSE)
            }
        )
    }
})
