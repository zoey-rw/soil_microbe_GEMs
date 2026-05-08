#!/usr/bin/env Rscript
#
# tests/_shim_equivalence_dump.R
#
# Dump slot vectors from a single SBML file via either the patched sybilSBML
# back-end or the cobra-via-reticulate shim. Output is JSON on stdout — the
# Python harness in tests/test_shim_equivalence.py invokes this helper twice
# per species (once per backend) and diffs the JSON.
#
# Usage:
#   Rscript tests/_shim_equivalence_dump.R <backend> <input_sbml> <output_json>
#
# Where <backend> is one of:
#   sybilsbml   — load library(sybilSBML) and call readSBMLmod() directly.
#   shim        — source pipeline/sbml_io_cobra.R and call its readSBMLmod().
#
# The helper is deliberately a separate Rscript (not reticulate-from-pytest)
# so that pytest never has to load the Python interpreter that R itself
# embeds via reticulate. Backend processes are isolated; their crashes don't
# leak into the test session.
#
# Exits non-zero on any read failure. Writes a JSON object to <output_json>
# with keys:
#   met_id, met_name, react_id, react_rev, gpr,
#   mod_compart, met_comp,
#   met_compart      — derived: mod_compart[met_comp] (string per metabolite)
#   met_annotation   — character vector aligned with met_id
#   backend, source_file, n_mets, n_reacts

suppressPackageStartupMessages({
    library(jsonlite)
    library(methods)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("usage: Rscript _shim_equivalence_dump.R <backend> <input_sbml> <output_json>")
}
backend     <- args[[1]]
input_sbml  <- normalizePath(args[[2]], mustWork = TRUE)
output_json <- args[[3]]

# Resolve repo root from the script's own location so we can find the shim.
script_path <- (function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
    }
    NA_character_
})()
repo_root <- if (!is.na(script_path)) {
    normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
    getwd()
}

read_via_sybilsbml <- function(path) {
    suppressPackageStartupMessages(library(sybilSBML))
    suppressWarnings(readSBMLmod(path))
}

read_via_shim <- function(path) {
    Sys.setenv(SOIL_MICROBE_GEMS_USE_COBRA_SHIM = "1")
    shim_file <- file.path(repo_root, "pipeline", "sbml_io_cobra.R")
    if (!file.exists(shim_file)) {
        stop("cobra shim not found at: ", shim_file)
    }
    suppressWarnings(source(shim_file, local = FALSE))
    suppressWarnings(readSBMLmod(path))
}

backend <- tolower(backend)
if (backend == "sybilsbml") {
    model <- read_via_sybilsbml(input_sbml)
} else if (backend == "shim" || backend == "cobra" || backend == "cobra_shim") {
    model <- read_via_shim(input_sbml)
} else {
    stop("unknown backend: ", backend, " (expected 'sybilsbml' or 'shim')")
}

# --- Pull annotation vector. Layout differs between back-ends: sybilSBML
# stores it as a column in a data.frame at @met_attr; the shim stores it as
# a character vector under met_attr$annotation. Handle both.
extract_annotation <- function(model) {
    ma <- model@met_attr
    if (is.null(ma)) return(rep("", length(model@met_id)))
    if (is.data.frame(ma) && "annotation" %in% colnames(ma)) {
        return(as.character(ma$annotation))
    }
    if (is.list(ma) && "annotation" %in% names(ma)) {
        return(as.character(ma$annotation))
    }
    rep("", length(model@met_id))
}

annotation <- extract_annotation(model)

mod_compart <- as.character(model@mod_compart)
met_comp    <- as.integer(model@met_comp)
# Build a string compartment for each metabolite. If met_comp index is 0 or
# out-of-range, emit empty string so downstream comparison still aligns.
safe_idx <- function(i) {
    if (is.na(i) || i < 1L || i > length(mod_compart)) "" else mod_compart[i]
}
met_compart <- vapply(met_comp, safe_idx, character(1))

result <- list(
    backend        = backend,
    source_file    = input_sbml,
    n_mets         = length(model@met_id),
    n_reacts       = length(model@react_id),
    met_id         = as.character(model@met_id),
    met_name       = as.character(model@met_name),
    met_comp       = met_comp,
    met_compart    = met_compart,
    mod_compart    = mod_compart,
    met_annotation = annotation,
    react_id       = as.character(model@react_id),
    react_rev      = as.logical(model@react_rev),
    gpr            = as.character(model@gpr)
)

# auto_unbox FALSE so single-element vectors stay arrays (prevents [["x"]] vs
# "x" surprises on the Python side).
out_dir <- dirname(output_json)
if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
writeLines(jsonlite::toJSON(result, auto_unbox = FALSE, null = "null",
                            na = "string", pretty = FALSE), output_json)
cat(sprintf("[dump] %s: %d mets, %d reacts -> %s\n",
            backend, result$n_mets, result$n_reacts, output_json))
