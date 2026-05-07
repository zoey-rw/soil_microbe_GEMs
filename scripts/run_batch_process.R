#!/usr/bin/env Rscript
# CLI entry point for the SBML standardization pipeline.
#
# Usage:
#   Rscript scripts/run_batch_process.R <command> [options]
#
# Commands:
#   curated  Process the one-directory-per-species layout under species/
#   flat     Process a flat directory of SBML files (e.g. carvefungi_species/input/)
#
# Options:
#   --species-dir=PATH   Where to read species directories (curated mode).
#                        Default: $(repo)/species
#   --input-dir=PATH     Where to read SBML files (flat mode).
#                        Default: $(repo)/carvefungi_species/input
#   --output-dir=PATH    Where to write processed SBML (flat mode).
#                        Default: $(repo)/carvefungi_species/processed
#   --filter=PATTERN     Optional regex filter; restricts which species/files
#                        are processed. Curated mode: matches species directory
#                        name. Flat mode: matches filename.
#   --ref-data=PATH      Path to MetanetX reference RDS. Required.
#   --deprecated=PATH    Path to deprecated_recode_mets.rds. Required.
#
# Examples:
#   Rscript scripts/run_batch_process.R curated \
#       --filter='nitrobacter|nitrosomonas' \
#       --ref-data=/data/refs/metanetx_reference_data.rds \
#       --deprecated=reference_data/deprecated_recode_mets.rds
#
#   Rscript scripts/run_batch_process.R flat \
#       --input-dir=carvefungi_species/input \
#       --output-dir=carvefungi_species/processed \
#       --filter='^Aureobasidium' \
#       --ref-data=/data/refs/metanetx_reference_data.rds \
#       --deprecated=reference_data/deprecated_recode_mets.rds

# Resolve repo root without requiring any packages (so --help works on a
# fresh checkout before scripts/install_r_deps.R has been run).
script_path <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(normalizePath(sub("^--file=", "", file_arg[1])))
    }
    # Fallback for source()-from-R use
    normalizePath(sys.frame(1)$ofile %||% "scripts/run_batch_process.R")
})()
`%||%` <- function(a, b) if (is.null(a)) b else a
repo_root <- dirname(dirname(script_path))

parse_arg <- function(args, key, default = NULL) {
    pat <- paste0("^", key, "=")
    hit <- grep(pat, args, value = TRUE)
    if (length(hit) == 0) return(default)
    sub(pat, "", hit[length(hit)])
}

main <- function(args) {
    if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
        # Print the leading comment block as the help text.
        lines <- readLines(script_path, n = 38)
        cat(lines, sep = "\n"); cat("\n")
        quit(status = 0)
    }
    cmd <- args[1]
    rest <- args[-1]

    ref_data_path    <- parse_arg(rest, "--ref-data")
    deprecated_path  <- parse_arg(rest, "--deprecated",
                                  default = file.path(repo_root,
                                                      "reference_data",
                                                      "deprecated_recode_mets.rds"))
    filter_pattern   <- parse_arg(rest, "--filter")

    if (is.null(ref_data_path) || !file.exists(ref_data_path)) {
        stop("--ref-data must point to a readable RDS file. ",
             "See README \"Reference data\" section for download instructions.")
    }
    if (!file.exists(deprecated_path)) {
        stop("--deprecated path does not exist: ", deprecated_path)
    }

    setwd(repo_root)
    source(file.path(repo_root, "pipeline", "process_sbml_species.R"))
    source(file.path(repo_root, "pipeline", "sbml_processing_utils.R"))
    source(file.path(repo_root, "pipeline", "processing_utils.R"))
    source(file.path(repo_root, "pipeline", "batch_process_all.R"))

    cat("Loading reference data from", ref_data_path, "\n")
    ref_data <- readRDS(ref_data_path)
    deprecated_recode <- readRDS(deprecated_path)

    if (cmd == "curated") {
        species_dir <- parse_arg(rest, "--species-dir",
                                 default = file.path(repo_root, "species"))
        if (!dir.exists(species_dir)) {
            stop("species-dir does not exist: ", species_dir)
        }
        species_filter <- NULL
        if (!is.null(filter_pattern)) {
            all_dirs <- basename(list.dirs(species_dir, recursive = FALSE))
            species_filter <- all_dirs[grepl(filter_pattern, all_dirs, perl = TRUE)]
            if (length(species_filter) == 0) {
                stop("No species matched --filter pattern: ", filter_pattern)
            }
            cat("Filter matched", length(species_filter), "species\n")
        }
        results <- batch_process_remote(species_dir, ref_data, deprecated_recode,
                                        species_filter = species_filter)
    } else if (cmd == "flat") {
        input_dir  <- parse_arg(rest, "--input-dir",
                                default = file.path(repo_root, "carvefungi_species", "input"))
        output_dir <- parse_arg(rest, "--output-dir",
                                default = file.path(repo_root, "carvefungi_species", "processed"))
        if (!dir.exists(input_dir)) stop("input-dir does not exist: ", input_dir)

        results <- process_flat_directory(
            input_dir         = input_dir,
            output_dir        = output_dir,
            ref_data          = ref_data,
            deprecated_recode = deprecated_recode,
            file_filter       = filter_pattern
        )
    } else {
        stop("Unknown command: ", cmd, ". Use 'curated' or 'flat'.")
    }

    if (!is.null(results)) {
        cat("\nDone. Processed", nrow(results), "items;",
            sum(results$processing_success, na.rm = TRUE), "successful.\n")
    }
}

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    main(args)
}
