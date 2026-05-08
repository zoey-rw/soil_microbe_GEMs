#' Cobra-via-reticulate SBML I/O shim.
#'
#' Drop-in replacement for the small slice of sybilSBML that the pipeline
#' actually uses: readSBMLmod(), writeSBML(), findExchReact(), and the
#' modelorg slot accessors (@met_id, @met_name, @met_comp, @met_attr,
#' @mod_compart, @react_id, @react_rev, @gpr).
#'
#' Backed by pipeline/sbml_io/__init__.py via the reticulate package.
#' This shim is transitional; see issue #2 for the eventual full Python port.
#'
#' To enable, source this file *instead of* `library(sybilSBML)` at the top
#' of pipeline/process_sbml_species.R, processing_utils.R, and
#' check_exchange_metabolites.R. Or set SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1
#' before sourcing those files (see the conditional pattern at the bottom).

suppressPackageStartupMessages({
    library(reticulate)
    library(methods)
})

# ----------------------------------------------------------------------------
# S4 classes mirroring sybilSBML's modelorg / exchReact for slot compatibility.
# ----------------------------------------------------------------------------

setClass(
    "modelorg_cobra",
    representation(
        handle      = "integer",
        met_id      = "character",
        met_name    = "character",
        met_comp    = "integer",
        met_attr    = "list",       # holds $annotation = character vector
        mod_compart = "character",
        react_id    = "character",
        react_name  = "character",
        react_rev   = "logical",
        obj_coef    = "numeric",    # per-reaction objective coefficient
        lowbnd      = "numeric",    # per-reaction lower bound
        uppbnd      = "numeric",    # per-reaction upper bound
        gpr         = "character",
        source_file = "character"
    )
)

setClass(
    "exchReact_cobra",
    representation(
        react_id    = "character",
        met_id      = "character",
        uptake      = "logical",
        lower_bound = "numeric",
        upper_bound = "numeric"
    )
)

# Soft-define a "modelorg" virtual class so `is(x, "modelorg")` checks pass
# without requiring sybilSBML to be loaded.
if (!isClass("modelorg")) {
    setClass("modelorg", representation("VIRTUAL"))
}
setIs("modelorg_cobra", "modelorg",
      coerce  = function(from) from,
      replace = function(from, value) value)

# ----------------------------------------------------------------------------
# Lazy-loaded Python module.
# ----------------------------------------------------------------------------

.sbml_io_module <- NULL

.load_sbml_io <- function() {
    if (is.null(.sbml_io_module)) {
        repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                                   mustWork = FALSE)
        if (!dir.exists(file.path(repo_root, "pipeline", "sbml_io"))) {
            # Fallback: assume CWD is repo root
            repo_root <- getwd()
        }
        reticulate::use_python(Sys.which("python3"), required = FALSE)
        py_path <- file.path(repo_root, "pipeline")
        # Prepend pipeline/ to sys.path so `import sbml_io` resolves
        sys <- reticulate::import("sys", convert = FALSE)
        if (!(py_path %in% reticulate::py_to_r(sys$path))) {
            sys$path$insert(0L, py_path)
        }
        .sbml_io_module <<- reticulate::import("sbml_io", convert = TRUE)
    }
    .sbml_io_module
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# Drop-in replacements for sybilSBML public API.
# ----------------------------------------------------------------------------

#' Read an SBML model via cobra. Signature accepts (and ignores) the
#' sybilSBML-specific tuning flags so callers don't need to change their
#' do.call() invocations.
readSBMLmod <- function(file,
                        validateSBML = TRUE,
                        bndCond      = TRUE,
                        mergeMet     = TRUE,
                        balanceReact = TRUE,
                        def_bnd      = 1000,
                        ...) {
    mod <- .load_sbml_io()
    state <- mod$read_sbml(file)

    obj <- new(
        "modelorg_cobra",
        handle      = as.integer(state$handle),
        met_id      = as.character(state$met_id),
        met_name    = as.character(state$met_name),
        met_comp    = as.integer(state$met_comp),
        met_attr    = list(annotation = as.character(state$met_annotation)),
        mod_compart = as.character(state$mod_compart),
        react_id    = as.character(state$react_id),
        react_name  = as.character(state$react_name %||% state$react_id),
        react_rev   = as.logical(state$react_rev),
        obj_coef    = as.numeric(state$obj_coef %||% rep(0, length(state$react_id))),
        lowbnd      = as.numeric(state$lowbnd   %||% rep(-1000, length(state$react_id))),
        uppbnd      = as.numeric(state$uppbnd   %||% rep( 1000, length(state$react_id))),
        gpr         = as.character(state$gpr),
        source_file = as.character(file)
    )

    # Note: no GC finalizer attached — S4 objects can't accept reg.finalizer.
    # The Python-side cache grows during a batch run; call sbml_io_clear()
    # explicitly between large batches if memory matters.
    obj
}

#' Drop all cached cobra models. Use between large batches if memory is a
#' concern; otherwise the cache is freed when R exits.
sbml_io_clear <- function() {
    mod <- .load_sbml_io()
    n_before <- mod$cache_size()
    py <- reticulate::import("sbml_io", convert = TRUE)
    # Directly clear the cache by releasing each handle we know about.
    sys <- reticulate::import("sys", convert = FALSE)
    cache <- reticulate::py_get_attr(py, "_CACHE")
    keys <- names(reticulate::py_to_r(cache))
    if (!is.null(keys)) {
        for (k in keys) py$release(as.integer(k))
    }
    cat("sbml_io cache: cleared", n_before, "models\n")
    invisible(n_before)
}

#' Write a (possibly mutated) model back to SBML.
writeSBML <- function(model, level = 3, filename, ...) {
    if (!inherits(model, "modelorg_cobra")) {
        stop("writeSBML (cobra shim) expects a modelorg_cobra; got ", class(model))
    }
    mod <- .load_sbml_io()
    state <- list(
        handle   = as.integer(model@handle),
        met_id   = as.character(model@met_id),
        react_id = as.character(model@react_id),
        gpr      = as.character(model@gpr),
        obj_coef = as.numeric(model@obj_coef),
        lowbnd   = as.numeric(model@lowbnd),
        uppbnd   = as.numeric(model@uppbnd)
    )
    res <- mod$write_sbml(state, filename, as.integer(level))
    invisible(res)
}

#' Find exchange reactions, returning an exchReact-compatible S4 object.
findExchReact <- function(model) {
    if (!inherits(model, "modelorg_cobra")) {
        stop("findExchReact (cobra shim) expects a modelorg_cobra; got ", class(model))
    }
    mod <- .load_sbml_io()
    state <- list(handle = as.integer(model@handle))
    ex <- mod$exchanges(state)
    new(
        "exchReact_cobra",
        react_id    = as.character(ex$react_id),
        met_id      = as.character(ex$met_id),
        uptake      = as.logical(ex$uptake),
        lower_bound = as.numeric(ex$lower_bound),
        upper_bound = as.numeric(ex$upper_bound)
    )
}
