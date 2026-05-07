# Vendored installer for sybilSBML

This directory contains a self-contained installer for `sybilSBML 3.1.2` that
patches the upstream tarball to compile against modern libSBML (5.20+).

`sybilSBML` was archived from CRAN in 2020. Without these patches, the
upstream source fails on libSBML ≥ 5.19 with two issues:

1. `unknown type name 'FbcVariableType_t'` — libSBML's `FluxObjective.h`
   started referencing this enum unconditionally; it lives in
   `<sbml/packages/fbc/extension/FbcExtension.h>` which sybilSBML doesn't
   pull in.
2. `incompatible pointer types` — the fbc-plugin C-API's species/reaction
   setters now take `FbcSBasePlugin_t *` instead of `SBasePlugin_t *`.

`sybilSBML.patch` adds the missing include and inserts explicit
`(FbcSBasePlugin_t *)` casts at the relevant call sites in `src/sybilSBML.c`.
It does not change any other behavior. Verified on R 4.3.3 + libSBML 5.20.2:
loads cleanly and reads `species/nitrobacter_winogradskyi_iFC579/iFC579_input.xml`
(1007 metabolites, 1128 reactions).

## Usage

```bash
# From the repo root:
bash pipeline/vendor/install_sybilSBML.sh
```

Requires `R`, `libsbml5-dev`, `libglpk-dev`, plus internet access to fetch
the source tarball. Falls back from the CRAN archive URL to the
`github.com/cran/sybilSBML` mirror if needed.

## Relationship to the cobra shim

The repo also ships a `library(sybilSBML)`-compatible alternative built
on top of cobrapy (`pipeline/sbml_io_cobra.R`, enabled via
`SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1`). Both back-ends are supported during
the transition window (see issue #2 and #5). Pick whichever is easier in
your environment:

- The patched native sybilSBML is faster and handles the existing test
  cohort byte-identically (it was the original implementation).
- The cobra shim is portable to systems where libSBML/glpk linkage is
  hard, and is on the path to the eventual pure-Python pipeline.
