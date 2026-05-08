# Contributing & maintainer notes

This is a working notebook for contributors and the project lead. It catalogues
what changed in the May 2026 cleanup pass, documents the env vars and entry
points that landed, and lists the deferred items that are tracked under
GitHub issues.

## Backends

The pipeline reads + writes SBML through one of two interchangeable backends:

| Backend | How to enable | When to use |
|---|---|---|
| Patched native sybilSBML | (default) — install via `bash pipeline/vendor/install_sybilSBML.sh` | Day-to-day. Faster, byte-identical with historical outputs. |
| Cobra shim (reticulate) | `export SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1` | When sybilSBML can't compile (modern libSBML), or as an A/B against the patched one. Closer to the eventual full Python port. |

The selection happens in `pipeline/sbml_io_loader.R`, which is sourced by
`process_sbml_species.R`, `processing_utils.R`, and `check_exchange_metabolites.R`
in place of the old `library(sybilSBML)` line. To swap permanently in either
direction, edit the loader; the call sites stay unchanged.

## Environment variables

| Variable | Purpose | Default |
|---|---|---|
| `SOIL_MICROBE_GEMS_USE_COBRA_SHIM` | Switch backend (`1` = cobra shim) | unset = patched sybilSBML |
| `SOIL_MICROBE_GEMS_HAVE_SYBILSBML` | Tells `tests/test_shim_equivalence.py` whether to run | `1` locally; `0` on stock GHA |
| `SOIL_MICROBE_GEMS_VALIDATE_TIMEOUT_SEC` | Cobra-validation timeout per model | 600 |
| `SOIL_MICROBE_DB` | Path to the upstream `soil_microbe_db` data dir (used by `comets_shinyapp_example/01_merge_env_data.r`) | unset (script errors with a pointer) |
| `MNX_VERSION` | MetaNetX release pin for `scripts/fetch_metanetx_refdata.sh` | unset = latest |
| `MNX_RELEASE_TAG` / `MNX_RELEASE_REPO` | GitHub-release fallback for sandboxed environments where `metanetx.org` is firewalled | unset |
| `POSTGRES_USER` / `POSTGRES_PASSWORD` / `POSTGRES_DB` | Shiny app's Postgres creds (read from `.env`) | none — startup fails if `POSTGRES_PASSWORD` is unset |

## Entry points

| Command | Purpose |
|---|---|
| `Rscript scripts/install_r_deps.R` | Install all R dependencies including the patched sybilSBML |
| `bash scripts/fetch_metanetx_refdata.sh` | Download MetaNetX TSVs and build `metanetx_reference_data.rds` |
| `Rscript scripts/run_batch_process.R --help` | Process species (curated or flat) end-to-end |
| `Rscript pipeline/batch_validate_all.R` | Run COBRApy validation across all processed species |
| `python pipeline/reproduce_validation.py` | Independent Python validation, writes `pipeline/reproduction_report.json` |
| `python pipeline/convert_to_comets.py <input.xml>` | Convert one processed SBML → COMETS `.cmd` |
| `pytest tests/` | Run the regression harness |
| `python tests/regenerate_golden.py` | Refresh `tests/golden/smoke_cohort.json` after a deliberate pipeline change |

## Tests

- `tests/test_pipeline_smoke.py` — five fast species, asserts read+growth equivalence to `tests/golden/smoke_cohort.json`.
- `tests/test_processed_validity.py` — static-XML regression net for the structural bug Plan D found. Maintains a `KNOWN_BROKEN` set of species whose committed `*_processed.xml` still has bugs (xfailed). Removing a name from the set is the signal that the underlying fix has been verified end-to-end.
- `tests/test_convert_to_comets.py` — COMETS converter tests on the smoke cohort.
- `tests/test_shim_equivalence.py` — slot-by-slot diff between patched sybilSBML and the cobra-shim backend (issue #5). Skips automatically on environments without sybilSBML; gate via `SOIL_MICROBE_GEMS_HAVE_SYBILSBML`.
- `tests/test_repair_fbc_objective.py` — synthetic-input coverage of `repair_fbc_objective()` (issue #6). All six policy branches plus the no-op case. Doesn't require either backend installed.
- `tests/test_sid_reencode.py` — synthetic-input coverage of `enforce_sbml_sids()` (issue #7). Doesn't require either backend.

CI:

- `.github/workflows/test.yml` — runs `pytest tests/` on every push / PR.
- `.github/workflows/full-validation.yml` — manually triggered (`workflow_dispatch`); runs `pipeline/reproduce_validation.py` across the full curated cohort and uploads the JSON report.

## Open issues (and what blocks them)

| # | Title | Status |
|---|---|---|
| [#2](https://github.com/zoey-rw/soil_microbe_GEMs/issues/2) | Plan a full Python port | Tracking issue. Stage 1 (shim) landed. Stage 2 (full port) blocked on the equivalence harness in #5. |
| [#3](https://github.com/zoey-rw/soil_microbe_GEMs/issues/3) | shinyapp scripts 02-04 still depend on hardcoded CWD | Untouched. Blocked on access to the upstream `soil_microbe_db` dataset. |
| [#4](https://github.com/zoey-rw/soil_microbe_GEMs/issues/4) | Verify the `merged_df` alias in `01_merge_env_data.r` | Workaround alias landed; verification pending end-to-end run with real data. |
| [#5](https://github.com/zoey-rw/soil_microbe_GEMs/issues/5) | Per-species shim equivalence harness | Landed (`tests/test_shim_equivalence.py`). Local: 34 passed / 3 skipped / 3 xfailed for documented divergences. |
| [#6](https://github.com/zoey-rw/soil_microbe_GEMs/issues/6) | Translate legacy `OBJECTIVE_COEFFICIENT` to `fbc:fluxObjective` | Landed (`repair_fbc_objective()` in `process_sbml_species.R`; 7 unit tests). End-to-end verification on iAF987 still pending real-MetaNetX run. |
| [#7](https://github.com/zoey-rw/soil_microbe_GEMs/issues/7) | Re-encode SBML-SId-illegal characters before writeSBML | Landed (`enforce_sbml_sids()` in `process_sbml_species.R`; 5 unit tests). End-to-end verification on iJDZ836 still pending real-MetaNetX run. |

Bound clamping (Plan D's third bug) is intentionally not filed — per the project lead, the eventual fix must be configurable per species via `config/species_registry.yaml`. Until that design is locked in, the current ±1000 clamping behaviour is unchanged.

## Maintainer action checklist (one-time, on a non-sandboxed machine)

To finish closing Plan D:

```bash
git pull
bash scripts/fetch_metanetx_refdata.sh           # downloads ~150 MB
Rscript scripts/run_batch_process.R curated \
    --filter='^(cesiribacter|halomonas|hansenula|methylorubrum|nitrobacter_winogradskyi|pseudomonas_putida|paenarthrobacter)' \
    --ref-data=reference_data/metanetx_reference_data.rds \
    --deprecated=reference_data/deprecated_recode_mets.rds
git add species/    # review the regenerated *_processed.xml first
git commit -m "Regenerate processed XMLs with Plan D fix"

# Then in tests/test_processed_validity.py, remove the 7 verified species
# from KNOWN_BROKEN (cesiribacter, halomonas, hansenula, methylorubrum,
# nitrobacter_winogradskyi, pseudomonas_putida, paenarthrobacter).
pytest tests/   # confirm only neurospora_crassa_iJDZ836 xfails
```

Aureobasidium re-run (35 unprocessed CarveFungi inputs):

```bash
Rscript scripts/run_batch_process.R flat \
    --input-dir=carvefungi_species/input \
    --output-dir=carvefungi_species/processed \
    --filter='^Aureobasidium' \
    --ref-data=reference_data/metanetx_reference_data.rds \
    --deprecated=reference_data/deprecated_recode_mets.rds
```

## Layout reminder

See README.md → "Repository layout" for the full table of paths and what
lives where.
