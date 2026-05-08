# Archived scripts

Historical / deprecated code, kept in git but not part of the active pipeline.

| File | Status | Replacement |
|---|---|---|
| `source_functions.R` | Superseded by `pipeline/sbml_processing_utils.R` | Use the active pipeline file. |
| `create_media_database.r` | Reference-only — uses paths that no longer exist | No direct replacement; was a one-off media-curation helper. |
| `evaluate_carbon_usage.ipynb` | Reference notebook | Hardcoded `/projectnb2/...` paths and depends on `memote` (now in `requirements.txt`). Use as a starting template if you need carbon-usage scans; will not run as-is. |

These files are kept for traceability — searching git history is not always the
fastest way to find a previous approach. **Nothing here is on the path of any
current entry point** (`scripts/run_batch_process.R`, `pipeline/*.R`, the
Shiny app, or the test harness). Feel free to delete an entry once you're
confident its replacement is sufficient.
