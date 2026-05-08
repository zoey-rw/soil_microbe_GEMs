#!/usr/bin/env bash
# Fetch the MetaNetX reference TSVs and build the cached
# reference_data/metanetx_reference_data.rds that the pipeline reads.
#
# Usage:
#   bash scripts/fetch_metanetx_refdata.sh           # default (most recent)
#   MNX_VERSION=4.5 bash scripts/fetch_metanetx_refdata.sh   # pin a release
#
# Notes:
# - These files are gitignored due to size (~150 MB total uncompressed).
# - metanetx.org may rate-limit; if so, retry with a delay.
# - The script is idempotent: skips files already present unless FORCE=1.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$REPO_ROOT/reference_data"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

MNX_VERSION="${MNX_VERSION:-}"
FORCE="${FORCE:-0}"

if [[ -n "$MNX_VERSION" ]]; then
    BASE="https://www.metanetx.org/ftp/${MNX_VERSION}"
    echo "Using pinned MetaNetX version: $MNX_VERSION"
else
    BASE="https://www.metanetx.org/cgi-bin/mnxget/mnxref"
    echo "Using latest MetaNetX (unpinned). Set MNX_VERSION=4.5 (or similar) for reproducibility."
fi

FILES=(
    chem_xref.tsv
    chem_prop.tsv
    reac_xref.tsv
    reac_prop.tsv
)

for f in "${FILES[@]}"; do
    if [[ -s "$f" && "$FORCE" != "1" ]]; then
        echo "  $f exists ($(stat -c%s "$f") bytes); skipping (FORCE=1 to re-download)"
        continue
    fi
    echo "  fetching $f from $BASE/$f ..."
    if ! curl -fL --retry 3 --retry-delay 5 -o "$f.tmp" "$BASE/$f"; then
        echo "  ERROR: failed to fetch $f. Check network access to metanetx.org."
        rm -f "$f.tmp"
        exit 1
    fi
    mv "$f.tmp" "$f"
    echo "    -> $(stat -c%s "$f") bytes"
done

echo
echo "All TSVs in place. Building cached metanetx_reference_data.rds ..."
echo "(this calls pipeline/sbml_processing_utils.R::get_reference_data())"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found. Install R first (apt install r-base-core)."
    exit 1
fi

cd "$REPO_ROOT"
Rscript -e '
suppressPackageStartupMessages({
    library(here)
})
source(here::here("pipeline", "sbml_processing_utils.R"))
cat("Calling get_reference_data() ...\n")
ref <- get_reference_data()
out <- here::here("reference_data", "metanetx_reference_data.rds")
saveRDS(ref, out)
cat(sprintf("Wrote %s\n", out))
cat(sprintf("  chem_xref: %d rows, chem_prop: %d rows\n",
            nrow(ref$chem_xref), nrow(ref$chem_prop)))
'

echo
echo "Done. The pipeline can now use:"
echo "  Rscript scripts/run_batch_process.R curated \\"
echo "    --ref-data=reference_data/metanetx_reference_data.rds \\"
echo "    --deprecated=reference_data/deprecated_recode_mets.rds \\"
echo "    --filter=<pattern>"
