#!/usr/bin/env bash
# Fetch the MetaNetX reference TSVs and build the cached
# reference_data/metanetx_reference_data.rds that the pipeline reads.
#
# Usage:
#   bash scripts/fetch_metanetx_refdata.sh                     # default
#   MNX_VERSION=4.5 bash scripts/fetch_metanetx_refdata.sh     # pin a release
#   MNX_RELEASE_TAG=metanetx-4.5 bash scripts/fetch_metanetx_refdata.sh
#                                                              # use a github
#                                                              # release as the
#                                                              # source instead
#                                                              # of metanetx.org
#
# Notes:
# - The TSVs are gitignored due to size (~150 MB total).
# - Sources tried, in order:
#     1. github release on $MNX_RELEASE_REPO (default zoey-rw/soil_microbe_GEMs)
#        if MNX_RELEASE_TAG is set. Use this in sandboxed environments where
#        metanetx.org is firewalled but github.com is reachable. To create the
#        release, upload chem_xref.tsv, chem_prop.tsv, reac_xref.tsv,
#        reac_prop.tsv as assets on a release named e.g. "metanetx-4.5".
#     2. metanetx.org/cgi-bin/mnxget/mnxref/ (or /ftp/$MNX_VERSION/ if pinned).
# - The script is idempotent: skips files already present unless FORCE=1.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$REPO_ROOT/reference_data"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

MNX_VERSION="${MNX_VERSION:-}"
MNX_RELEASE_TAG="${MNX_RELEASE_TAG:-}"
MNX_RELEASE_REPO="${MNX_RELEASE_REPO:-zoey-rw/soil_microbe_GEMs}"
FORCE="${FORCE:-0}"

if [[ -n "$MNX_RELEASE_TAG" ]]; then
    BASE="https://github.com/${MNX_RELEASE_REPO}/releases/download/${MNX_RELEASE_TAG}"
    SOURCE_LABEL="github release ${MNX_RELEASE_REPO}@${MNX_RELEASE_TAG}"
elif [[ -n "$MNX_VERSION" ]]; then
    BASE="https://www.metanetx.org/ftp/${MNX_VERSION}"
    SOURCE_LABEL="MetaNetX FTP, version ${MNX_VERSION}"
else
    BASE="https://www.metanetx.org/cgi-bin/mnxget/mnxref"
    SOURCE_LABEL="MetaNetX (latest, unpinned)"
fi
echo "Source: $SOURCE_LABEL"

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
        echo "  ERROR: failed to fetch $f from $SOURCE_LABEL."
        if [[ -z "$MNX_RELEASE_TAG" ]]; then
            echo "  In sandboxed environments where metanetx.org is firewalled,"
            echo "  upload the four TSVs as assets on a github release on"
            echo "  $MNX_RELEASE_REPO and re-run with MNX_RELEASE_TAG=<tag>."
        fi
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
