"""Static-validity tests: regression net for the structural bugs Plan D found.

These run on the already-committed *_processed.xml files without invoking
COBRApy or the pipeline, so they are fast (~2s for the curated cohort).

Bugs guarded against:
  - Duplicate species ids (cesiri, halom, iRP911 family).
  - compartment="NA" (same family).
  - Missing required `id` attribute on species (iJDZ836 family).
  - Invalid SBML SId characters in species/reaction id attributes.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from conftest import REPO_ROOT

# Species whose committed *_processed.xml files still contain known bugs.
#
# The Plan D fix in commit 7fa0d64 (handle_duplicates uniqueness +
# @met_comp recovery) was verified to resolve 7 of these 8 entries
# end-to-end (see issue #5 for the equivalence harness; agent ran the
# pipeline on a stub MetanetX, dedup safety net cleared 265 collisions
# on iRZ1179, NA-compartment recovery cleared 513 on halom, etc.). To
# ACTUALLY remove these entries from KNOWN_BROKEN, the maintainer must
# re-run the pipeline against the real MetanetX reference data and
# commit the regenerated *_processed.xml files. Until then the existing
# committed outputs reflect the pre-fix pipeline and should xfail.
#
# iJDZ836 is the one survivor whose bug is NOT covered by the Plan D
# fix. See issue #7: 5 metal species decode to ids containing `+` and
# `[]` (e.g. `Fe+2[CCO-EXTRACELLULAR]`), which are not legal SBML SIds.
KNOWN_BROKEN = {
    # 7 species: Plan D fix verified to resolve. Remove after re-running
    # the pipeline with real MetanetX data and re-committing the XMLs.
    "cesiribacter_andamanensis_Xu_cesiri",
    "halomonas_stevensii_Xu_halom",
    "hansenula_polymorpha_hanpo_z",
    "methylorubrum_extorquens_iRP911",
    "nitrobacter_winogradskyi_iFC579",
    "pseudomonas_putida_iJN1462",
    "paenarthrobacter_aurescens_iRZ1179",
    # Surviving issue, separate fix needed (see issue #7):
    "neurospora_crassa_iJDZ836",
}

SID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
ID_ATTR_RE = re.compile(r'<(species|reaction)\s[^>]*\bid="([^"]*)"')
SPECIES_NA_COMPARTMENT_RE = re.compile(r'<species\s[^>]*\bcompartment="NA"')
SPECIES_NO_ID_RE = re.compile(r'<species\s(?![^>]*\bid=)[^>]*>')


def _processed_files():
    """Yield (species_dirname, processed_xml_path) tuples."""
    for sp in sorted((REPO_ROOT / "species").glob("*")):
        if not sp.is_dir():
            continue
        for proc in sp.glob("*_processed.xml"):
            yield sp.name, proc


PROCESSED_PARAMS = list(_processed_files())


@pytest.mark.parametrize(
    "species_name,processed_path",
    PROCESSED_PARAMS,
    ids=[name for name, _ in PROCESSED_PARAMS],
)
def test_no_duplicate_species_ids(species_name, processed_path):
    if species_name in KNOWN_BROKEN:
        pytest.xfail(f"known broken: {species_name} (Plan D bug, tracked)")
    text = processed_path.read_text(errors="replace")
    species_ids = re.findall(r'<species\s[^>]*\bid="([^"]+)"', text)
    seen = set()
    dupes = []
    for sid in species_ids:
        if sid in seen:
            dupes.append(sid)
        else:
            seen.add(sid)
    assert not dupes, f"duplicate species ids in {processed_path}: {dupes[:5]}"


@pytest.mark.parametrize(
    "species_name,processed_path",
    PROCESSED_PARAMS,
    ids=[name for name, _ in PROCESSED_PARAMS],
)
def test_no_na_compartment(species_name, processed_path):
    if species_name in KNOWN_BROKEN:
        pytest.xfail(f"known broken: {species_name} (Plan D bug, tracked)")
    text = processed_path.read_text(errors="replace")
    matches = SPECIES_NA_COMPARTMENT_RE.findall(text)
    assert not matches, (
        f"{processed_path.name} contains {len(matches)} species blocks "
        f'with compartment="NA" (invalid SBML, see Plan D notes)'
    )


@pytest.mark.parametrize(
    "species_name,processed_path",
    PROCESSED_PARAMS,
    ids=[name for name, _ in PROCESSED_PARAMS],
)
def test_all_species_have_id(species_name, processed_path):
    if species_name in KNOWN_BROKEN:
        pytest.xfail(f"known broken: {species_name} (Plan D bug, tracked)")
    text = processed_path.read_text(errors="replace")
    bad = SPECIES_NO_ID_RE.findall(text)
    assert not bad, (
        f"{processed_path.name} has {len(bad)} <species> blocks without id"
    )


@pytest.mark.parametrize(
    "species_name,processed_path",
    PROCESSED_PARAMS,
    ids=[name for name, _ in PROCESSED_PARAMS],
)
def test_ids_are_valid_sbml_sids(species_name, processed_path):
    text = processed_path.read_text(errors="replace")
    illegal = []
    for tag, value in ID_ATTR_RE.findall(text):
        if not SID_RE.match(value):
            illegal.append((tag, value))
            if len(illegal) >= 5:  # cap reporting
                break
    assert not illegal, (
        f"{processed_path.name} has {len(illegal)}+ invalid SBML SIds "
        f"(must match [A-Za-z_][A-Za-z0-9_]*): {illegal}"
    )
