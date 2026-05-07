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

# Species with KNOWN-BAD processed files. These are the regressions Plan D
# identified, plus 3 more (iFC579, iJN1462, iRZ1179) the harness uncovered
# with finer granularity than Plan D's manual analysis. All are the same
# root cause: the species-merge / "additional reactions" code path emits
# duplicate ids, id-less species, or compartment="NA". They will be removed
# from this list as Plan D's pipeline fix lands.
KNOWN_BROKEN = {
    # Original Plan D findings (unreadable by cobra after processing)
    "cesiribacter_andamanensis_Xu_cesiri",
    "halomonas_stevensii_Xu_halom",
    "hansenula_polymorpha_hanpo_z",
    "methylorubrum_extorquens_iRP911",
    "neurospora_crassa_iJDZ836",
    # Discovered by this test harness (still readable but malformed)
    "nitrobacter_winogradskyi_iFC579",       # 1x compartment="NA"
    "pseudomonas_putida_iJN1462",            # 2x compartment="NA"
    "paenarthrobacter_aurescens_iRZ1179",    # 265 duplicate species ids
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
