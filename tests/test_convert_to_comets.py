"""Tests for pipeline/convert_to_comets.py.

Runs the COMETS converter against the smoke cohort and asserts each
produced .cmd file is well-formed (has the expected section headers and
matching SMATRIX dimensions).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

from conftest import REPO_ROOT, GOLDEN, species_dir

# Make pipeline/ importable so we can call convert() directly.
sys.path.insert(0, str(REPO_ROOT / "pipeline"))

import convert_to_comets  # noqa: E402

REQUIRED_SECTIONS = (
    "SMATRIX",
    "BOUNDS",
    "OBJECTIVE",
    "OBJECTIVE_STYLE",
    "OPTIMIZER",
    "EXCHANGE_REACTIONS",
)


@pytest.fixture(params=sorted(GOLDEN.keys()), ids=lambda s: s)
def species_name(request):
    return request.param


def test_convert_produces_well_formed_cmd(species_name, tmp_path):
    expected = GOLDEN[species_name]
    proc = species_dir(species_name) / expected["processed_file"]
    assert proc.exists(), f"processed SBML missing: {proc}"

    out = convert_to_comets.convert(
        input_path=proc,
        output_dir=tmp_path,
        warn_no_biomass=False,
    )
    assert out.exists() and out.stat().st_size > 0, f"empty .cmd output for {species_name}"

    text = out.read_text()
    for section in REQUIRED_SECTIONS:
        assert re.search(rf"^{section}\b", text, re.MULTILINE), (
            f"{species_name}: missing required COMETS section {section!r}"
        )

    # SMATRIX header should declare counts that match the cobra model
    smatrix_header = re.search(r"^SMATRIX\s+(\d+)\s+(\d+)", text, re.MULTILINE)
    assert smatrix_header, f"{species_name}: SMATRIX header malformed"
    n_mets, n_rxns = int(smatrix_header.group(1)), int(smatrix_header.group(2))
    assert n_mets == expected["processed"]["metabolites"], (
        f"{species_name}: SMATRIX metabolite count {n_mets} != "
        f"cobra count {expected['processed']['metabolites']}"
    )
    assert n_rxns == expected["processed"]["reactions"], (
        f"{species_name}: SMATRIX reaction count {n_rxns} != "
        f"cobra count {expected['processed']['reactions']}"
    )


def test_convert_glpk_option(tmp_path):
    """GLPK optimizer override should appear in the output."""
    name = next(iter(GOLDEN))
    proc = species_dir(name) / GOLDEN[name]["processed_file"]
    out = convert_to_comets.convert(
        input_path=proc,
        output_dir=tmp_path,
        warn_no_biomass=False,
        optimizer="GLPK",
    )
    text = out.read_text()
    optimizer_match = re.search(r"^OPTIMIZER\s+(\w+)", text, re.MULTILINE)
    assert optimizer_match, "OPTIMIZER line missing"
    assert optimizer_match.group(1) == "GLPK", (
        f"expected GLPK, got {optimizer_match.group(1)}"
    )
