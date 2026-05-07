"""Smoke tests: input + processed SBML for the cohort load and grow as expected.

Each species in tests/golden/smoke_cohort.json is parametrised so failures
identify the species directly. Tolerance for growth-rate equivalence is
1e-4 (matches the golden file's 6-decimal rounding plus solver float noise).
"""
from __future__ import annotations

import math

from conftest import silence_cobra, species_dir

GROWTH_TOL = 1e-4


def _read_growth(cobra, path):
    with silence_cobra():
        m = cobra.io.read_sbml_model(str(path))
        g = float(m.slim_optimize())
    return m, g


def test_input_readable_and_grows(golden_entry, cobra_module):
    name, expected = golden_entry
    inp = species_dir(name) / expected["input_file"]
    assert inp.exists(), f"input SBML missing: {inp}"
    m, g = _read_growth(cobra_module, inp)
    assert len(m.metabolites) == expected["input"]["metabolites"], \
        f"input metabolite count drift for {name}"
    assert len(m.reactions) == expected["input"]["reactions"], \
        f"input reaction count drift for {name}"
    assert math.isclose(g, expected["input"]["growth_rate"], abs_tol=GROWTH_TOL), \
        f"input growth changed for {name}: got {g}, expected {expected['input']['growth_rate']}"


def test_processed_readable_and_grows(golden_entry, cobra_module):
    name, expected = golden_entry
    proc = species_dir(name) / expected["processed_file"]
    assert proc.exists(), f"processed SBML missing: {proc}"
    m, g = _read_growth(cobra_module, proc)
    assert len(m.metabolites) == expected["processed"]["metabolites"], \
        f"processed metabolite count drift for {name}"
    assert len(m.reactions) == expected["processed"]["reactions"], \
        f"processed reaction count drift for {name}"
    assert math.isclose(g, expected["processed"]["growth_rate"], abs_tol=GROWTH_TOL), \
        f"processed growth changed for {name}: got {g}, expected {expected['processed']['growth_rate']}"


def test_input_and_processed_growth_equivalent(golden_entry, cobra_module):
    name, expected = golden_entry
    g_in = expected["input"]["growth_rate"]
    g_pr = expected["processed"]["growth_rate"]
    assert math.isclose(g_in, g_pr, abs_tol=GROWTH_TOL), \
        f"pipeline changed growth for {name}: input {g_in}, processed {g_pr}"
