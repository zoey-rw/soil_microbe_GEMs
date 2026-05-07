"""Shared fixtures and helpers for pipeline regression tests.

The smoke cohort is a small set of fast-loading species committed to
tests/golden/smoke_cohort.json. Re-generate it with:

    python tests/regenerate_golden.py
"""
from __future__ import annotations

import contextlib
import json
import logging
import os
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN = json.loads((REPO_ROOT / "tests" / "golden" / "smoke_cohort.json").read_text())


@contextlib.contextmanager
def silence_cobra():
    """COBRApy emits a flood of XML/SBML warnings on import + read.
    Suppress them at the test boundary so failures stay readable.
    """
    warnings.filterwarnings("ignore")
    logging.getLogger("cobra").setLevel(logging.CRITICAL)
    with open(os.devnull, "w") as devnull, contextlib.redirect_stderr(devnull):
        yield


@pytest.fixture(scope="session")
def cobra_module():
    with silence_cobra():
        import cobra  # noqa: WPS433  (deferred import is intentional)
    return cobra


@pytest.fixture(params=sorted(GOLDEN.keys()), ids=lambda s: s)
def golden_entry(request):
    name = request.param
    return name, GOLDEN[name]


def species_dir(name: str) -> Path:
    return REPO_ROOT / "species" / name
