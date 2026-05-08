"""Unit tests for repair_fbc_objective() in pipeline/process_sbml_species.R.

Drives the function through an Rscript subprocess against synthetic
modelorg-shaped objects (no real SBML reads needed), so the test is fast
and doesn't depend on either backend being installed.

Coverage
--------
- Detects biomass by react_id pattern, sets objective to 1.
- Detects biomass by react_name when react_id doesn't say "biomass".
- Picks the longest-id biomass when multiple candidates exist.
- Unlocks (0, 0) biomass to (0, 1000) and emits a warning.
- Leaves an existing fbc:objective alone (warns instead of overriding).
- No-op when the backend doesn't expose the required slots.
- No biomass + no objective: emits a single warning.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "pipeline" / "process_sbml_species.R"

if shutil.which("Rscript") is None:
    pytest.skip("Rscript not available", allow_module_level=True)


def _run_repair(scenario_r: str) -> dict:
    """Invoke repair_fbc_objective on a synthetic S4 object built in R.

    `scenario_r` constructs a setClass-shaped object named `m`, then
    the harness calls repair_fbc_objective(m, list()) and prints the
    result as JSON.
    """
    harness = textwrap.dedent(f"""
        suppressPackageStartupMessages({{
            library(methods)
            library(jsonlite)
        }})

        # Define a minimal modelorg shape with the slots repair_fbc_objective uses.
        setClass("modelorg_test",
                 representation(react_id   = "character",
                                react_name = "character",
                                obj_coef   = "numeric",
                                lowbnd     = "numeric",
                                uppbnd     = "numeric"))

        # Source repair_fbc_objective and detect_biomass_reaction by extracting
        # them from the source file (the file's own dependencies are NOT loaded).
        source_lines <- readLines("{SCRIPT}")
        # Find the line range for both functions.
        defs <- which(grepl("^(detect_biomass_reaction|repair_fbc_objective)\\\\s*<-", source_lines))
        last_close <- max(which(grepl("^\\\\}}\\\\s*$", source_lines[defs[1]:length(source_lines)]))) + defs[1] - 1
        # Source just our two function definitions plus any helpers they need.
        # Eval the slice that contains them (from the first def to the second's closing brace).
        end_idx <- defs[length(defs)]
        # Find the closing brace at the lowest indent after end_idx.
        closing <- which(source_lines[end_idx:length(source_lines)] == "}}")[1] + end_idx - 1
        eval(parse(text = paste(source_lines[defs[1]:closing], collapse = "\\n")))

        # Build the scenario object.
        {scenario_r}

        log <- list(warnings = list())
        result <- repair_fbc_objective(m, log)
        out <- list(
            obj_coef          = as.numeric(result$model@obj_coef),
            lowbnd            = as.numeric(result$model@lowbnd),
            uppbnd            = as.numeric(result$model@uppbnd),
            warnings          = result$processing_log$warnings,
            objective_repaired = result$processing_log$objective_repaired,
            biomass_unlocked   = result$processing_log$biomass_unlocked
        )
        cat(toJSON(out, auto_unbox = TRUE, null = "null"))
    """)
    proc = subprocess.run(
        ["Rscript", "--vanilla", "-e", harness],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0, f"Rscript failed:\n{proc.stderr}"
    # The function uses cat() for human-readable warnings; parse only the trailing JSON.
    out = proc.stdout.strip()
    # Find the JSON object — last {...} block
    last_brace_open = out.rfind("{")
    json_block = out[last_brace_open:]
    return json.loads(json_block)


def test_detect_biomass_by_react_id_and_set_objective():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("rxn001", "BIOMASS_core", "rxn003"),
                 react_name = c("",       "",             ""),
                 obj_coef   = c(0, 0, 0),
                 lowbnd     = c(-1000, 0, -1000),
                 uppbnd     = c( 1000, 1000, 1000))
    """)
    assert r["obj_coef"] == [0.0, 1.0, 0.0]
    assert r["objective_repaired"] == "BIOMASS_core"
    assert r.get("biomass_unlocked") in (None, "null")


def test_detect_biomass_by_react_name_when_id_does_not_say():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("rxn001", "rxn002", "rxn003"),
                 react_name = c("",       "Biomass production reaction", ""),
                 obj_coef   = c(0, 0, 0),
                 lowbnd     = c(-1000, 0, -1000),
                 uppbnd     = c( 1000, 1000, 1000))
    """)
    assert r["obj_coef"] == [0.0, 1.0, 0.0]


def test_pick_longest_id_when_multiple_biomass_candidates():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("biomass0", "BIOMASS_Gm_GS15_core_79p20M", "rxn"),
                 react_name = c("",         "",                            ""),
                 obj_coef   = c(0, 0, 0),
                 lowbnd     = c(0, 0, -1000),
                 uppbnd     = c(1000, 1000, 1000))
    """)
    # Longest id wins.
    assert r["objective_repaired"] == "BIOMASS_Gm_GS15_core_79p20M"
    assert r["obj_coef"] == [0.0, 1.0, 0.0]


def test_unlock_locked_biomass():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("rxn001", "BIOMASS_locked", "rxn003"),
                 react_name = c("", "", ""),
                 obj_coef   = c(0, 0, 0),
                 lowbnd     = c(-1000, 0, -1000),
                 uppbnd     = c( 1000, 0, 1000))
    """)
    assert r["uppbnd"] == [1000.0, 1000.0, 1000.0], "biomass upper bound should unlock"
    assert r["lowbnd"] == [-1000.0, 0.0, -1000.0],  "biomass lower bound stays at 0"
    assert r["biomass_unlocked"] == "BIOMASS_locked"


def test_existing_objective_left_alone_with_warning():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("rxn001", "BIOMASS_core", "ATP_maintenance"),
                 react_name = c("", "", ""),
                 obj_coef   = c(0, 0, 1),  # objective already on ATP_maintenance, NOT biomass
                 lowbnd     = c(-1000, 0, 0),
                 uppbnd     = c( 1000, 1000, 1000))
    """)
    # Plan: existing objective preserved (not overwritten); a warning emitted.
    assert r["obj_coef"] == [0.0, 0.0, 1.0]
    assert any("BIOMASS_core" in w and "ATP_maintenance" in w for w in r["warnings"]), (
        f"expected a warning naming both biomass and the current objective; got {r['warnings']}"
    )


def test_no_biomass_and_no_objective_emits_one_warning():
    r = _run_repair("""
        m <- new("modelorg_test",
                 react_id   = c("rxn001", "rxn002"),
                 react_name = c("", ""),
                 obj_coef   = c(0, 0),
                 lowbnd     = c(-1000, -1000),
                 uppbnd     = c( 1000, 1000))
    """)
    assert r["obj_coef"] == [0.0, 0.0]
    assert any("growth" in w.lower() and "0" in w for w in r["warnings"]), (
        f"expected a warning that growth will be 0; got {r['warnings']}"
    )


def test_no_op_when_required_slots_are_missing():
    """If the backend doesn't expose obj_coef / lowbnd / uppbnd, the function
    should leave the model unchanged and not error out."""
    harness = textwrap.dedent(f"""
        suppressPackageStartupMessages({{
            library(methods); library(jsonlite)
        }})
        setClass("modelorg_thin",
                 representation(react_id = "character", react_name = "character"))
        source_lines <- readLines("{SCRIPT}")
        defs <- which(grepl("^(detect_biomass_reaction|repair_fbc_objective)\\\\s*<-", source_lines))
        end_idx <- defs[length(defs)]
        closing <- which(source_lines[end_idx:length(source_lines)] == "}}")[1] + end_idx - 1
        eval(parse(text = paste(source_lines[defs[1]:closing], collapse = "\\n")))

        m <- new("modelorg_thin",
                 react_id   = c("rxn001", "BIOMASS_core"),
                 react_name = c("", ""))
        log <- list(warnings = list())
        result <- repair_fbc_objective(m, log)
        # If we got here, the function bailed cleanly. Confirm the returned
        # model is the same object we passed in.
        cat(toJSON(list(
            ok = TRUE,
            same_class = identical(class(result$model)[1], "modelorg_thin"),
            warnings = result$processing_log$warnings
        ), auto_unbox = TRUE, null = "null"))
    """)
    proc = subprocess.run(
        ["Rscript", "--vanilla", "-e", harness],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0, f"Rscript failed:\n{proc.stderr}"
    out = proc.stdout.strip()
    last_brace_open = out.rfind("{")
    r = json.loads(out[last_brace_open:])
    assert r["ok"] is True
    assert r["same_class"] is True
    # No warnings should have been emitted (early return path).
    assert r["warnings"] in ([], None, "null", [], {})
