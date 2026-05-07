#!/usr/bin/env python3
"""Reproduce the R pipeline's validation results in pure Python.

For every species directory under <repo>/species/:
  * locate the input SBML (any *.xml or *.sbml that is not *_processed.xml)
  * locate the processed SBML (*_processed.xml)
  * try cobra.io.read_sbml_model on each (capturing warnings/exceptions)
  * for any readable model, run slim_optimize() and capture growth, n_mets, n_rxns
  * compare against the existing per-species validation_results.json
  * emit a fresh JSON report at pipeline/reproduction_report.json

This script is read-only: it does not modify anything in the species/
directories or the R pipeline. It only writes the final JSON report.
"""

from __future__ import annotations

import json
import os
import sys
import time
import traceback
import warnings
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
SPECIES_DIR = REPO_ROOT / "species"
OUT_PATH = REPO_ROOT / "pipeline" / "reproduction_report.json"
GROWTH_TOL = 1e-6


def _find_xml_pair(species_path: Path) -> tuple[Path | None, Path | None]:
    """Return (input_xml, processed_xml) inside a species directory.

    The R pipeline writes <id>_processed.xml. The input may be named anything
    else with an .xml or .sbml extension. Heuristics:
      * processed = first file matching *_processed.xml
      * input = the first .xml/.sbml that is NOT *_processed.xml; we prefer
        names ending with _input or matching the modelid in processed (the
        prefix before _processed.xml), then anything else, deterministic by
        sorted name to keep results reproducible.
    """
    candidates = [
        p
        for p in species_path.iterdir()
        if p.is_file() and p.suffix.lower() in {".xml", ".sbml"}
    ]
    processed = next((p for p in candidates if p.name.endswith("_processed.xml")), None)

    non_processed = sorted(
        [p for p in candidates if not p.name.endswith("_processed.xml")],
        key=lambda p: p.name,
    )

    # Heuristic ranking
    def rank(p: Path) -> tuple[int, str]:
        name = p.name
        if name.endswith("_input.xml"):
            return (0, name)
        if processed is not None:
            modelid = processed.name[: -len("_processed.xml")]
            # Prefer files whose stem == modelid (e.g. iAF987.xml)
            if p.stem == modelid:
                return (1, name)
        if name.endswith("_modified.xml"):
            return (3, name)
        return (2, name)

    non_processed.sort(key=rank)
    inp = non_processed[0] if non_processed else None
    return inp, processed


def _try_read(path: Path) -> dict[str, Any]:
    """Attempt to read an SBML file. Return a dict with status + metrics."""
    record: dict[str, Any] = {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else None,
        "readable": False,
        "error": None,
        "warnings": [],
        "metabolites": None,
        "reactions": None,
        "growth_rate": None,
        "objective": None,
        "read_seconds": None,
        "optimize_seconds": None,
    }
    if not path.exists():
        record["error"] = "file does not exist"
        return record

    # Late import so cobra optimize problems propagate cleanly
    import cobra  # noqa: WPS433

    t0 = time.time()
    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            model = cobra.io.read_sbml_model(str(path))
        record["read_seconds"] = round(time.time() - t0, 3)
        record["warnings"] = [
            f"{w.category.__name__}: {str(w.message)[:300]}" for w in ws
        ]
        record["readable"] = True
        record["metabolites"] = len(model.metabolites)
        record["reactions"] = len(model.reactions)
        try:
            obj_expr = model.objective.expression
            record["objective"] = str(obj_expr)[:300]
        except Exception as exc:  # noqa: BLE001
            record["objective"] = f"<error reading objective: {exc}>"
    except Exception as exc:  # noqa: BLE001
        record["read_seconds"] = round(time.time() - t0, 3)
        record["error"] = f"{type(exc).__name__}: {exc}"
        return record

    t1 = time.time()
    try:
        sol = model.slim_optimize(error_value=None)
        record["growth_rate"] = (
            float(sol) if sol is not None and sol == sol else None  # NaN check
        )
    except Exception as exc:  # noqa: BLE001
        record["growth_rate"] = None
        record["error_optimize"] = f"{type(exc).__name__}: {exc}"
    record["optimize_seconds"] = round(time.time() - t1, 3)
    return record


def _load_old_validation(species_path: Path) -> dict[str, Any] | None:
    p = species_path / "validation_results.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception as exc:  # noqa: BLE001
        return {"_load_error": f"{type(exc).__name__}: {exc}"}


def _approx(a: float | None, b: float | None, tol: float = GROWTH_TOL) -> bool:
    if a is None or b is None:
        return False
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return False


def _compare(old: dict[str, Any] | None, new_input: dict[str, Any], new_proc: dict[str, Any]) -> dict[str, Any]:
    """Compare new run with the previously recorded validation_results.json."""
    cmp: dict[str, Any] = {
        "old_present": old is not None,
        "old_input_readable": None,
        "old_processed_readable": None,
        "old_growth_input": None,
        "old_growth_processed": None,
        "input_readable_match": None,
        "processed_readable_match": None,
        "input_growth_match": None,
        "processed_growth_match": None,
        "growth_diff_input": None,
        "growth_diff_processed": None,
        "old_input_metab": None,
        "old_input_rxn": None,
        "old_proc_metab": None,
        "old_proc_rxn": None,
        "input_size_match": None,
        "processed_size_match": None,
    }
    if not old:
        return cmp

    cmp["old_input_readable"] = old.get("input_readable_cobra")
    cmp["old_processed_readable"] = old.get("processed_readable_cobra")

    def _num(v: Any) -> float | None:
        if isinstance(v, (int, float)):
            return float(v)
        return None

    cmp["old_growth_input"] = _num(old.get("growth_rate_input"))
    cmp["old_growth_processed"] = _num(old.get("growth_rate_processed"))

    cmp["input_readable_match"] = (
        bool(cmp["old_input_readable"]) == bool(new_input.get("readable"))
    )
    cmp["processed_readable_match"] = (
        bool(cmp["old_processed_readable"]) == bool(new_proc.get("readable"))
    )

    if cmp["old_growth_input"] is not None and new_input.get("growth_rate") is not None:
        cmp["growth_diff_input"] = abs(
            cmp["old_growth_input"] - new_input["growth_rate"]
        )
        # The R pipeline rounds growth to 4 decimal places in the report; the
        # JSON itself stores rounded values too. Match within either 1e-6
        # OR within rounding of the recorded precision.
        cmp["input_growth_match"] = (
            cmp["growth_diff_input"] <= GROWTH_TOL
            or cmp["growth_diff_input"] <= 5e-5  # half ULP at 4 decimals
        )
    if cmp["old_growth_processed"] is not None and new_proc.get("growth_rate") is not None:
        cmp["growth_diff_processed"] = abs(
            cmp["old_growth_processed"] - new_proc["growth_rate"]
        )
        cmp["processed_growth_match"] = (
            cmp["growth_diff_processed"] <= GROWTH_TOL
            or cmp["growth_diff_processed"] <= 5e-5
        )

    old_in_size = old.get("input_model_size") or {}
    old_proc_size = old.get("processed_model_size") or {}
    cmp["old_input_metab"] = old_in_size.get("metabolites")
    cmp["old_input_rxn"] = old_in_size.get("reactions")
    cmp["old_proc_metab"] = old_proc_size.get("metabolites")
    cmp["old_proc_rxn"] = old_proc_size.get("reactions")

    if cmp["old_input_metab"] is not None and new_input.get("metabolites") is not None:
        cmp["input_size_match"] = (
            cmp["old_input_metab"] == new_input["metabolites"]
            and cmp["old_input_rxn"] == new_input["reactions"]
        )
    if cmp["old_proc_metab"] is not None and new_proc.get("metabolites") is not None:
        cmp["processed_size_match"] = (
            cmp["old_proc_metab"] == new_proc["metabolites"]
            and cmp["old_proc_rxn"] == new_proc["reactions"]
        )
    return cmp


def _run() -> dict[str, Any]:
    species_dirs = sorted([p for p in SPECIES_DIR.iterdir() if p.is_dir()])
    per_species: list[dict[str, Any]] = []

    for sp in species_dirs:
        print(f"[{sp.name}]", flush=True)
        inp, proc = _find_xml_pair(sp)
        rec: dict[str, Any] = {
            "species": sp.name,
            "input_file": inp.name if inp else None,
            "processed_file": proc.name if proc else None,
            "input": _try_read(inp) if inp else {"error": "no input file located"},
            "processed": _try_read(proc) if proc else {"error": "no processed file located"},
        }
        rec["old_validation"] = _load_old_validation(sp)
        rec["compare"] = _compare(rec["old_validation"], rec["input"], rec["processed"])
        # Cheap derived fields
        in_g = rec["input"].get("growth_rate")
        pr_g = rec["processed"].get("growth_rate")
        rec["growth_input"] = in_g
        rec["growth_processed"] = pr_g
        rec["growth_diff_input_processed"] = (
            abs(in_g - pr_g) if (in_g is not None and pr_g is not None) else None
        )
        rec["both_readable"] = bool(rec["input"].get("readable")) and bool(
            rec["processed"].get("readable")
        )
        rec["readability_change"] = _classify_readability(rec)
        per_species.append(rec)

    summary = _summarize(per_species)
    return {"summary": summary, "species": per_species}


def _classify_readability(rec: dict[str, Any]) -> str:
    ir = bool(rec["input"].get("readable"))
    pr = bool(rec["processed"].get("readable"))
    if ir and pr:
        return "both_readable"
    if ir and not pr:
        return "broken_by_pipeline"
    if not ir and pr:
        return "fixed_by_pipeline"
    return "both_unreadable"


def _summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(records)
    has_old = sum(1 for r in records if r["compare"]["old_present"])
    in_readable = sum(1 for r in records if r["input"].get("readable"))
    proc_readable = sum(1 for r in records if r["processed"].get("readable"))
    both_readable = sum(1 for r in records if r["both_readable"])
    broken = sum(1 for r in records if r["readability_change"] == "broken_by_pipeline")
    fixed = sum(1 for r in records if r["readability_change"] == "fixed_by_pipeline")
    both_unreadable = sum(
        1 for r in records if r["readability_change"] == "both_unreadable"
    )

    # Compare growth: only count species where the OLD recorded a numeric
    # growth value AND we re-computed a numeric growth value.
    in_growth_compared = 0
    in_growth_match = 0
    pr_growth_compared = 0
    pr_growth_match = 0
    mismatches: list[dict[str, Any]] = []
    readability_mismatches: list[dict[str, Any]] = []

    for r in records:
        cmp = r["compare"]
        if cmp.get("input_growth_match") is not None:
            in_growth_compared += 1
            if cmp["input_growth_match"]:
                in_growth_match += 1
            else:
                mismatches.append(
                    {
                        "species": r["species"],
                        "side": "input",
                        "old": cmp["old_growth_input"],
                        "new": r["input"].get("growth_rate"),
                        "diff": cmp["growth_diff_input"],
                    }
                )
        if cmp.get("processed_growth_match") is not None:
            pr_growth_compared += 1
            if cmp["processed_growth_match"]:
                pr_growth_match += 1
            else:
                mismatches.append(
                    {
                        "species": r["species"],
                        "side": "processed",
                        "old": cmp["old_growth_processed"],
                        "new": r["processed"].get("growth_rate"),
                        "diff": cmp["growth_diff_processed"],
                    }
                )

        if cmp.get("input_readable_match") is False:
            readability_mismatches.append(
                {
                    "species": r["species"],
                    "side": "input",
                    "old": cmp["old_input_readable"],
                    "new": r["input"].get("readable"),
                }
            )
        if cmp.get("processed_readable_match") is False:
            readability_mismatches.append(
                {
                    "species": r["species"],
                    "side": "processed",
                    "old": cmp["old_processed_readable"],
                    "new": r["processed"].get("readable"),
                }
            )

    growth_diverges_proc_vs_input = [
        {
            "species": r["species"],
            "input": r["growth_input"],
            "processed": r["growth_processed"],
            "diff": r["growth_diff_input_processed"],
        }
        for r in records
        if r["growth_diff_input_processed"] is not None
        and r["growth_diff_input_processed"] > 1e-6
    ]

    return {
        "total_species": total,
        "species_with_old_validation": has_old,
        "input_readable": in_readable,
        "processed_readable": proc_readable,
        "both_readable": both_readable,
        "broken_by_pipeline": broken,
        "fixed_by_pipeline": fixed,
        "both_unreadable": both_unreadable,
        "input_growth_compared": in_growth_compared,
        "input_growth_match": in_growth_match,
        "processed_growth_compared": pr_growth_compared,
        "processed_growth_match": pr_growth_match,
        "growth_mismatches": mismatches,
        "readability_mismatches": readability_mismatches,
        "growth_diverges_input_vs_processed": growth_diverges_proc_vs_input,
        "tolerance_used": GROWTH_TOL,
    }


def main() -> int:
    if not SPECIES_DIR.exists():
        print(f"missing species dir: {SPECIES_DIR}", file=sys.stderr)
        return 2
    try:
        report = _run()
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        return 1
    OUT_PATH.write_text(json.dumps(report, indent=2, default=str))
    s = report["summary"]
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    for k in (
        "total_species",
        "species_with_old_validation",
        "input_readable",
        "processed_readable",
        "both_readable",
        "broken_by_pipeline",
        "fixed_by_pipeline",
        "both_unreadable",
        "input_growth_compared",
        "input_growth_match",
        "processed_growth_compared",
        "processed_growth_match",
    ):
        print(f"  {k}: {s[k]}")
    print(f"  growth_mismatches: {len(s['growth_mismatches'])}")
    for m in s["growth_mismatches"]:
        print(f"    - {m['species']} [{m['side']}] old={m['old']} new={m['new']} diff={m['diff']}")
    print(f"  readability_mismatches: {len(s['readability_mismatches'])}")
    for m in s["readability_mismatches"]:
        print(f"    - {m['species']} [{m['side']}] old={m['old']} new={m['new']}")
    print(f"\nWrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
