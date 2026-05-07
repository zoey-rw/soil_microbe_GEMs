"""Re-generate tests/golden/smoke_cohort.json.

Run from the repo root:
    python tests/regenerate_golden.py

Adds a species:
    Edit COHORT below, then re-run.
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
import contextlib
import logging
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

COHORT = [
    "nitrobacter_winogradskyi_iFC579",
    "methanosarcina_barkeri_iMG746",
    "methanococcus_maripaludis_iMR557",
    "nitrosomonas_europaea_iGC535",
    "nitrosopumilus_maritimus_NmrFL413",
]


@contextlib.contextmanager
def silence_cobra():
    warnings.filterwarnings("ignore")
    logging.getLogger("cobra").setLevel(logging.CRITICAL)
    with open(os.devnull, "w") as devnull, contextlib.redirect_stderr(devnull):
        yield


def main() -> int:
    with silence_cobra():
        import cobra

    golden = {}
    for name in COHORT:
        sp = REPO_ROOT / "species" / name
        inputs = [p for p in sp.glob("*.xml") if "_processed" not in p.name]
        procs = list(sp.glob("*_processed.xml"))
        if not inputs or not procs:
            print(f"  SKIP {name}: missing input or processed", file=sys.stderr)
            continue
        inp = sorted(inputs, key=lambda p: len(p.name))[0]
        proc = procs[0]
        t0 = time.time()
        with silence_cobra():
            m_in = cobra.io.read_sbml_model(str(inp))
            m_pr = cobra.io.read_sbml_model(str(proc))
            g_in = float(m_in.slim_optimize())
            g_pr = float(m_pr.slim_optimize())
        dt = time.time() - t0
        golden[name] = {
            "input_file": inp.name,
            "processed_file": proc.name,
            "input": {
                "metabolites": len(m_in.metabolites),
                "reactions": len(m_in.reactions),
                "growth_rate": round(g_in, 6),
            },
            "processed": {
                "metabolites": len(m_pr.metabolites),
                "reactions": len(m_pr.reactions),
                "growth_rate": round(g_pr, 6),
            },
            "load_time_seconds": round(dt, 2),
        }
        print(f"  {name}: input g={g_in:.4f}, processed g={g_pr:.4f}, {dt:.1f}s")

    out_path = REPO_ROOT / "tests" / "golden" / "smoke_cohort.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(golden, indent=2, sort_keys=True) + "\n")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
