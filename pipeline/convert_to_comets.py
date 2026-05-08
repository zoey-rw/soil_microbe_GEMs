#!/usr/bin/env python3
"""Convert a processed SBML model into COMETS-format files.

The pipeline produces standardised SBML at species/<name>/<id>_processed.xml.
COMETS dynamic-FBA simulations want a different on-disk format:
a `.cmd` model file (S-matrix + bounds + objective + exchange list) and,
optionally, a `.layout` file describing the spatial / multi-species setup.

This script automates the per-model side: read the processed SBML via
cobra, build a cometspy.model, and write the .cmd file. Layout assembly
(combining multiple .cmd files into a multi-species simulation) is a
separate concern; see scripts/build_comets_layout.py once that lands.

Usage:
    python pipeline/convert_to_comets.py <input.xml> [--output <dir>]
                                         [--id <model_id>]
                                         [--no-warn-no-biomass]

Examples:
    # Single species
    python pipeline/convert_to_comets.py \
        species/nitrobacter_winogradskyi_iFC579/iFC579_processed.xml \
        --output build/comets/

    # Batch — bash loop is fine; one model per call keeps memory bounded
    for f in species/*/*_processed.xml; do
        python pipeline/convert_to_comets.py "$f" --output build/comets/
    done
"""
from __future__ import annotations

import argparse
import contextlib
import logging
import os
import re
import shutil
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Optional

# Suppress import-time noise.
warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)


@contextlib.contextmanager
def _silence():
    with open(os.devnull, "w") as devnull, contextlib.redirect_stderr(devnull):
        yield


_BIOMASS_RE = re.compile(r"biomass", re.IGNORECASE)


def _detect_biomass(model) -> Optional[str]:
    """Return the id of the most plausible biomass reaction, or None."""
    candidates = []
    for r in model.reactions:
        if _BIOMASS_RE.search(r.id) or _BIOMASS_RE.search(r.name or ""):
            candidates.append(r)
    if not candidates:
        return None
    # Prefer the longest id (more specific names like BIOMASS_Gm_GS15_core)
    return max(candidates, key=lambda r: len(r.id)).id


def _ensure_objective(cobra_model, biomass_id: Optional[str]) -> Optional[str]:
    """If the model's objective is empty or pointing at an exchange, redirect
    it to the detected biomass reaction. Returns the chosen objective id (or
    None if we couldn't find a sensible target)."""
    obj_rxns = [r for r in cobra_model.reactions if r.objective_coefficient != 0]
    if obj_rxns and not all(_is_exchangeish(r) for r in obj_rxns):
        return obj_rxns[0].id  # already a non-exchange objective; trust it

    if biomass_id is None:
        return None

    cobra_model.objective = biomass_id
    return biomass_id


def _is_exchangeish(rxn) -> bool:
    """Is this reaction an exchange/sink/demand? cobra exposes .boundary."""
    try:
        return bool(rxn.boundary)
    except Exception:
        return rxn.id.startswith("EX_") or rxn.id.startswith("DM_") or rxn.id.startswith("SK_")


def convert(
    input_path: Path,
    output_dir: Path,
    model_id: Optional[str] = None,
    warn_no_biomass: bool = True,
    optimizer: Optional[str] = None,
) -> Path:
    """Read input_path with cobra, build a cometspy model, write .cmd to
    output_dir. Returns the path of the written .cmd file."""
    import cobra
    import cometspy

    if not input_path.exists():
        raise FileNotFoundError(f"input SBML does not exist: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    if model_id is None:
        # Strip "_processed" or any other "_*" suffix that follows the id
        stem = input_path.stem.replace("_processed", "")
        # If the directory name contains a "_<modelid>" suffix, prefer that
        parent_id = input_path.parent.name.split("_")[-1]
        model_id = parent_id if parent_id else stem

    with _silence():
        cobra_model = cobra.io.read_sbml_model(str(input_path))
    cobra_model.id = model_id

    biomass_id = _detect_biomass(cobra_model)
    objective_id = _ensure_objective(cobra_model, biomass_id)

    if objective_id is None and warn_no_biomass:
        print(
            f"WARNING: {input_path.name} has no detectable biomass reaction; "
            "COMETS model will lack a meaningful objective. Fix upstream "
            "(see issue #6) or pass --no-warn-no-biomass.",
            file=sys.stderr,
        )
    elif (
        objective_id is not None
        and biomass_id is not None
        and objective_id != biomass_id
        and warn_no_biomass
    ):
        print(
            f"NOTE: {input_path.name} kept its existing objective "
            f"({objective_id!r}) rather than overriding to detected "
            f"biomass {biomass_id!r}.",
            file=sys.stderr,
        )

    # Build cometspy model. cometspy appends the cobra_model.id to whatever
    # path we pass (so passing "foo.cmd" produces "foo.cmdiFC579.cmd"). Work
    # around by writing into a temp dir and then renaming.
    cm = cometspy.model(cobra_model)
    if optimizer is not None:
        cm.change_optimizer(optimizer)
    with tempfile.TemporaryDirectory() as tmp:
        cm.write_comets_model(os.path.join(tmp, ""))
        produced = [p for p in os.listdir(tmp) if p.endswith(".cmd")]
        if not produced:
            raise RuntimeError(
                f"cometspy.write_comets_model wrote no .cmd in {tmp}; "
                f"directory contents: {os.listdir(tmp)}"
            )
        # Rename to the canonical "<model_id>.cmd"
        out_path = output_dir / f"{model_id}.cmd"
        shutil.move(os.path.join(tmp, produced[0]), out_path)

    return out_path


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("input", type=Path, help="processed SBML file")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("build/comets"),
        help="output directory (default: build/comets/)",
    )
    p.add_argument(
        "--id",
        dest="model_id",
        default=None,
        help="override the model_id used for the output filename",
    )
    p.add_argument(
        "--no-warn-no-biomass",
        action="store_true",
        help="suppress warnings about missing biomass / objective",
    )
    p.add_argument(
        "--optimizer",
        default=None,
        choices=[None, "GUROBI", "GLPK"],
        help="LP solver for the COMETS run (default: cometspy's default, GUROBI). "
             "Use GLPK for license-free runs.",
    )
    args = p.parse_args(argv)

    out = convert(
        args.input,
        args.output,
        model_id=args.model_id,
        warn_no_biomass=not args.no_warn_no_biomass,
        optimizer=args.optimizer,
    )
    print(f"wrote {out} ({out.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
