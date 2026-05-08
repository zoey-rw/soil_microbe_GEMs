"""Python helper for the cobra-via-reticulate SBML shim.

This module is the bridge that the R shim (pipeline/sbml_io_cobra.R) uses
in place of sybilSBML's readSBMLmod/writeSBML/findExchReact. It exposes
just what the pipeline actually touches:

  read_sbml(path)       -> dict of slots that R wraps in an S4 modelorg_cobra
  write_sbml(state)     -> applies slot deltas back to the cached cobra model
                           and writes SBML to disk
  exchanges(state)      -> dict mirroring sybilSBML's `exchReact` slots

The cobra model itself is held in an in-memory cache (keyed by an integer
"handle") so we never round-trip the full S-matrix through R/Python. R only
sees the slot vectors it needs to mutate.
"""
from __future__ import annotations

import contextlib
import logging
import os
import threading
import warnings
from itertools import count
from typing import Any, Dict, List, Optional, Tuple

# Suppress libsbml/cobra warnings at import time.
warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)

import cobra  # noqa: E402

_CACHE: Dict[int, "_ModelEntry"] = {}
_CACHE_LOCK = threading.Lock()
_HANDLE_COUNTER = count(1)


class _ModelEntry:
    __slots__ = (
        "model",
        "original_met_ids",
        "original_react_ids",
        "original_gprs",
        "original_path",
    )

    def __init__(self, model: cobra.Model, path: str) -> None:
        self.model = model
        self.original_met_ids = [m.id for m in model.metabolites]
        self.original_react_ids = [r.id for r in model.reactions]
        self.original_gprs = [r.gene_reaction_rule for r in model.reactions]
        self.original_path = path


@contextlib.contextmanager
def _silence():
    with open(os.devnull, "w") as devnull, contextlib.redirect_stderr(devnull):
        yield


def _compartment_index_map(model: cobra.Model) -> Tuple[List[str], Dict[str, int]]:
    """Return (compartment_codes, code -> 1-based index)."""
    codes = sorted(model.compartments.keys())
    return codes, {code: i + 1 for i, code in enumerate(codes)}


# Map cobra's namespace keys to the identifiers.org-style namespaces that
# the R-side detect_annotation_pattern() regex hunts for.
# Anything not in this map is passed through verbatim.
_NS_ALIASES = {
    "chebi": "chebi.compound",        # regex wants chebi.compound/ or CHEBI:
    "kegg":  "kegg.compound",
    "seed":  "seed.compound",
    "bigg":  "bigg.metabolite",
    "metanetx": "metanetx.chemical",
    "pubchem": "pubchem.compound",
    "hmdb":  "hmdb",
}


def _serialize_annotation(ann: Any) -> str:
    """Render a cobra metabolite.annotation dict to a string the R-side
    pattern detector recognises.

    The pipeline's detect_annotation_pattern() looks for substrings like
    'metanetx.chemical/', 'bigg.metabolite/', 'CHEBI:', 'chebi.compound/',
    etc. cobra normalises annotations into a dict keyed by namespace. We
    render each (k, v) pair in the form `<namespace>/<id>` joined by
    spaces. The CHEBI value is also emitted with its 'CHEBI:' prefix
    intact so the pipeline's `CHEBI:` alternation matches.
    """
    if not ann:
        return ""
    if isinstance(ann, str):
        return ann
    if not isinstance(ann, dict):
        return str(ann)

    parts: List[str] = []
    for namespace, value in ann.items():
        ns = str(namespace).lower()
        canonical_ns = _NS_ALIASES.get(ns, ns)
        values = value if isinstance(value, list) else [value]
        for v in values:
            sval = str(v)
            parts.append(f"{canonical_ns}/{sval}")
            # Emit CHEBI: prefix form too so the regex `CHEBI:` matches.
            if canonical_ns == "chebi.compound" and not sval.upper().startswith("CHEBI:"):
                parts.append(f"CHEBI:{sval}")
    return " ".join(parts)


def read_sbml(path: str) -> Dict[str, Any]:
    """Read an SBML file via cobra, register it in the cache, and return a
    dict the R shim can wrap in an S4 modelorg_cobra.

    The dict has the slots:
      handle               int       — cache key
      met_id               list[str] — metabolite ids (sybilSBML @met_id)
      met_name             list[str] — names (@met_name)
      met_comp             list[int] — 1-based compartment indices (@met_comp)
      met_annotation       list[str] — flattened annotations (@met_attr$annotation)
      mod_compart          list[str] — compartment codes by index (@mod_compart)
      react_id             list[str] — reaction ids (@react_id)
      react_rev            list[bool] — reversibility (@react_rev)
      gpr                  list[str] — GPR strings (@gpr)
    """
    with _silence():
        model = cobra.io.read_sbml_model(path)

    handle = next(_HANDLE_COUNTER)
    with _CACHE_LOCK:
        _CACHE[handle] = _ModelEntry(model, path)

    codes, code_to_idx = _compartment_index_map(model)
    met_comp = [code_to_idx.get(m.compartment, 0) for m in model.metabolites]

    return {
        "handle": handle,
        "met_id": [m.id for m in model.metabolites],
        "met_name": [m.name or "" for m in model.metabolites],
        "met_comp": met_comp,
        "met_annotation": [_serialize_annotation(m.annotation) for m in model.metabolites],
        "mod_compart": codes,
        "react_id": [r.id for r in model.reactions],
        "react_name": [r.name or r.id for r in model.reactions],
        "react_rev": [r.lower_bound < 0 for r in model.reactions],
        "obj_coef": [float(r.objective_coefficient) for r in model.reactions],
        "lowbnd": [float(r.lower_bound) for r in model.reactions],
        "uppbnd": [float(r.upper_bound) for r in model.reactions],
        "gpr": [r.gene_reaction_rule for r in model.reactions],
    }


def _entry(handle: int) -> _ModelEntry:
    with _CACHE_LOCK:
        try:
            return _CACHE[int(handle)]
        except KeyError:
            raise KeyError(
                f"sbml_io: model handle {handle} is not in the cache. "
                "Call read_sbml() first or check that the R shim hasn't dropped the reference."
            )


def _apply_id_renames(
    objects, original_ids: List[str], new_ids: List[str], object_kind: str
) -> int:
    """Apply id renames to a cobra DictList (metabolites or reactions) by
    iterating with a temporary token. Returns the number of renames applied.
    """
    if len(original_ids) != len(new_ids):
        raise ValueError(
            f"{object_kind} count changed: was {len(original_ids)}, now {len(new_ids)}"
        )
    if len(objects) != len(original_ids):
        raise ValueError(
            f"{object_kind} length mismatch: cobra has {len(objects)} but "
            f"original_ids has {len(original_ids)}"
        )

    renames = 0
    # Two-pass to avoid id collisions when swaps cross.
    tmp_token = "__sbml_io_pending__"
    for i, (old, new) in enumerate(zip(original_ids, new_ids)):
        if old != new:
            objects[i].id = f"{tmp_token}{i}"
    for i, (old, new) in enumerate(zip(original_ids, new_ids)):
        if old != new:
            objects[i].id = new
            renames += 1
    return renames


def write_sbml(state: Dict[str, Any], filename: str, level: int = 3) -> Dict[str, Any]:
    """Apply the slot vectors in ``state`` back to the cached cobra model
    and write SBML to ``filename``.

    Slots applied:
      met_id            -> rename metabolites
      react_id          -> rename reactions
      gpr               -> reaction.gene_reaction_rule
      obj_coef          -> reaction.objective_coefficient
      lowbnd / uppbnd   -> reaction.lower_bound / upper_bound

    Returns a small dict with rename counts so the R caller can log them.
    """
    handle = int(state["handle"])
    entry = _entry(handle)
    model = entry.model

    new_met_ids = list(state.get("met_id", entry.original_met_ids))
    new_react_ids = list(state.get("react_id", entry.original_react_ids))
    new_gprs = list(state.get("gpr", entry.original_gprs))

    met_renames = _apply_id_renames(
        model.metabolites, entry.original_met_ids, new_met_ids, "metabolite"
    )
    react_renames = _apply_id_renames(
        model.reactions, entry.original_react_ids, new_react_ids, "reaction"
    )

    gpr_updates = 0
    for i, (old, new) in enumerate(zip(entry.original_gprs, new_gprs)):
        if (old or "") != (new or ""):
            model.reactions[i].gene_reaction_rule = new or ""
            gpr_updates += 1

    # Apply objective + bound updates if provided. These are only present
    # when the R caller is using the modelorg_cobra slots that were added
    # in the issue #6 work — older callers omit them.
    obj_updates = 0
    bound_updates = 0
    new_obj = state.get("obj_coef")
    if new_obj is not None and len(new_obj) == len(model.reactions):
        for i, coef in enumerate(new_obj):
            current = float(model.reactions[i].objective_coefficient)
            if abs(current - float(coef)) > 1e-12:
                model.reactions[i].objective_coefficient = float(coef)
                obj_updates += 1

    new_lb = state.get("lowbnd")
    new_ub = state.get("uppbnd")
    if (new_lb is not None and new_ub is not None
            and len(new_lb) == len(model.reactions)
            and len(new_ub) == len(model.reactions)):
        for i, (lb, ub) in enumerate(zip(new_lb, new_ub)):
            r = model.reactions[i]
            cur_lb, cur_ub = float(r.lower_bound), float(r.upper_bound)
            new_lb_f, new_ub_f = float(lb), float(ub)
            if abs(cur_lb - new_lb_f) > 1e-12 or abs(cur_ub - new_ub_f) > 1e-12:
                r.bounds = (new_lb_f, new_ub_f)
                bound_updates += 1

    # Refresh the entry's "originals" so a second write_sbml call applies
    # only the new deltas.
    entry.original_met_ids = list(new_met_ids)
    entry.original_react_ids = list(new_react_ids)
    entry.original_gprs = list(new_gprs)

    out_dir = os.path.dirname(os.path.abspath(filename))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with _silence():
        cobra.io.write_sbml_model(model, filename)

    return {
        "filename": filename,
        "metabolite_renames": met_renames,
        "reaction_renames": react_renames,
        "gpr_updates": gpr_updates,
        "objective_updates": obj_updates,
        "bound_updates": bound_updates,
    }


def exchanges(state: Dict[str, Any]) -> Dict[str, Any]:
    """Return a dict mirroring sybilSBML's exchReact slot layout.

    Keys: react_id, met_id, uptake (boolean), lower_bound, upper_bound.
    """
    handle = int(state["handle"])
    entry = _entry(handle)
    model = entry.model

    react_ids: List[str] = []
    met_ids: List[str] = []
    uptakes: List[bool] = []
    lbs: List[float] = []
    ubs: List[float] = []

    for r in model.exchanges:
        react_ids.append(r.id)
        # cobra exchanges have exactly one metabolite; pick first key
        met_ids.append(next(iter(r.metabolites)).id if r.metabolites else "")
        # An exchange is an "uptake" candidate if it can carry negative flux
        uptakes.append(r.lower_bound < 0)
        lbs.append(float(r.lower_bound))
        ubs.append(float(r.upper_bound))

    return {
        "react_id": react_ids,
        "met_id": met_ids,
        "uptake": uptakes,
        "lower_bound": lbs,
        "upper_bound": ubs,
    }


def release(handle: int) -> bool:
    """Drop a model from the cache (called by the R shim's finalizer)."""
    with _CACHE_LOCK:
        return _CACHE.pop(int(handle), None) is not None


def cache_size() -> int:
    """For diagnostics/tests."""
    with _CACHE_LOCK:
        return len(_CACHE)
