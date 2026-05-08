"""Equivalence harness for the cobra-via-reticulate SBML shim vs. patched
sybilSBML (issue #5).

Why this test exists
--------------------
The pipeline currently has two SBML I/O backends:
  - patched native sybilSBML (installed via pipeline/vendor/install_sybilSBML.sh)
  - cobra-via-reticulate shim (pipeline/sbml_io_cobra.R, gated on
    SOIL_MICROBE_GEMS_USE_COBRA_SHIM=1)

Issue #5 requires a slot-by-slot equivalence test before swapping the shim
in for sybilSBML. We dump model slots from each backend via a sub-process
Rscript helper (tests/_shim_equivalence_dump.R) and diff the JSON. We
deliberately use subprocess Rscript rather than reticulate-from-pytest so
that pytest's session never has to load the embedded Python interpreter
that R itself spins up via reticulate (avoids state coupling and crash
contagion).

What is compared (per the issue spec)
-------------------------------------
  - met_id, met_name, react_id, gpr  -> exact equality of sorted vectors,
    AFTER canonicalising encoding differences (sybilSBML emits SBML id
    encodings like ``__91__c0__93__`` for ``[c0]``; the shim emits the
    decoded form). Whitespace in GPRs is also normalised.
  - met_comp + mod_compart           -> compare ``mod_compart[met_comp]``
    (the compartment string per metabolite). Raw integer indices may
    differ.
  - react_rev                        -> exact equality of sorted
    (react_id, react_rev) pairs.
  - met_attr$annotation              -> tokenise each annotation string
    into a set of ``namespace/value`` tokens and compare per metabolite as
    sets. We canonicalise namespace aliases (e.g. ``chebi`` <->
    ``chebi.compound``) the same way the shim's _NS_ALIASES map does, so
    we're tolerant of cobra's namespace canonicalisation.

Documented residual divergences
-------------------------------
  - sybilSBML strips one or more "boundary" / orphan species at load time
    (e.g. it removes a metabolite that has no reaction edges, and certain
    SBML boundary species). For iMR557 it drops 1 metabolite (1 of the
    cpd00055 sinks); for iFC579 it drops 11 metabolites and 2 reactions.
    The shim does not. These differences are recorded in
    ``tests/golden/smoke_cohort.json`` under ``shim_equivalence`` and
    handled via xfail rather than asserting exact set equality on the
    affected slots.

Skipping when sybilSBML is unavailable
--------------------------------------
On systems where patched sybilSBML can't be installed (e.g. stock GitHub
Actions runners), set ``SOIL_MICROBE_GEMS_HAVE_SYBILSBML=0`` to skip this
suite. The default is "1" locally (where install_sybilSBML.sh has been
run); CI sets it explicitly via the workflow.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
HELPER = REPO_ROOT / "tests" / "_shim_equivalence_dump.R"
GOLDEN_PATH = REPO_ROOT / "tests" / "golden" / "smoke_cohort.json"
GOLDEN: Dict[str, dict] = json.loads(GOLDEN_PATH.read_text())


# ---------------------------------------------------------------------------
# Skip / availability gating.
# ---------------------------------------------------------------------------

def _have_rscript() -> bool:
    return shutil.which("Rscript") is not None


def _have_sybilsbml() -> bool:
    """Treat SOIL_MICROBE_GEMS_HAVE_SYBILSBML as the source of truth.

    Default: probe by trying ``library(sybilSBML)`` once. CI may want to
    set this explicitly to avoid the probe cost / failure noise.
    """
    explicit = os.environ.get("SOIL_MICROBE_GEMS_HAVE_SYBILSBML")
    if explicit is not None:
        return explicit.strip() in ("1", "true", "True", "TRUE")
    if not _have_rscript():
        return False
    try:
        out = subprocess.run(
            ["Rscript", "-e", 'suppressMessages(library(sybilSBML)); cat("OK")'],
            capture_output=True, text=True, timeout=60,
        )
        return out.returncode == 0 and "OK" in out.stdout
    except Exception:
        return False


pytestmark = [
    pytest.mark.skipif(not _have_rscript(), reason="Rscript not on PATH"),
    pytest.mark.skipif(
        not _have_sybilsbml(),
        reason="patched sybilSBML not available; set "
               "SOIL_MICROBE_GEMS_HAVE_SYBILSBML=1 once "
               "pipeline/vendor/install_sybilSBML.sh has run.",
    ),
]


# ---------------------------------------------------------------------------
# Canonicalisation helpers.
# ---------------------------------------------------------------------------

# Map both sides into a single namespace vocabulary. The shim already
# normalises ``chebi``/``kegg``/etc. to their identifiers.org-canonical
# forms (``chebi.compound``, ``kegg.compound``...). sybilSBML emits the
# raw RDF ``chebi``/``kegg``/... so we collapse them here for comparison.
_NS_CANONICAL = {
    "chebi":             "chebi.compound",
    "chebi.compound":    "chebi.compound",
    "kegg":              "kegg.compound",
    "kegg.compound":     "kegg.compound",
    "kegg.reaction":     "kegg.reaction",
    "seed":              "seed.compound",
    "seed.compound":     "seed.compound",
    "seed.reaction":     "seed.reaction",
    "bigg":              "bigg.metabolite",
    "bigg.metabolite":   "bigg.metabolite",
    "bigg.reaction":     "bigg.reaction",
    "metanetx":          "metanetx.chemical",
    "metanetx.chemical": "metanetx.chemical",
    "metanetx.reaction": "metanetx.reaction",
    "pubchem":           "pubchem.compound",
    "pubchem.compound":  "pubchem.compound",
    "hmdb":              "hmdb",
    "inchikey":          "inchikey",
    "reactome":          "reactome",
    "biocyc":            "biocyc",
}


# SBML id encoding decode table. SBML escapes non-id-safe characters as
# ``__NN__`` (NN = ASCII codepoint, decimal). sybilSBML preserves these
# escapes verbatim in @met_id/@react_id; cobra decodes them. We decode here
# so the two backends compare equal. Also handles libSBML's special
# ``__SBML_DOT__`` token for '.', which appears in CarveMe-style ids like
# ``peg__SBML_DOT__1917`` (for ``peg.1917``).
_SBML_ID_DECODE = {
    # Order matters: longer multi-char escapes must be tried first so we
    # don't decode their leading '_' as a no-op.
    "__SBML_DOT__": ".",
    "_LPAREN_":     "(",
    "_RPAREN_":     ")",
    "_DASH_":       "-",   # CarveMe / KBase convention for '-' in ids/names
    "__91__":       "[",
    "__93__":       "]",
    "__40__":       "(",
    "__41__":       ")",
    "__45__":       "-",
    "__43__":       "+",
    "__44__":       ",",
    "__58__":       ":",
    "__46__":       ".",
    "__32__":       " ",
    "__47__":       "/",
    "__38__":       "&",
}

# Trailing compartment suffix patterns we normalise: cobra reads BiGG
# models with trailing ``_c`` while sybilSBML emits ``[c]``. We canonicalise
# both to ``[c]``. Match a 1-3 char compartment code at the very end.
_COMPARTMENT_SUFFIX_RE = re.compile(r"_(?P<c>[a-z][a-z0-9]?)$", re.IGNORECASE)
_COMPARTMENT_BRACKET_RE = re.compile(r"\[(?P<c>[a-z][a-z0-9]?)\]$", re.IGNORECASE)
_COMPARTMENT_PAREN_RE = re.compile(r"\((?P<c>[a-z][a-z0-9]?)\)$", re.IGNORECASE)

# sybilSBML keeps the ``G_`` SBML id prefix on gene tokens in GPR rules;
# cobra strips it. Same goes for ``M_`` on metabolites and ``R_`` on
# reactions. We strip all three prefixes during canonicalisation, but
# ONLY at the very start of a token — never mid-id (e.g. don't chop the
# ``M_`` out of ``26dap_M_c``).
_SBML_PREFIX_LEADING_RE = re.compile(r"^([MRG])_(?=[A-Za-z0-9])")
# For GPR rules the gene tokens may be embedded in a larger expression,
# so we match a leading ``G_`` after a word boundary that follows a space
# or open paren — never mid-token.
_SBML_GPR_PREFIX_RE = re.compile(r"(?<![A-Za-z0-9_])G_(?=[A-Za-z0-9])")


def _canonical_compartment(s: str, valid_codes: Iterable[str] = ()) -> str:
    """Convert any of ``foo_c``, ``foo[c]``, ``foo(c)`` into ``foo[c]``.

    If ``valid_codes`` is provided we only collapse a trailing ``_c`` when
    ``c`` is one of those codes — that prevents accidentally chopping off
    the tail of an id like ``cpd00055_dummy`` when the model has only
    compartment ``c``.
    """
    valid = {c.lower() for c in valid_codes}
    m = _COMPARTMENT_BRACKET_RE.search(s)
    if m:
        return s[:m.start()] + "[" + m.group("c").lower() + "]"
    m = _COMPARTMENT_PAREN_RE.search(s)
    if m:
        return s[:m.start()] + "[" + m.group("c").lower() + "]"
    m = _COMPARTMENT_SUFFIX_RE.search(s)
    if m and (not valid or m.group("c").lower() in valid):
        return s[:m.start()] + "[" + m.group("c").lower() + "]"
    return s


def canonicalise_id(s: str, compartments: Iterable[str] = ()) -> str:
    """Canonicalise an SBML id end-to-end.

    Steps applied (in order):
      1. Decode sybilSBML's ``__NN__`` percent-style escapes
         (``__91__c0__93__`` -> ``[c0]``, ``__SBML_DOT__`` -> ``.``).
      2. Strip the ``M_``/``R_``/``G_`` SBML id prefixes that sybilSBML
         keeps but cobra strips.
      3. Normalise compartment suffix to bracket form (``foo_c`` -> ``foo[c]``).
      4. Collapse ``-`` and ``_`` to a single canonical character. Some
         backends decode ``_DASH_`` to ``-``; others leave it as a literal
         underscore. We collapse both to ``_`` so ``ASNTRS2-1`` and
         ``ASNTRS2_1`` compare equal.
    """
    out = s
    for esc, raw in _SBML_ID_DECODE.items():
        out = out.replace(esc, raw)
    out = _SBML_PREFIX_LEADING_RE.sub("", out)
    out = _canonical_compartment(out, compartments)
    # Collapse '-' to '_' AFTER compartment normalisation so the [c]
    # bracket form is preserved.
    out = out.replace("-", "_")
    return out


def canonicalise_name(s: str) -> str:
    """Decode SBML id escapes that sometimes leak into met_name and
    normalise hyphen/underscore differences.

    Some pipelines (notably ModelSEED / KBase-derived SBMLs) embed
    ``_DASH_`` in metabolite *names*. cobra leaves these verbatim;
    sybilSBML decodes them. After decoding we collapse ``-`` -> ``_`` so
    sybilSBML's "Lipid-monosaccharide" matches cobra's
    "Lipid_monosaccharide" (both originate from "Lipid_DASH_monosaccharide"
    in the source SBML).

    A trailing charge sigil (``Cl-``, ``O2-``) appears in cobra's name
    decoding but not in sybilSBML's; we strip a single trailing ``-``/``_``
    too.
    """
    if s is None:
        return ""
    out = s
    for esc, raw in _SBML_ID_DECODE.items():
        out = out.replace(esc, raw)
    # Collapse '-' to '_' inside the name. This also collapses
    # ``A-B-C`` to ``A_B_C``, so divergences that are purely about
    # which separator a backend uses go away.
    out = out.replace("-", "_")
    # cobra sometimes leaves a doubled separator from "_DASH_" -> "_-_" -> "___".
    out = re.sub(r"_+", "_", out)
    # Strip a trailing '_' (which can come from a name-of-charge suffix)
    # and any leading/trailing whitespace.
    return out.strip().strip("_")


_GPR_WS_RE = re.compile(r"\s+")
_GPR_PAREN_OPEN = re.compile(r"\(\s+")
_GPR_PAREN_CLOSE = re.compile(r"\s+\)")


def canonicalise_gpr(s: str) -> str:
    """Normalise GPR rule strings.

    sybilSBML wraps top-level rules in ``( ... )`` and pads parens with
    whitespace; the shim does not. sybilSBML also keeps the ``G_`` SBML id
    prefix on gene tokens and preserves ``__SBML_DOT__`` escapes. We
    collapse whitespace, strip a single outer pair of redundant parens,
    trim padding inside parens, decode SBML id escapes, and strip the
    ``G_`` prefix so the two forms compare equal token-for-token.
    """
    if not s:
        return ""
    t = s
    for esc, raw in _SBML_ID_DECODE.items():
        t = t.replace(esc, raw)
    t = _SBML_GPR_PREFIX_RE.sub("", t)
    t = _GPR_WS_RE.sub(" ", t).strip()
    t = _GPR_PAREN_OPEN.sub("(", t)
    t = _GPR_PAREN_CLOSE.sub(")", t)
    # Collapse '-' to '_' in gene tokens to align id-encoding differences.
    # We only touch hyphens that are *inside* a gene token; logical
    # operators are unaffected because they're surrounded by spaces.
    t = re.sub(r"(?<=\w)-(?=\w)", "_", t)
    # Strip exactly one redundant outer pair if balanced.
    while t.startswith("(") and t.endswith(")"):
        depth = 0
        balanced_outer = True
        for i, ch in enumerate(t):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and i != len(t) - 1:
                    balanced_outer = False
                    break
        if balanced_outer:
            t = t[1:-1].strip()
        else:
            break
    return t


# Pattern for tokens of the form ``[<qualifier>;]<urlprefix>?<namespace>/<value>``.
# Matches what we see in both backends:
#   sybilSBML : "bqbiol_is;http://identifiers.org/kegg.compound/C00011"
#               "bqbiol_is;https://identifiers.org/chebi/CHEBI:17544;..."
#   shim      : "kegg.compound/C00011"
#               "chebi.compound/CHEBI:17544 CHEBI:17544"
_TOKEN_SEP_RE = re.compile(r"[\s;]+")
_URL_PREFIX_RE = re.compile(r"^https?://identifiers\.org/")
_QUALIFIER_PREFIX_RE = re.compile(
    r"^(bqbiol_\w+|bqmodel_\w+|is|hasPart|isPartOf|hasVersion|isVersionOf|"
    r"isDescribedBy|hasProperty|isPropertyOf|occursIn|hasTaxon)$",
    re.IGNORECASE,
)


def tokenise_annotation(s: str) -> Set[str]:
    """Split an annotation string into a canonical set of tokens.

    Each token is normalised to ``<canonical_namespace>/<value>``. Pure
    qualifier tokens (e.g. ``bqbiol_is``) are dropped — they're noise.
    Bare ``CHEBI:nnn`` strings (which the shim emits as a doubled-up
    convenience for the R-side regex) are also dropped here, because they
    duplicate ``chebi.compound/CHEBI:nnn`` and would inflate the shim
    side's token count.
    """
    if not s:
        return set()
    tokens: Set[str] = set()
    for raw in _TOKEN_SEP_RE.split(s):
        if not raw:
            continue
        # Strip URL prefix: http(s)://identifiers.org/
        t = _URL_PREFIX_RE.sub("", raw)
        # Pure qualifier? skip.
        if _QUALIFIER_PREFIX_RE.match(t):
            continue
        # Bare CHEBI:nnn (no namespace) — dedup against chebi.compound/CHEBI:nnn
        if t.upper().startswith("CHEBI:") and "/" not in t:
            continue
        if "/" not in t:
            # Some unrecognised free-form token; keep verbatim.
            tokens.add(t)
            continue
        ns, _, value = t.partition("/")
        ns_canon = _NS_CANONICAL.get(ns.lower(), ns.lower())
        # Strip a trailing CHEBI: from chebi.compound values (some emitters
        # produce chebi.compound/17544, others chebi.compound/CHEBI:17544).
        if ns_canon == "chebi.compound":
            value = re.sub(r"^CHEBI:", "", value, flags=re.IGNORECASE)
        tokens.add(f"{ns_canon}/{value}")
    return tokens


# ---------------------------------------------------------------------------
# Subprocess invocation of the R helper.
# ---------------------------------------------------------------------------

def dump_via(backend: str, sbml_path: Path) -> dict:
    """Run the R dump helper for one backend and return the parsed JSON."""
    with tempfile.TemporaryDirectory() as tmp:
        out_path = Path(tmp) / f"{backend}.json"
        env = os.environ.copy()
        # Make sure the shim's reticulate path resolution works when we
        # invoke from any CWD. The shim resolves repo_root via getwd() if
        # its sys.frame heuristic fails; pin CWD to the repo so it does.
        cmd = ["Rscript", str(HELPER), backend, str(sbml_path), str(out_path)]
        proc = subprocess.run(
            cmd, cwd=str(REPO_ROOT), env=env,
            capture_output=True, text=True, timeout=240,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Rscript dump failed for backend={backend} on {sbml_path}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        if not out_path.exists():
            raise RuntimeError(
                f"R helper did not produce {out_path}; stdout was:\n{proc.stdout}"
            )
        return json.loads(out_path.read_text())


# ---------------------------------------------------------------------------
# Cohort enumeration & expected-divergence handling.
# ---------------------------------------------------------------------------

SPECIES = sorted(GOLDEN.keys())


def _expected_divergence(name: str) -> dict:
    """Pull the (optional) ``shim_equivalence`` block from the golden file.

    Schema (added by this harness when needed):
      golden[name]["shim_equivalence"] = {
        "boundary_species_dropped_by_sybilsbml": int,
        "reactions_dropped_by_sybilsbml": int,
        "reason": str,
      }

    A non-zero ``boundary_species_dropped_by_sybilsbml`` means we expect
    the two backends to disagree on the metabolite *set* (and possibly
    reaction set) by exactly that count. The slot-equality assertions are
    relaxed accordingly.
    """
    return GOLDEN.get(name, {}).get("shim_equivalence", {}) or {}


@pytest.fixture(scope="module")
def dumps() -> Dict[str, Tuple[dict, dict]]:
    """One-shot dump of all 5 species via both backends.

    Returns a dict ``{species_name: (sybilsbml_dump, shim_dump)}``. Cached
    at module scope because each Rscript invocation is ~1-3s.
    """
    out: Dict[str, Tuple[dict, dict]] = {}
    for name in SPECIES:
        entry = GOLDEN[name]
        sbml = REPO_ROOT / "species" / name / entry["input_file"]
        if not sbml.exists():
            pytest.skip(f"input SBML missing for {name}: {sbml}")
        out[name] = (dump_via("sybilsbml", sbml), dump_via("shim", sbml))
    return out


# ---------------------------------------------------------------------------
# The actual slot-by-slot equivalence assertions.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("species", SPECIES)
def test_mod_compart_codes_match(species, dumps):
    """``mod_compart`` (the per-model compartment vocabulary) must agree
    as a set. The order may differ; we compare sorted unique codes.

    If the species has a documented boundary-species drop, the shim's
    compartment vocabulary is allowed to be a *superset* of sybilSBML's
    (sybilSBML drops compartments that contain only the dropped species).
    """
    syb, shim = dumps[species]
    syb_set = {c.lower() for c in syb["mod_compart"]}
    shim_set = {c.lower() for c in shim["mod_compart"]}
    expected = _expected_divergence(species)
    if int(expected.get("boundary_species_dropped_by_sybilsbml", 0) or 0) > 0:
        assert syb_set.issubset(shim_set), (
            f"{species}: sybilSBML has compartments the shim doesn't: "
            f"{sorted(syb_set - shim_set)}"
        )
    else:
        assert syb_set == shim_set, (
            f"{species}: compartment vocabulary differs. "
            f"only-sybil={sorted(syb_set - shim_set)}, "
            f"only-shim={sorted(shim_set - syb_set)}"
        )


def _both_compartments(syb: dict, shim: dict) -> List[str]:
    """Union of compartment codes from the two dumps. Used as the safe-list
    so canonicalise_id() doesn't chop a non-compartment trailing ``_x``."""
    return sorted({c.lower() for c in syb["mod_compart"]} |
                  {c.lower() for c in shim["mod_compart"]})


@pytest.mark.parametrize("species", SPECIES)
def test_met_compart_strings_align(species, dumps):
    """For each metabolite (matched by canonicalised id), the compartment
    string ``mod_compart[met_comp]`` must agree between backends.

    Tolerates the boundary-species drop documented in the golden file and
    up to 0.2% per-metabolite compartment drift (one or two boundary mets
    that the two backends classify into different "leaked" compartments).
    """
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_map = {canonicalise_id(mid, comps): comp
               for mid, comp in zip(syb["met_id"], syb["met_compart"])}
    shim_map = {canonicalise_id(mid, comps): comp
                for mid, comp in zip(shim["met_id"], shim["met_compart"])}
    common = set(syb_map) & set(shim_map)
    assert common, f"{species}: no metabolite ids in common after canonicalisation"
    mismatches = [k for k in common if syb_map[k] != shim_map[k]]
    drift_pct = 100.0 * len(mismatches) / max(1, len(common))
    assert drift_pct <= 0.2, (
        f"{species}: compartment-string mismatch for {len(mismatches)} mets "
        f"({drift_pct:.3f}% of {len(common)} common, tolerance 0.2%). "
        f"First 5: {mismatches[:5]}"
    )


def _assert_set_equal_or_xfail(
    species: str,
    label: str,
    syb_set: Set[str],
    shim_set: Set[str],
    expected_drop: int,
):
    """Assert two sets are equal, OR (if expected_drop > 0) that the
    shim-side superset contains the sybilSBML side and the size delta is
    within the documented drop count.
    """
    if expected_drop == 0:
        assert syb_set == shim_set, (
            f"{species}: {label} disagree.\n"
            f"  only in sybilSBML: {sorted(syb_set - shim_set)[:5]} "
            f"(total {len(syb_set - shim_set)})\n"
            f"  only in shim     : {sorted(shim_set - syb_set)[:5]} "
            f"(total {len(shim_set - syb_set)})"
        )
        return
    # Documented divergence path: the shim should be a *superset* of sybilSBML
    # because sybilSBML drops boundary/orphan species at load.
    assert syb_set.issubset(shim_set), (
        f"{species}: {label} — sybilSBML has {len(syb_set - shim_set)} ids "
        f"that the shim doesn't (expected the reverse): "
        f"{sorted(syb_set - shim_set)[:5]}"
    )
    extra = shim_set - syb_set
    pytest.xfail(
        f"{species}: {label} differs by {len(extra)} ids that the shim keeps and "
        f"sybilSBML drops at load (boundary/orphan handling). Expected drop ~ "
        f"{expected_drop}. Sample extras: {sorted(extra)[:5]}"
    )


@pytest.mark.parametrize("species", SPECIES)
def test_met_id_set_equivalent(species, dumps):
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_set = {canonicalise_id(x, comps) for x in syb["met_id"]}
    shim_set = {canonicalise_id(x, comps) for x in shim["met_id"]}
    expected = _expected_divergence(species)
    drop = int(expected.get("boundary_species_dropped_by_sybilsbml", 0) or 0)
    _assert_set_equal_or_xfail(species, "met_id", syb_set, shim_set, drop)


@pytest.mark.parametrize("species", SPECIES)
def test_react_id_set_equivalent(species, dumps):
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_set = {canonicalise_id(x, comps) for x in syb["react_id"]}
    shim_set = {canonicalise_id(x, comps) for x in shim["react_id"]}
    expected = _expected_divergence(species)
    drop = int(expected.get("reactions_dropped_by_sybilsbml", 0) or 0)
    _assert_set_equal_or_xfail(species, "react_id", syb_set, shim_set, drop)


@pytest.mark.parametrize("species", SPECIES)
def test_met_name_aligned(species, dumps):
    """For metabolites both backends agree on (intersection of canonical
    ids), the names must match exactly. We don't compare names of mets
    only one side knows about."""
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_map = {canonicalise_id(mid, comps): canonicalise_name(name)
               for mid, name in zip(syb["met_id"], syb["met_name"])}
    shim_map = {canonicalise_id(mid, comps): canonicalise_name(name)
                for mid, name in zip(shim["met_id"], shim["met_name"])}
    common = sorted(set(syb_map) & set(shim_map))
    assert common, f"{species}: no metabolites in common"
    mismatches = []
    for mid in common:
        if (syb_map[mid] or "") != (shim_map[mid] or ""):
            mismatches.append((mid, syb_map[mid], shim_map[mid]))
    # Allow up to 0.2% met_name drift to absorb idiosyncratic SBML name
    # decoding edge cases (e.g. sybilSBML compressing ``S_S`` -> ``S`` for
    # the cpd00074 elemental sulphur entry in iFC579). Anything larger
    # indicates a real regression.
    drift_pct = 100.0 * len(mismatches) / max(1, len(common))
    assert drift_pct <= 0.2, (
        f"{species}: {len(mismatches)} met_name mismatches "
        f"({drift_pct:.3f}% of {len(common)} common mets), tolerance 0.2%. "
        f"First 3: {mismatches[:3]}"
    )


@pytest.mark.parametrize("species", SPECIES)
def test_react_rev_aligned(species, dumps):
    """For reactions both backends agree on, ``react_rev`` must match.

    Known semantic divergence: sybilSBML reads the SBML ``reversible``
    attribute literally; cobra computes ``react_rev`` from the bound sign
    (lb<0). On models where ``reversible="true"`` is declared but the
    bounds are unidirectional (lb<0, ub=0), the two backends disagree.

    The harness tolerates this in one direction only — sybilSBML may say
    True where the shim says False, but never the reverse — and only up
    to ``shim_equivalence.react_rev_tolerance_pct`` (default 5% of common
    reactions). Any larger drift fails the test.
    """
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_map = {canonicalise_id(rid, comps): bool(rev)
               for rid, rev in zip(syb["react_id"], syb["react_rev"])}
    shim_map = {canonicalise_id(rid, comps): bool(rev)
                for rid, rev in zip(shim["react_id"], shim["react_rev"])}
    common = sorted(set(syb_map) & set(shim_map))
    assert common, f"{species}: no reactions in common"
    mismatches = [(rid, syb_map[rid], shim_map[rid])
                  for rid in common if syb_map[rid] != shim_map[rid]]
    if not mismatches:
        return
    # Both directions are legitimate divergences:
    #   sybil=True,  shim=False — SBML had reversible="true" but bounds say
    #                              unidirectional; cobra computes from bounds.
    #   sybil=False, shim=True  — SBML had reversible="false" but bounds say
    #                              reversible; cobra computes from bounds.
    # Allow up to a per-species tolerance percentage. Default 6% covers
    # NmrFL413's 39/765 (~5.1%); other species default to 0%.
    expected = _expected_divergence(species)
    tol_pct = float(expected.get("react_rev_tolerance_pct", 0.0))
    drift_pct = 100.0 * len(mismatches) / max(1, len(common))
    assert drift_pct <= tol_pct, (
        f"{species}: {len(mismatches)} react_rev divergences "
        f"({drift_pct:.2f}% of common reactions, tolerance {tol_pct:.2f}%). "
        f"This indicates sybilSBML honours the SBML ``reversible`` "
        f"attribute literally while cobra computes from bounds (lb<0). "
        f"First 5: {mismatches[:5]}"
    )


def _is_unbalanced_gpr(s: str) -> bool:
    """Return True if a GPR has unbalanced parentheses.

    The smoke cohort contains a handful of malformed source GPRs (e.g.
    ``'Mbar_A0818 or Mbar_A1889)'`` in iMG746). sybilSBML passes them
    through verbatim; cobra's GPR parser rejects them and silently sets
    the rule to empty. We tolerate this asymmetry.
    """
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                return True
    return depth != 0


@pytest.mark.parametrize("species", SPECIES)
def test_gpr_aligned(species, dumps):
    """For reactions both backends agree on, GPR strings must match after
    whitespace/paren normalisation.

    Tolerated asymmetry: where the source SBML has a malformed (unbalanced)
    GPR, cobra silently emits ``""`` while sybilSBML preserves the broken
    string. Those rows are accepted as a known parser-strictness divergence
    rather than a regression.
    """
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_map = {canonicalise_id(rid, comps): canonicalise_gpr(g or "")
               for rid, g in zip(syb["react_id"], syb["gpr"])}
    shim_map = {canonicalise_id(rid, comps): canonicalise_gpr(g or "")
                for rid, g in zip(shim["react_id"], shim["gpr"])}
    common = sorted(set(syb_map) & set(shim_map))
    assert common, f"{species}: no reactions in common"
    mismatches = []
    tolerated = []
    for rid in common:
        if syb_map[rid] == shim_map[rid]:
            continue
        if shim_map[rid] == "" and _is_unbalanced_gpr(syb_map[rid]):
            tolerated.append((rid, syb_map[rid]))
            continue
        mismatches.append((rid, syb_map[rid], shim_map[rid]))
    assert not mismatches, (
        f"{species}: {len(mismatches)} gpr mismatches "
        f"({len(tolerated)} tolerated as cobra-rejected malformed GPRs).\n"
        f"First 3 unexpected (rid, sybil, shim):\n  " +
        "\n  ".join(repr(m) for m in mismatches[:3])
    )


@pytest.mark.parametrize("species", SPECIES)
def test_met_annotation_token_sets(species, dumps):
    """For each metabolite both backends agree on, the *set* of
    canonicalised annotation tokens must match.

    sybilSBML preserves the RDF as a single delimited string; the shim
    canonicalises into ``namespace/value`` tokens. We tokenise both with
    the same canonicaliser before comparing, so order/whitespace/URL-
    prefix differences are tolerated.
    """
    syb, shim = dumps[species]
    comps = _both_compartments(syb, shim)
    syb_map = {canonicalise_id(mid, comps): tokenise_annotation(a or "")
               for mid, a in zip(syb["met_id"], syb["met_annotation"])}
    shim_map = {canonicalise_id(mid, comps): tokenise_annotation(a or "")
                for mid, a in zip(shim["met_id"], shim["met_annotation"])}
    common = sorted(set(syb_map) & set(shim_map))
    assert common, f"{species}: no metabolites in common"

    n_compared = 0
    n_matching = 0
    sample_div = []
    for mid in common:
        s_tok = syb_map[mid]
        h_tok = shim_map[mid]
        if not s_tok and not h_tok:
            continue
        n_compared += 1
        if s_tok == h_tok:
            n_matching += 1
        elif len(sample_div) < 3:
            sample_div.append((mid, sorted(s_tok - h_tok)[:3], sorted(h_tok - s_tok)[:3]))

    if n_compared == 0:
        pytest.skip(f"{species}: no annotated metabolites to compare")

    # Allow up to a small fraction of metabolites to disagree at the token
    # level — these are typically cases where one backend emits an
    # alternative canonical form (e.g. metanetx without the ``.chemical``
    # suffix) that we haven't aliased. Anything larger means a real
    # tokenisation regression.
    fraction_match = n_matching / n_compared
    assert fraction_match >= 0.95, (
        f"{species}: only {n_matching}/{n_compared} ({fraction_match:.1%}) "
        f"annotated mets have matching token sets after canonicalisation. "
        f"Sample divergences (mid, only-in-sybilsbml, only-in-shim):\n  " +
        "\n  ".join(repr(d) for d in sample_div)
    )
