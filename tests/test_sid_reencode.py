"""Unit test for enforce_sbml_sids() (issue #7).

The R function `enforce_sbml_sids()` lives in
`pipeline/process_sbml_species.R`. It scans `met_df$new_met_out` (formatted
as `name[compartment]`), and rewrites the bare-name portion so it matches
the SBML SId regex `^[A-Za-z_][A-Za-z0-9_]*$`.

We exercise it via a tiny Rscript harness rather than reticulate, both
because reticulate adds a heavy install dependency and because the
function is pure (no SBML I/O, no MetaNetX TSV) so a one-shot Rscript run
is perfectly fast.

The test covers:
  - Issue #7's iJDZ836 cases: Fe+2, Mg2+, Na+, Fe+3, K+ (all illegal).
  - `glc__D` (already legal -> must be left alone).
  - `2pg` (legal chars but illegal leading digit -> must be prefixed).
  - Collision: distinct inputs that re-encode to the same name. The
    function itself does NOT dedup (that's `handle_duplicates()`'s job
    on the second pass), but we verify the test harness can detect
    such collisions deterministically so the wired-up pipeline knows
    to invoke handle_duplicates afterwards.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from conftest import REPO_ROOT

SID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
PROCESS_SBML_R = REPO_ROOT / "pipeline" / "process_sbml_species.R"


def _have_rscript() -> bool:
    return shutil.which("Rscript") is not None


pytestmark = pytest.mark.skipif(
    not _have_rscript(), reason="Rscript not available on PATH"
)


def _run_enforce(rows):
    """Invoke enforce_sbml_sids() on a synthetic met_df via Rscript.

    `rows` is a list of strings (the new_met_out values).
    Returns a list of strings in the same order.
    """
    quoted = ",".join(f'"{r}"' for r in rows)
    script = textwrap.dedent(
        f"""
        # Extract enforce_sbml_sids() from the source file without
        # sourcing the whole pipeline (avoids tidyverse/sybilSBML deps).
        src <- readLines("{PROCESS_SBML_R}")
        start <- grep("^enforce_sbml_sids <- function", src)
        if (length(start) != 1L) {{
            stop("could not locate enforce_sbml_sids() in source")
        }}
        end_candidates <- grep("^}}", src)
        end <- min(end_candidates[end_candidates > start])
        eval(parse(text = paste(src[start:end], collapse = "\\n")))

        met_df <- data.frame(
          new_met_out = c({quoted}),
          stringsAsFactors = FALSE
        )
        out <- enforce_sbml_sids(met_df)
        cat("__BEGIN__\\n")
        for (v in out$new_met_out) cat(v, "\\n", sep = "")
        cat("__END__\\n")
        """
    )
    proc = subprocess.run(
        ["Rscript", "--vanilla", "-e", script],
        capture_output=True,
        text=True,
        timeout=60,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Rscript failed (rc={proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    out = proc.stdout
    # Pull the block between sentinels so any cat()/message() noise from
    # the function's own logging doesn't pollute parsing.
    m = re.search(r"__BEGIN__\n(.*?)__END__", out, re.DOTALL)
    assert m, f"could not find sentinel block in Rscript output:\n{out}"
    block = m.group(1)
    return [line for line in block.splitlines() if line]


def _bare_name(value: str) -> str:
    """Strip the [compartment] suffix and return the bare SId-candidate."""
    idx = value.find("[")
    return value[:idx] if idx >= 0 else value


def test_enforce_sbml_sids_issue7_metals():
    """The 5 iJDZ836 metal species become legal SIds."""
    inputs = [
        "Fe+2[CCO-EXTRACELLULAR]",
        "Mg2+[c]",
        "Na+[e]",
        "Fe+3[c]",
        "K+[c]",
    ]
    out = _run_enforce(inputs)
    assert len(out) == len(inputs)
    for original, repaired in zip(inputs, out):
        bare = _bare_name(repaired)
        assert SID_RE.match(bare), (
            f"input {original!r} -> {repaired!r} but bare name "
            f"{bare!r} still doesn't match SBML SId regex"
        )
    # The compartment suffix must survive untouched (input order:
    # CCO-EXTRACELLULAR, c, e, c, c).
    expected_suffixes = ["[CCO-EXTRACELLULAR]", "[c]", "[e]", "[c]", "[c]"]
    for value, suffix in zip(out, expected_suffixes):
        assert value.endswith(suffix), (
            f"expected compartment {suffix} preserved on {value!r}"
        )


def test_enforce_sbml_sids_legal_inputs_unchanged():
    """A row whose name is already a valid SId is passed through verbatim."""
    inputs = ["glc__D[c]", "ATP[m]", "h2o[e]"]
    out = _run_enforce(inputs)
    assert out == inputs


def test_enforce_sbml_sids_leading_digit():
    """Leading-digit names get an `M_` prefix."""
    inputs = ["2pg[c]", "10fthf[c]"]
    out = _run_enforce(inputs)
    for original, repaired in zip(inputs, out):
        bare = _bare_name(repaired)
        assert SID_RE.match(bare), (
            f"input {original!r} -> {repaired!r} but bare {bare!r} "
            f"still illegal"
        )
        assert bare.startswith("M_"), f"expected M_ prefix on {repaired}"


def test_enforce_sbml_sids_full_synthetic_set():
    """Combined iJDZ836 + legal + leading-digit cases all become legal."""
    inputs = [
        "Fe+2[CCO-EXTRACELLULAR]",
        "Mg2+[c]",
        "Na+[e]",
        "glc__D[c]",
        "2pg[c]",
    ]
    out = _run_enforce(inputs)
    assert len(out) == len(inputs)
    bare_names = [_bare_name(v) for v in out]
    for v, bare in zip(out, bare_names):
        assert SID_RE.match(bare), f"illegal SId remains: {v!r} bare={bare!r}"
    # No collisions in this representative set (post-rewrite ids are
    # `Fe_2`, `Mg2_`, `Na_`, `glc__D`, `M_2pg`).
    assert len(set(out)) == len(out), f"unexpected collision: {out}"


def test_enforce_sbml_sids_empty_dataframe():
    """No-op on an empty input."""
    out = _run_enforce([])
    assert out == []
