"""Generate ``tests/golden/<module>.json`` by driving R fda 6.3.0.

Usage::

    .venv/bin/python tools/make_golden.py <module>

Drives R through ``rpy2`` when importable, else falls back to an
``Rscript`` subprocess running the same driver snippet. Either path
sources ``tools/golden_r/common.R`` (helpers + seeding + JSON writer) and
then ``tools/golden_r/<module>.R`` (the module-specific golden cases), so
the R code itself is identical regardless of which path drives it.

See ``docs/dev/conventions.md`` for the golden-file JSON schema.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

SEED = 20260906
REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN_R_DIR = REPO_ROOT / "tools" / "golden_r"
GOLDEN_OUT_DIR = REPO_ROOT / "tests" / "golden"


def _driver_snippet(module: str) -> str:
    """Build the R snippet that seeds RNG and sources common.R + <module>.R."""
    common_r = GOLDEN_R_DIR / "common.R"
    module_r = GOLDEN_R_DIR / f"{module}.R"
    if not module_r.exists():
        raise FileNotFoundError(f"no R script for module {module!r}: {module_r}")
    out_json = GOLDEN_OUT_DIR / f"{module}.json"
    return (
        f'GOLDEN_OUT <- "{out_json.as_posix()}"\n'
        f"GOLDEN_SEED <- {SEED}\n"
        f'source("{common_r.as_posix()}")\n'
        f'source("{module_r.as_posix()}")\n'
    )


def _run_via_rpy2(snippet: str) -> None:
    """Execute the driver snippet through rpy2's embedded R interpreter."""
    import rpy2.robjects as robjects

    robjects.r(snippet)


def _run_via_rscript(snippet: str) -> None:
    """Execute the driver snippet via an ``Rscript`` subprocess."""
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as fh:
        fh.write(snippet)
        driver_path = fh.name
    try:
        result = subprocess.run(
            ["Rscript", "--vanilla", driver_path],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            sys.stderr.write(result.stdout)
            sys.stderr.write(result.stderr)
            raise RuntimeError(f"Rscript failed with exit code {result.returncode}")
    finally:
        Path(driver_path).unlink(missing_ok=True)


def generate(module: str) -> Path:
    """Generate ``tests/golden/<module>.json`` and return its path."""
    GOLDEN_OUT_DIR.mkdir(parents=True, exist_ok=True)
    snippet = _driver_snippet(module)
    try:
        import rpy2.robjects  # noqa: F401
    except ImportError:
        _run_via_rscript(snippet)
    else:
        _run_via_rpy2(snippet)
    out_path = GOLDEN_OUT_DIR / f"{module}.json"
    if not out_path.exists():
        raise RuntimeError(f"expected output not written: {out_path}")
    return out_path


def main() -> int:
    """CLI entry point: generate the golden file for one module."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("module", help="module name, e.g. 'basis' or 'core'")
    args = parser.parse_args()

    out_path = generate(args.module)
    data = json.loads(out_path.read_text())
    n_cases = len(data.get("cases", []))
    meta = data.get("meta", {})
    print(f"wrote {out_path} ({n_cases} cases)")
    print(f"meta: {meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
