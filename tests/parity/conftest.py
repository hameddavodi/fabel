"""Golden-file loading helpers for parity tests.

Golden files live at ``tests/golden/<module>.json`` and are generated only
by ``tools/make_golden.py`` (see ``docs/dev/conventions.md`` for the schema).
Parity tests (``tests/parity/test_<module>.py``) load expected R output via
``golden`` / ``golden_cases`` and compare against ``fabel`` with
``numpy.testing.assert_allclose(actual, expected, rtol=case_rtol, atol=1e-12)``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

GOLDEN_DIR = Path(__file__).resolve().parent.parent / "golden"


def golden(module: str) -> dict[str, Any]:
    """Load and return the full parsed ``tests/golden/<module>.json`` document."""
    path = GOLDEN_DIR / f"{module}.json"
    with path.open(encoding="utf-8") as fh:
        doc: dict[str, Any] = json.load(fh)
        return doc


def golden_cases(module: str) -> list[dict[str, Any]]:
    """Return the ``cases`` list from ``tests/golden/<module>.json``."""
    cases: list[dict[str, Any]] = golden(module)["cases"]
    return cases


def case_rtol(module: str, case: dict[str, Any]) -> float:
    """Return the tolerance to use for one case.

    Honours a per-case ``"rtol"`` override (e.g. iterative fits like
    ``smooth.monotone``/``smooth.pos`` need a looser tolerance than the
    module default); falls back to ``meta["rtol"]`` when absent.
    """
    if "rtol" in case:
        rtol: float = case["rtol"]
        return rtol
    meta_rtol: float = golden(module)["meta"]["rtol"]
    return meta_rtol
