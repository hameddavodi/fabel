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

import numpy as np

from fabel import BSpline, Constant, Exponential, Fourier, Monomial
from fabel.basis import Basis

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


def build_basis(spec: dict[str, Any]) -> Basis:
    """Construct the Fabel basis described by a golden basis specification.

    Golden files record a basis as R's ``basisfd`` summary: ``type``,
    ``rangeval``, ``nbasis`` and the type-specific ``params`` (interior knots
    for a spline, the period for a Fourier basis, the rates for an exponential
    one).  A B-spline's order is recovered as ``nbasis - n_interior_knots``.
    """
    kind = spec["type"]
    domain = (float(spec["rangeval"][0]), float(spec["rangeval"][1]))
    params: Any = spec.get("params")
    if kind == "bspline":
        interior = [] if params is None else np.atleast_1d(np.asarray(params, dtype=float)).tolist()
        breaks = [domain[0], *interior, domain[1]]
        return BSpline(domain=domain, order=spec["nbasis"] - len(interior), breaks=breaks)
    if kind == "fourier":
        return Fourier(domain=domain, n_basis=spec["nbasis"], period=float(params))
    if kind == "monomial":
        return Monomial(domain=domain, n_basis=spec["nbasis"])
    if kind == "expon":
        return Exponential(domain=domain, rates=np.atleast_1d(params).tolist())
    if kind == "const":
        return Constant(domain=domain)
    raise AssertionError(f"unknown golden basis type {kind!r}")


def compare(actual: Any, expected: Any, rtol: float) -> None:
    """Compare two arrays, scaling ``atol`` by the magnitude of the expected values."""
    want = np.asarray(expected, dtype=float)
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), want, rtol=rtol, atol=1e-12 * max(1.0, scale)
    )
