"""Parity of :mod:`fabel.density` against golden output from R ``fda`` 6.3.0.

Every golden case is an ``intensity.fd`` fit (fda 6.3.0 no longer ships
``density.fd``).  The density cases use the exact identity between the two
problems: when the basis spans the constants and the penalty annihilates them,
the intensity optimum has ``∫ exp W = n`` and ``exp(W) / n`` is the penalised
density estimate.  So :func:`fit_density` is compared with R's
``exp(W) / n`` and :func:`fit_intensity` with R's ``exp(W)``.

Cases are parametrised by ``(case, field)``: R's quadrature of ``∫ exp W`` is
accurate to about 1e-5 (see ``R_FDA_DEFECTS``), which is inside the iterative
tolerance for every intensity and density value but not for the few fields
that are relative comparisons of numbers close to zero.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import Monomial
from fabel.basis import Basis
from fabel.density import _Problem, fit_density, fit_intensity

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "density"

_INTENSITY_FIELDS = ("coefs", "log_intensity", "intensity", "f")
_DENSITY_FIELDS = ("density", "log_density")

_QUADRATURE = (
    "R's intensity.fd integrates exp(W) with a quadrature that is off by up to "
    "1.6e-4 in its reported criterion f; its coefficients are therefore not "
    "stationary for the exact criterion, and Fabel's (exact Gauss-Legendre "
    "integral, gradient norm < 2e-13) give a lower exact criterion in every "
    "case (test_fabel_improves_on_r_criterion). "
)

#: ``(case, field)`` pairs where R is the less accurate side.  Measured with
#: this golden file: the absolute error is at R's quadrature level, but the
#: entry it lands on is close to zero, so the entrywise relative error exceeds
#: the iterative tolerance 1e-5.  Strict: a regeneration that fixes one fails.
R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    ("intensity_sine_bspline11_order3_L1_lam0", "coefs"): _QUADRATURE
    + "Exact gradient norm at R's coefficients 4.07e-5.  Worst entry: "
    "coefficient 7, R 0.11139673524, Fabel 0.11140196492 (abs 5.2e-6, rel 4.7e-5).",
    ("density_regina_bspline13_L2_lam0.1", "coefs"): _QUADRATURE
    + "Exact gradient norm at R's coefficients 4.18e-5.  Worst entry: "
    "coefficient 6, R -0.037409071057, Fabel -0.037405675073 (abs 3.4e-6, rel 9.1e-5).",
    ("density_regina_bspline13_L2_lam0.1", "log_intensity"): _QUADRATURE
    + "Worst entry: W(42.85), R 0.033795023457, Fabel 0.033793807749 "
    "(abs 1.2e-6, rel 3.6e-5); the intensity itself agrees to 4.3e-6.",
    ("density_gamma_bspline9_L2_lam1", "log_intensity"): _QUADRATURE
    + "Exact gradient norm at R's coefficients 1.60e-5.  Worst entry: "
    "W(8.7), R 0.010768158864, Fabel 0.010768044513 (abs 1.1e-7, rel 1.1e-5).",
    ("density_truncnormal_monomial3_L1_lam0", "log_intensity"): _QUADRATURE
    + "Exact gradient norm at R's coefficients 5.16e-5.  Worst entry: "
    "W(-2.1), R -0.033503115309, Fabel -0.033503752936 (abs 6.4e-7, rel 1.9e-5).",
}


def _basis(spec: dict[str, Any]) -> Basis:
    """Golden basis, including R's monomial type (not in the shared helper)."""
    if spec["type"] == "monom":
        domain = (float(spec["rangeval"][0]), float(spec["rangeval"][1]))
        return Monomial(domain=domain, n_basis=spec["nbasis"])
    return build_basis(spec)


def _params() -> list[Any]:
    params = []
    for case in golden_cases(MODULE):
        fields: tuple[str, ...] = _INTENSITY_FIELDS
        if case["input"]["kind"] == "density":
            fields = fields + _DENSITY_FIELDS
        for field in fields:
            key = (case["name"], field)
            marks = []
            if key in R_FDA_DEFECTS:
                marks.append(pytest.mark.xfail(reason=R_FDA_DEFECTS[key], strict=True))
            params.append(pytest.param(case, field, id=f"{case['name']}-{field}", marks=marks))
    return params


def _fit_args(case: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    inp = case["input"]
    x = np.asarray(inp["x"], dtype=float)
    return x, {"basis": _basis(inp["basis"]), "lam": float(inp["lambda"]), "penalty": inp["lfd"]}


@pytest.mark.parametrize(("case", "field"), _params())
def test_density_parity(case: dict[str, Any], field: str) -> None:
    x, kwargs = _fit_args(case)
    out = case["output"]
    t = np.asarray(out["t"], dtype=float)
    rtol = case_rtol(MODULE, case)
    if field in _DENSITY_FIELDS:
        dens = fit_density(x, **kwargs)
        actual = dens(t) if field == "density" else dens.log_density(t)
    else:
        res = fit_intensity(x, **kwargs)
        if field == "coefs":
            actual = res.fd.coefs[:, 0]
        elif field == "log_intensity":
            actual = res.fd(t)[:, 0]
        elif field == "intensity":
            actual = res(t)
        else:
            actual = res.criterion
    compare(actual, out[field], rtol)


@pytest.mark.parametrize("case", golden_cases(MODULE), ids=lambda c: c["name"])
def test_fabel_improves_on_r_criterion(case: dict[str, Any]) -> None:
    """R's coefficients are near-optimal, and never better, for the exact criterion."""
    x, kwargs = _fit_args(case)
    res = fit_intensity(x, **kwargs)
    basis = kwargs["basis"]
    problem = _Problem("intensity", basis, x, basis.penalty(kwargs["penalty"]), kwargs["lam"])
    r_value = problem.value(np.asarray(case["output"]["coefs"], dtype=float))
    gap = r_value - res.criterion
    assert 0.0 <= gap <= 1e-9 * max(1.0, abs(res.criterion))


def test_regina_sample_is_the_dataset() -> None:
    """The golden Regina sample is ReginaPrecip restricted to (2, 45] mm, sorted."""
    datasets = pytest.importorskip("fabel.datasets")
    try:
        value = np.asarray(datasets.load_regina_precip().value, dtype=float)
    except (OSError, ValueError) as exc:  # offline and not cached
        pytest.skip(f"ReginaPrecip is not available: {exc}")
    rain = np.sort(value[(value > 2.0) & (value <= 45.0)])
    case = next(c for c in golden_cases(MODULE) if c["input"].get("dataset") == "ReginaPrecip")
    np.testing.assert_array_equal(rain, np.asarray(case["input"]["x"], dtype=float))
