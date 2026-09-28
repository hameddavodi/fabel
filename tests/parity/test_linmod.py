"""Parity of :func:`fdatools.regression.linmod` against R ``fda`` 6.3.0's ``linmod``.

``tests/golden/linmod.json`` (generator ``tools/golden_r/linmod.R``) records
R's intercept, regression surface and fitted curves on the weather data
(log10 precipitation on temperature) and on three synthetic designs.  Every
field is asserted on its own.

R evaluates the integrals of the model numerically.  On the monomial designs
those rules are exact, and fdatools matches R to rounding on every field.  On the
Fourier and B-spline designs R's integrals are off in the 6th and 5th digit;
the affected fields are strict xfails carrying the measured error.  The
synthetic cases also record the integral and penalty matrices R's ``inprod``
and ``eval.penalty`` return, and feeding *those* into fdatools' normal equations
reproduces R's coefficients to 1e-8 -- so the model is the same and only the
quadrature differs.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache
from typing import Any

import numpy as np
import pytest

from fdatools import LDO, BSpline, FData, Monomial, inprod
from fdatools._backend import default_namespace
from fdatools.basis import Basis
from fdatools.regression import (
    LinmodResult,
    _linmod_normal_equations,
    _linmod_solve,
    _LinmodIntegrals,
    linmod,
)

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "linmod"
CASES = {case["name"]: case for case in golden_cases(MODULE)}
FIELDS = ("beta0estfd_coefs", "beta1estbifd_coefs", "yhatfdobj_coefs")

#: Points R's ``linmod`` evaluates the fitted curves on before refitting them in
#: the response basis by least squares (found by matching the golden output).
R_YHAT_GRID = 201

_WEATHER_QUADRATURE = (
    "R's integrals over the Fourier bases are approximate.  Its harmonic-"
    "accelerator penalty of the Fourier(11) coefficient basis (eval.penalty) is "
    "1.2e-4 (relative) off the closed form k^2 w^6 (k^2 - 1)^2, and its Fourier "
    "Gram inprod(fb11, fb11) 1.4e-6 off the identity.  fdatools' exact fit is "
    "{alpha} (intercept), {beta} (surface) and {yhat} (fitted curves) off R, "
    "relative to the largest entry."
)
_BSPLINE_QUADRATURE = (
    "R's inprod of cubic B-spline bases is 2.2e-5 (Gram) and 1.1e-4 (the "
    "covariate integrals int x_i theta_s) off the exact values; fdatools' exact "
    "intercept and surface are 5.9e-5 and 9.6e-5 off R.  With R's own matrices "
    "(r_integrals) fdatools reproduces R's coefficients to 1e-15 "
    "(test_normal_equations_with_r_integrals).  R's fitted curves are further a "
    "201-point least-squares fit, not the L2 projection: 1.4e-3 off fdatools "
    "(test_r_fitted_curves_are_a_grid_least_squares_fit)."
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    **{
        ("linmod_weather_logprecip_on_temp", field): _WEATHER_QUADRATURE.format(
            alpha="4.7e-7", beta="1.9e-6", yhat="3.7e-7"
        )
        for field in FIELDS
    },
    **{("linmod_synthetic_bspline", field): _BSPLINE_QUADRATURE for field in FIELDS},
}


# --------------------------------------------------------------------------- #
# model construction
# --------------------------------------------------------------------------- #


def _basis(spec: dict[str, Any]) -> Basis:
    """Build a golden basis; R names the monomial type ``monom`` (conftest knows ``monomial``)."""
    if spec["type"] == "monom":
        return Monomial(domain=(spec["rangeval"][0], spec["rangeval"][1]), n_basis=spec["nbasis"])
    return build_basis(spec)


def _penalty(inputs: dict[str, Any]) -> LDO:
    if inputs["penalty"] == "harmonic":
        return LDO.harmonic(period=float(inputs["period"]))
    return LDO(int(inputs["penalty"]))


def _curves(inputs: dict[str, Any], which: str) -> FData:
    return FData(
        np.asarray(inputs[f"{which}_coefs"], dtype=float), _basis(inputs[f"{which}_basis"])
    )


@cache
def _model(name: str) -> LinmodResult:
    inputs = CASES[name]["input"]
    op = _penalty(inputs)
    return linmod(
        _curves(inputs, "y"),
        _curves(inputs, "x"),
        alpha_basis=_basis(inputs["alpha_basis"]),
        s_basis=_basis(inputs["s_basis"]),
        t_basis=_basis(inputs["t_basis"]),
        lam_alpha=inputs["lambda_alpha"],
        lam_s=inputs["lambda_s"],
        lam_t=inputs["lambda_t"],
        penalty_alpha=op,
        penalty_s=op,
        penalty_t=op,
    )


_ACCESSORS: dict[str, Callable[[LinmodResult], Any]] = {
    "beta0estfd_coefs": lambda model: model.alpha.coefs[:, 0],
    "beta1estbifd_coefs": lambda model: model.beta.coefs,
    "yhatfdobj_coefs": lambda model: model.fitted.coefs,
}


def _params() -> list[Any]:
    params = []
    for name in CASES:
        for field in FIELDS:
            marks = []
            reason = R_FDA_DEFECTS.get((name, field))
            if reason is not None:
                marks.append(pytest.mark.xfail(reason=reason, strict=True))
            params.append(pytest.param(name, field, id=f"{name}-{field}", marks=marks))
    return params


# --------------------------------------------------------------------------- #
# golden parity
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "field"), _params())
def test_golden_field(name: str, field: str) -> None:
    case = CASES[name]
    want = np.asarray(case["output"][field], dtype=float)
    actual = np.asarray(_ACCESSORS[field](_model(name)), dtype=float)
    compare(np.reshape(actual, want.shape), want, case_rtol(MODULE, case))


def test_every_defect_names_a_real_field() -> None:
    known = {(name, field) for name in CASES for field in FIELDS}
    assert set(R_FDA_DEFECTS) <= known


def _relative(actual: Any, expected: Any) -> float:
    want = np.asarray(expected, dtype=float)
    got = np.reshape(np.asarray(actual, dtype=float), want.shape)
    return float(np.max(np.abs(got - want)) / np.max(np.abs(want)))


#: Relative error of fdatools against R quoted in the xfail reasons, as (low, high).
MEASURED = {
    ("linmod_weather_logprecip_on_temp", "beta0estfd_coefs"): (4e-7, 6e-7),
    ("linmod_weather_logprecip_on_temp", "beta1estbifd_coefs"): (1.5e-6, 2.5e-6),
    ("linmod_weather_logprecip_on_temp", "yhatfdobj_coefs"): (3e-7, 5e-7),
    ("linmod_synthetic_bspline", "beta0estfd_coefs"): (5e-5, 7e-5),
    ("linmod_synthetic_bspline", "beta1estbifd_coefs"): (8e-5, 1.2e-4),
    ("linmod_synthetic_bspline", "yhatfdobj_coefs"): (1e-3, 2e-3),
}


@pytest.mark.parametrize(("name", "field"), sorted(MEASURED))
def test_defect_is_as_measured(name: str, field: str) -> None:
    """Pin the size of R's quadrature error quoted in the xfail reasons."""
    low, high = MEASURED[name, field]
    error = _relative(_ACCESSORS[field](_model(name)), CASES[name]["output"][field])
    assert low < error < high


# --------------------------------------------------------------------------- #
# the model without the quadrature: R's own integrals
# --------------------------------------------------------------------------- #

_WITH_INTEGRALS = sorted(name for name, case in CASES.items() if "r_integrals" in case["output"])


def _r_integrals(name: str) -> tuple[np.ndarray, _LinmodIntegrals]:
    raw = {
        key: np.asarray(value, dtype=float)
        for key, value in CASES[name]["output"]["r_integrals"].items()
    }
    z = raw.pop("Z")
    return z, _LinmodIntegrals(**{key.lower(): value for key, value in raw.items()})


@pytest.mark.parametrize("name", _WITH_INTEGRALS)
def test_normal_equations_with_r_integrals(name: str) -> None:
    """fdatools' normal equations on R's integral matrices give R's coefficients."""
    case = CASES[name]
    inputs = case["input"]
    xp = default_namespace()
    z, integrals = _r_integrals(name)
    ycoefs = np.asarray(inputs["y_coefs"], dtype=float)
    lams = (inputs["lambda_alpha"], inputs["lambda_s"], inputs["lambda_t"])
    weights = np.ones(ycoefs.shape[1])
    cmat, dmat = _linmod_normal_equations(z, ycoefs, weights, integrals, lams, xp)
    shape = (integrals.gss.shape[0], integrals.gtt.shape[0])
    alpha, surface = _linmod_solve(cmat, dmat, integrals.gaa.shape[0], shape, xp)
    rtol = case_rtol(MODULE, case)
    compare(alpha, case["output"]["beta0estfd_coefs"], rtol)
    compare(surface, case["output"]["beta1estbifd_coefs"], rtol)


@pytest.mark.parametrize("name", [name for name in _WITH_INTEGRALS if "monomial" in name])
def test_r_integrals_match_fdatools_where_r_is_exact(name: str) -> None:
    """On the monomial designs R's integrals are exact, so fdatools' agree."""
    inputs = CASES[name]["input"]
    z, integrals = _r_integrals(name)
    alpha, sbasis, tbasis = (_basis(inputs[f"{k}_basis"]) for k in ("alpha", "s", "t"))
    op = _penalty(inputs)
    exact = _LinmodIntegrals.exact(
        _basis(inputs["y_basis"]), alpha, sbasis, tbasis, (op, op, op), default_namespace()
    )
    x = _curves(inputs, "x")
    compare(np.asarray(x.coefs).T @ np.asarray(inprod(x.basis, sbasis)), z, 1e-12)
    for field in ("gaa", "gat", "gtt", "gss", "gay", "gty", "ra", "rs", "rt"):
        compare(getattr(exact, field), getattr(integrals, field), 1e-12)


def test_r_fitted_curves_are_a_grid_least_squares_fit() -> None:
    """R's yhat is its fit evaluated on 201 points and refitted by least squares.

    Rebuilding it from R's own coefficients and integrals that way reproduces
    R's recorded fitted curves, which pins the second half of the B-spline xfail
    reason.
    """
    name = "linmod_synthetic_bspline"
    case = CASES[name]
    inputs = case["input"]
    z, _ = _r_integrals(name)
    ybasis, alpha_basis, tbasis = (_basis(inputs[f"{k}_basis"]) for k in ("y", "alpha", "t"))
    assert isinstance(ybasis, BSpline)
    grid = np.linspace(*ybasis.domain, R_YHAT_GRID)
    alpha = np.asarray(case["output"]["beta0estfd_coefs"], dtype=float)
    surface = np.asarray(case["output"]["beta1estbifd_coefs"], dtype=float)
    values = (np.asarray(alpha_basis(grid)) @ alpha)[:, None] + np.asarray(tbasis(grid)) @ (
        surface.T @ z.T
    )
    coefs = np.linalg.lstsq(np.asarray(ybasis(grid)), values, rcond=None)[0]
    compare(coefs, case["output"]["yhatfdobj_coefs"], case_rtol(MODULE, case))
