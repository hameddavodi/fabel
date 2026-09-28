"""Parity of :mod:`fabel.sparse` against golden output from R ``fda`` 6.3.0.

The golden file ``tests/golden/pace.json`` holds seeded sparse subsamples of
two real datasets (6 random ages for each of the 54 girls of ``growth$hgtf``;
10 random days for each of the 35 ``CanadianWeather`` stations) together with
R's ``smooth.sparse.mean``, ``covPACE``, ``pcaPACE`` and ``scoresPACE``
output.  Cases are parametrised by ``(case, field)``, as in
``test_decomposition.py``, so a defect in one R field cannot hide agreement in
another.

Signs
-----
R's ``pcaPACE`` signs each harmonic so that its coefficients in orthonormal
coordinates, ``U b`` with ``W + λR = UᵀU`` (upper Cholesky factor), sum to a
positive number.  Measured on 90 harmonics of random sparse weather subsamples
this held every time, while a positive plain coefficient sum held for only
76-79 %; the weather case of the golden file is one of the exceptions.  Fabel
uses the same rule, so harmonics are compared as they are.
"""

from __future__ import annotations

import warnings
from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel import LDO
from fabel.basis import Basis
from fabel.sparse import PACE, _harmonics, sparse_mean

from .conftest import build_basis as _build_basis
from .conftest import case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "pace"

# --------------------------------------------------------------------------- #
# R defects
# --------------------------------------------------------------------------- #

_INPROD_HARMONICS_REASON = (
    "R's pcaPACE builds the covariance operator J C J' with J = inprod("
    "harmonic basis, covariance basis), a Romberg quadrature good to 4-5 "
    "digits, and normalises with the exact eval.penalty Gram.  Measured "
    "max|inprod(b, b) - eval.penalty(b, 0)|: 4.0e-5 for the order-4 spline "
    "with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on "
    "[1, 365]; it is exactly 0 for monomial bases.  With R's own inprod() J "
    "(recorded in the golden file) Fabel's eigenproblem reproduces R's values "
    "to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram).  "
    "With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 "
    "(values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 "
    "/ growth cov-lambda 10 / weather.  Fabel uses the exact J."
)

_INPROD_COV_PENALTY_REASON = (
    "R's covPACE penalty is lambda (W (x) P + P (x) W) with W = inprod(basis, "
    "basis) (Romberg, 4.0e-5 off the exact Gram for this order-4 spline on "
    "[1, 18]) and the exact P = eval.penalty(basis, 2).  With R's inprod() W "
    "(recorded in the golden file) Fabel reproduces R's surface to 1.9e-14 "
    "(test_same_math_with_r_inprod_gram); with the exact W the relative gap is "
    "1.2e-5.  Unpenalised surfaces (cov lambda 0) match to 6e-14."
)

_SCORES_REASON = (
    "R's scoresPACE does not compute the conditional expectation "
    "Lambda Xi_i' Sigma_i^-1 (y_i - mu_i).  Measured as a black box: every "
    "finite row of R's score matrix is exactly parallel (to 1e-12) to "
    "lambda_k * xi_k(tau) for ONE integer day tau (test_r_scores_use_one_point), "
    "so each curve's score vector has rank one in the harmonics, and 33 of the "
    "35 rows are NA.  On a dense regular probe (20 points, all curves at t = "
    "1..20) curve i used tau = i, and a unit change of the first observation "
    "moved the scores about 365 times more than the same change of the fifth; "
    "changing another curve's data moved them too.  The same function errors "
    "('evalarg contains 1 NA') for non-integer times such as growth ages.  "
    "Fabel implements the published PACE estimator (Yao, Mueller & Wang 2005)."
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    ("pace_growth_bspline6_covlambda0", "values"): _INPROD_HARMONICS_REASON,
    ("pace_growth_bspline6_covlambda0", "varprop"): _INPROD_HARMONICS_REASON,
    ("pace_growth_bspline6_covlambda0", "harmonics"): _INPROD_HARMONICS_REASON,
    ("pace_growth_bspline6_covlambda10", "cov_coefs"): _INPROD_COV_PENALTY_REASON,
    ("pace_growth_bspline6_covlambda10", "values"): _INPROD_HARMONICS_REASON,
    ("pace_growth_bspline6_covlambda10", "varprop"): _INPROD_HARMONICS_REASON,
    ("pace_growth_bspline6_covlambda10", "harmonics"): _INPROD_HARMONICS_REASON,
    ("pace_weather_fourier5_bspline7", "values"): _INPROD_HARMONICS_REASON,
    ("pace_weather_fourier5_bspline7", "varprop"): _INPROD_HARMONICS_REASON,
    ("pace_weather_fourier5_bspline7", "harmonics"): _INPROD_HARMONICS_REASON,
    ("pace_weather_fourier5_bspline7", "scores"): _SCORES_REASON,
}

_COMPARED = ("mean_coefs", "cov_coefs", "values", "varprop", "harmonics", "scores")

# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def build_basis(spec: dict[str, Any]) -> Basis:
    """Build a golden basis; R names a monomial basis ``"monom"`` here."""
    if spec["type"] == "monom":
        spec = {**spec, "type": "monomial"}
    return _build_basis(spec)


@cache
def _cases() -> dict[str, dict[str, Any]]:
    return {case["name"]: case for case in golden_cases(MODULE)}


def _data(case: dict[str, Any]) -> tuple[list[np.ndarray], list[np.ndarray]]:
    inp = case["input"]
    times = [np.atleast_1d(np.asarray(v, dtype=float)) for v in inp["t"]]
    values = [np.atleast_1d(np.asarray(v, dtype=float)) for v in inp["y"]]
    return times, values


@cache
def _fit(name: str) -> PACE:
    case = _cases()[name]
    inp = case["input"]
    times, values = _data(case)
    model = PACE(
        n=inp["nharm"],
        basis=build_basis(inp["cov_basis"]),
        mean_basis=build_basis(inp["mean_basis"]),
        harmonic_basis=build_basis(inp["harm_basis"]),
        lam_mean=inp["mean_lambda"],
        lam_cov=inp["cov_lambda"],
        lam=inp["harm_lambda"],
    )
    # These dense-ish, low-noise subsamples give a non-positive sigma^2
    # estimate; the floor warning is expected and irrelevant to R's fields.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return model.fit(values, t=times)


def _actual(name: str, field: str) -> np.ndarray:
    model = _fit(name)
    if field == "mean_coefs":
        return np.asarray(model.mean_fd.coefs)[:, 0]
    if field == "cov_coefs":
        return np.asarray(model.cov.coefs)
    if field == "harmonics":
        return np.asarray(model.harmonics.coefs)
    return np.asarray(getattr(model, field))


def _params(cases: list[dict[str, Any]]) -> list[Any]:
    params = []
    for case in cases:
        for field in _COMPARED:
            if field not in case["output"]:
                continue
            reason = R_FDA_DEFECTS.get((case["name"], field))
            marks = [pytest.mark.xfail(strict=True, reason=reason)] if reason else []
            params.append(pytest.param(case, field, marks=marks, id=f"{case['name']}-{field}"))
    return params


_PACE_CASES = [c for c in golden_cases(MODULE) if c["name"].startswith("pace_")]
_MEAN_CASES = [c for c in golden_cases(MODULE) if c["name"].startswith("smooth_sparse_mean")]

# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("case", _MEAN_CASES, ids=lambda c: c["name"])
def test_sparse_mean(case: dict[str, Any]) -> None:
    times, values = _data(case)
    basis = build_basis(case["output"]["basis"])
    mean = sparse_mean(values, times, basis, lam=case["input"]["lambda"])
    compare(np.asarray(mean.coefs)[:, 0], case["output"]["coefs"], case_rtol(MODULE, case))


def test_sparse_mean_accepts_r_list_form() -> None:
    case = _MEAN_CASES[0]
    times, values = _data(case)
    pairs = [np.column_stack([t, y]) for t, y in zip(times, values, strict=True)]
    basis = build_basis(case["output"]["basis"])
    mean = sparse_mean(pairs, None, basis, lam=case["input"]["lambda"])
    compare(np.asarray(mean.coefs)[:, 0], case["output"]["coefs"], case_rtol(MODULE, case))


@pytest.mark.parametrize(("case", "field"), _params(_PACE_CASES))
def test_pace_field(case: dict[str, Any], field: str) -> None:
    compare(_actual(case["name"], field), case["output"][field], case_rtol(MODULE, case))


@pytest.mark.parametrize("case", _PACE_CASES, ids=lambda c: c["name"])
def test_same_math_with_r_inprod_gram(case: dict[str, Any]) -> None:
    """With R's inprod() matrices in place of the exact ones, R's output is reproduced.

    This isolates the only difference behind the strict xfails above: the
    quadrature of R's ``inprod()``, not the estimator.
    """
    from fabel.sparse import _fit_cov, _residuals

    inp, out = case["input"], case["output"]
    rtol = case_rtol(MODULE, case)
    cov_basis = build_basis(inp["cov_basis"])
    times, values = _data(case)
    mean = sparse_mean(values, times, build_basis(inp["mean_basis"]), lam=inp["mean_lambda"])
    cov = _fit_cov(
        times,
        _residuals(times, values, mean),
        cov_basis,
        inp["cov_lambda"],
        LDO(2),
        gram=np.asarray(out["inprod_cov_cov"]),
    )
    compare(cov, out["cov_coefs"], rtol)
    vals, harm = _harmonics(
        np.asarray(out["cov_coefs"]),
        cov_basis,
        build_basis(inp["harm_basis"]),
        inp["harm_lambda"],
        LDO(2),
        inp["nharm"],
        cross_gram=np.asarray(out["inprod_harm_cov"]),
    )
    compare(vals, out["values"], rtol)
    compare(harm, out["harmonics"], rtol)


def test_r_scores_use_one_point() -> None:
    """Measure the scoresPACE defect: each finite R score row is λ ∘ ξ(τ) times a scalar."""
    case = _cases()["pace_weather_fourier5_bspline7"]
    out = case["output"]
    harm_basis = build_basis(case["input"]["harm_basis"])
    days = np.arange(1.0, 366.0)
    lam_xi = harm_basis(days) @ np.asarray(out["harmonics"]) * np.asarray(out["values"])
    scores = np.array([[np.nan if v is None else v for v in row] for row in out["scores"]])
    finite = np.all(np.isfinite(scores), axis=1)
    assert 0 < int(finite.sum()) < scores.shape[0] // 4
    for row in scores[finite]:
        # 2-D cross product with every candidate point: zero means parallel.
        cross = lam_xi[:, 0] * row[1] - lam_xi[:, 1] * row[0]
        scale = np.abs(lam_xi).max(axis=1) * np.abs(row).max()
        assert float(np.min(np.abs(cross) / scale)) < 1e-12
