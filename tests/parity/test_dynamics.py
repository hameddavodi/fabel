"""Parity of :mod:`fdatools.dynamics` against golden output from R ``fda`` 6.3.0.

``tests/golden/dynamics.json`` records ``pda.fd``'s weight functions and
residual functions for three single-equation problems.  The curves ``pda.fd``
was run on are read from ``tests/parity/data/dynamics_inputs.json`` (written by
``tools/export_dynamics_inputs.R``), not re-smoothed here: the lip fit is so
ill-conditioned (normal matrix condition number ~3.7e9) that two correct
double-precision solvers disagree in its 7th significant digit, which would
test the smoother, not the principal differential analysis.

Every case runs with :class:`~fdatools.dynamics.PDA`'s default ``n_grid=501``,
the trapezoidal discretisation ``pda.fd`` uses.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from fdatools import FData
from fdatools.dynamics import PDA
from fdatools.smoothing import smooth

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "dynamics"
INPUTS = Path(__file__).resolve().parent / "data" / "dynamics_inputs.json"
CASES = {case["name"]: case for case in golden_cases(MODULE)}

#: Which fixture curves each golden case was fitted to.
_CURVES = {
    "pda_fd_order1_analytic_exp_decay": "exp_decay_fd_coefs",
    "pda_fd_order2_lip_constant": "lip_fd_coefs",
    "pda_fd_order2_lip_bspline11": "lip_fd_coefs",
}


@cache
def _inputs() -> dict[str, Any]:
    with INPUTS.open(encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    return data


def _curves(name: str) -> FData:
    case = CASES[name]
    return FData(
        np.asarray(_inputs()[_CURVES[name]], dtype=float), build_basis(case["input"]["xbasis"])
    )


@cache
def _fit(name: str) -> PDA:
    case = CASES[name]
    curves = _curves(name)
    if name == "pda_fd_order1_analytic_exp_decay":
        return PDA(order=1).fit(curves)
    if name == "pda_fd_order2_lip_constant":
        return PDA(order=2).fit(curves)
    weight_basis = build_basis(case["input"]["wbasis"])
    return PDA(order=2, weight_basis=weight_basis, lam=1e-6, penalty=2).fit(curves)


def test_every_golden_case_is_covered() -> None:
    assert set(CASES) == set(_CURVES)


def test_input_curves_reproduce_the_golden_design() -> None:
    """The fixture curves are the golden raw data smoothed on the golden basis.

    This guards the fixture, not a parity claim: fdatools' own ``smooth`` agrees
    with R's ``smooth.basis`` to 7e-15 on the well-conditioned exponential fit
    and to 5.9e-7 on the lip fit, whose normal matrix has condition number
    3.7e9 (a 40-digit solve of fdatools' system lands 7e-8 from fdatools' result).
    """
    case = CASES["pda_fd_order1_analytic_exp_decay"]
    inp = case["input"]
    ours = smooth(
        np.asarray(inp["xvec"]), np.asarray(inp["tvec"]), basis=build_basis(inp["xbasis"]), lam=0.0
    ).fd
    compare(ours.coefs, _curves(case["name"]).coefs, 1e-12)
    lip = CASES["pda_fd_order2_lip_constant"]["input"]
    ours = smooth(
        np.asarray(lip["lip"]),
        np.asarray(lip["liptime"]),
        basis=build_basis(lip["xbasis"]),
        lam=1e-8,
        penalty=4,
    ).fd
    expected = np.asarray(_curves("pda_fd_order2_lip_constant").coefs)
    gap = np.max(np.abs(np.asarray(ours.coefs) - expected)) / np.max(np.abs(expected))
    assert gap < 1e-6


# --------------------------------------------------------------------------- #
# weight functions
# --------------------------------------------------------------------------- #

_WEIGHT_FIELDS = [
    ("pda_fd_order1_analytic_exp_decay", 0, "bwt_coefs", None),
    ("pda_fd_order2_lip_constant", 0, "bwt0_coefs", "bwt0_vals"),
    ("pda_fd_order2_lip_constant", 1, "bwt1_coefs", "bwt1_vals"),
    ("pda_fd_order2_lip_bspline11", 0, "bwt0_coefs", "bwt0_vals"),
    ("pda_fd_order2_lip_bspline11", 1, "bwt1_coefs", "bwt1_vals"),
]


@pytest.mark.parametrize(("name", "j", "coef_field", "vals_field"), _WEIGHT_FIELDS)
def test_weight_coefficients(name: str, j: int, coef_field: str, vals_field: str | None) -> None:
    case = CASES[name]
    weight = _fit(name).weights_[j]
    compare(weight.coefs, case["output"][coef_field], case_rtol(MODULE, case))


@pytest.mark.parametrize(("name", "j", "coef_field", "vals_field"), _WEIGHT_FIELDS)
def test_weight_values(name: str, j: int, coef_field: str, vals_field: str | None) -> None:
    case = CASES[name]
    output = case["output"]
    grid = np.asarray(output["bwt_grid"])
    expected = output["bwt_vals"] if vals_field is None else output[vals_field]
    compare(_fit(name).weights_[j](grid)[:, 0], expected, case_rtol(MODULE, case))


# --------------------------------------------------------------------------- #
# residual functions
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", sorted(_CURVES))
def test_residual_coefficients(name: str) -> None:
    case = CASES[name]
    compare(_fit(name).residuals_.coefs, case["output"]["resfd_coefs"], case_rtol(MODULE, case))


@pytest.mark.parametrize("name", sorted(_CURVES))
def test_transform_reproduces_the_residuals(name: str) -> None:
    case = CASES[name]
    residuals = _fit(name).transform(_curves(name))
    compare(residuals.coefs, case["output"]["resfd_coefs"], case_rtol(MODULE, case))


# --------------------------------------------------------------------------- #
# the fitted equation, solved
# --------------------------------------------------------------------------- #


def test_solution_matches_the_r_estimate() -> None:
    """``solve`` integrates ``Dx = -β x`` for R's estimated constant ``β``."""
    case = CASES["pda_fd_order1_analytic_exp_decay"]
    grid = np.asarray(case["output"]["bwt_grid"])
    beta = float(np.ravel(case["output"]["bwt_coefs"])[0])
    solution = _fit(case["name"]).solve(grid, [1.0])
    compare(solution, np.exp(-beta * grid), case_rtol(MODULE, case))


def test_solution_tracks_the_analytic_ground_truth() -> None:
    """The solution differs from ``exp(-4t)`` only by the estimation error of ``β``.

    ``|exp(-β̂ t) - exp(-β t)| / exp(-β t) = |expm1((β - β̂) t)|``, which is the
    whole of the admissible gap: the estimate itself (R's and fdatools' agree to
    7e-16) is 7.8e-8 above ``β = 4`` because the curve is an order-5 spline fit
    of ``exp(-4t)``, not ``exp(-4t)`` itself.
    """
    case = CASES["pda_fd_order1_analytic_exp_decay"]
    grid = np.asarray(case["output"]["bwt_grid"])
    truth = np.asarray(case["output"]["analytic_solution"])
    beta_true = float(case["input"]["beta_true"])
    beta_hat = float(_fit(case["name"]).weights_[0].coefs[0, 0])
    np.testing.assert_allclose(truth, np.exp(-beta_true * grid), rtol=1e-15)
    gap = np.abs(_fit(case["name"]).solve(grid, [1.0]) / truth - 1.0)
    bound = np.abs(np.expm1((beta_true - beta_hat) * grid)) + case_rtol(MODULE, case)
    assert np.all(gap <= bound)
