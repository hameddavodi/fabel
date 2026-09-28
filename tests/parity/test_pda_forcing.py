"""Parity of forced PDA and of its stability analysis against R ``fda`` 6.3.0.

``tests/golden/pda_forcing.json`` (from ``tools/golden_r/pda_forcing.R``)
records ``pda.fd`` with forcing functions (``awtlist``/``ufdlist``) and
``eigen.pda`` on the fitted equations.  Every case stores the exact curves and
forcing functions R used, so the fits below start from R's inputs.

R's ``eigen.pda`` orders each row of eigenvalues by decreasing modulus, but
eigenvalues of (nearly) equal modulus may come out of LAPACK in either order,
so rows are compared after sorting by real, then imaginary part.
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel import Constant, FData
from fabel.dynamics import PDA

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "pda_forcing"
CASES = {case["name"]: case for case in golden_cases(MODULE)}

FORCED_SINGLE = (
    "pda_forced_order1_constant_exp",
    "pda_forced_multicurve_bspline_weights",
    "pda_forced_refinery_constant",
    "pda_forced_refinery_two_forcings",
    "pda_forced_order2_constant",
)
SYSTEMS = ("pda_forced_system_order1", "eigen_pda_handwriting_system")
WITH_EIGEN = (
    "pda_forced_refinery_constant",
    "pda_forced_refinery_two_forcings",
    "pda_forced_order2_constant",
    "pda_forced_system_order1",
    "eigen_pda_lip_bspline11",
    "eigen_pda_handwriting_system",
)


def _fd(info: dict[str, Any]) -> FData:
    coefs = np.asarray(info["coefs"], dtype=float)
    return FData(coefs, build_basis(info["basis"]))


def _forcing(case: dict[str, Any]) -> list[FData]:
    return [_fd(info) for info in case["input"].get("u", [])]


@cache
def _fit(name: str) -> PDA:
    case = CASES[name]
    inp = case["input"]
    if name == "pda_forced_order1_constant_exp":
        return PDA(order=1).fit(_fd(inp["x"]), forcing=_forcing(case))
    if name == "pda_forced_multicurve_bspline_weights":
        wbasis = build_basis(inp["wbasis"])
        return PDA(
            order=1,
            weight_basis=wbasis,
            lam=inp["bwt_lambda"],
            penalty=inp["penalty"],
            forcing_basis=wbasis,
            forcing_lam=inp["awt_lambda"],
        ).fit(_fd(inp["x"]), forcing=_forcing(case))
    if name == "pda_forced_refinery_constant":
        return PDA(order=1).fit(_fd(inp["x"]), forcing=_forcing(case))
    if name == "pda_forced_refinery_two_forcings":
        wbasis = build_basis(inp["wbasis"])
        return PDA(
            order=1,
            weight_basis=wbasis,
            lam=inp["bwt_lambda"],
            penalty=inp["penalty"],
            forcing_basis=[wbasis, build_basis(inp["cbasis"])],
            forcing_lam=[inp["awt_lambda"], 0.0],
        ).fit(_fd(inp["x"]), forcing=_forcing(case))
    if name == "pda_forced_order2_constant":
        return PDA(order=2).fit(_fd(inp["x"]), forcing=_forcing(case))
    if name == "pda_forced_system_order1":
        u = _fd(inp["u"][0])
        return PDA(order=1).fit(_system_curves(case), forcing=[u, u])
    if name == "eigen_pda_lip_bspline11":
        return PDA(
            order=2, weight_basis=build_basis(inp["wbasis"]), lam=inp["lambda"], penalty=2
        ).fit(_fd(inp["x"]))
    return PDA(
        order=2,
        weight_basis=build_basis(inp["wbasis"]),
        lam=inp["lambda"],
        penalty=inp["penalty"],
        n_grid=inp["nfine"],
    ).fit(_system_curves(case))


def _system_curves(case: dict[str, Any]) -> FData:
    parts = [_fd(info) for info in case["input"]["x"]]
    coefs = np.stack([np.asarray(p.coefs).reshape(p.basis.n_basis, -1) for p in parts], axis=2)
    return FData(coefs, parts[0].basis)


# --------------------------------------------------------------------------- #
# forced single equations
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", FORCED_SINGLE)
def test_weight_coefficients(name: str) -> None:
    case = CASES[name]
    pda = _fit(name)
    for weight, expected in zip(pda.weights_, case["output"]["bwt_coefs"], strict=True):
        compare(weight.coefs, expected, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", FORCED_SINGLE)
def test_forcing_weight_coefficients(name: str) -> None:
    case = CASES[name]
    pda = _fit(name)
    expected_all = case["output"]["awt_coefs"]
    assert len(pda.forcing_weights_) == len(expected_all)
    for weight, expected in zip(pda.forcing_weights_, expected_all, strict=True):
        compare(weight.coefs, expected, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", FORCED_SINGLE)
def test_residual_coefficients(name: str) -> None:
    case = CASES[name]
    compare(_fit(name).residuals_.coefs, case["output"]["resfd_coefs"][0], case_rtol(MODULE, case))


@pytest.mark.parametrize("name", FORCED_SINGLE)
def test_transform_reproduces_the_residuals(name: str) -> None:
    case = CASES[name]
    residuals = _fit(name).transform(_fd(case["input"]["x"]), forcing=_forcing(case))
    compare(residuals.coefs, case["output"]["resfd_coefs"][0], case_rtol(MODULE, case))


def test_weight_values_on_a_grid() -> None:
    case = CASES["pda_forced_multicurve_bspline_weights"]
    out = case["output"]
    grid = np.asarray(out["grid"])
    pda = _fit(case["name"])
    compare(pda.weights_[0](grid)[:, 0], out["bwt_vals"], case_rtol(MODULE, case))
    compare(pda.forcing_weights_[0](grid)[:, 0], out["awt_vals"], case_rtol(MODULE, case))


# --------------------------------------------------------------------------- #
# systems
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", SYSTEMS)
def test_system_weights(name: str) -> None:
    case = CASES[name]
    pda = _fit(name)
    fitted = [
        weight.coefs for equation in pda.weights_ for variable in equation for weight in variable
    ]
    expected = case["output"]["bwt_coefs_flat"]
    for ours, theirs in zip(fitted, expected, strict=True):
        compare(ours, theirs, case_rtol(MODULE, case))


def test_system_forcing_weights() -> None:
    case = CASES["pda_forced_system_order1"]
    fitted = [w.coefs for equation in _fit(case["name"]).forcing_weights_ for w in equation]
    expected = case["output"]["awt_coefs_flat"]
    assert len(fitted) == len(expected) == 2
    for ours, theirs in zip(fitted, expected, strict=True):
        compare(ours, theirs, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", SYSTEMS)
def test_system_residuals(name: str) -> None:
    case = CASES[name]
    residuals = np.asarray(_fit(name).residuals_.coefs)
    for i, expected in enumerate(case["output"]["resfd_coefs"]):
        want = np.asarray(expected, dtype=float).reshape(residuals.shape[0], -1)
        compare(residuals[:, :, i], want, case_rtol(MODULE, case))


# --------------------------------------------------------------------------- #
# stability (eigen.pda)
# --------------------------------------------------------------------------- #


def _stability(name: str) -> Any:
    """Run ``PDA.stability`` on R's ``argvals``, with the forcing R was given.

    The second-order case was fitted with three copies of ``u = 1`` (one per
    curve); ``eigen.pda`` got a single ``u = 1``, and so does Fabel.
    """
    grid = np.asarray(CASES[name]["output"]["argvals"])
    pda = _fit(name)
    if name == "pda_forced_order2_constant":
        one = FData(np.array([1.0]), Constant(domain=pda.domain_))
        return pda.stability(n_points=grid.size, forcing=one)
    return pda.stability(n_points=grid.size)


def _canonical(values: np.ndarray) -> np.ndarray:
    order = np.lexsort((values.imag, values.real), axis=1)
    return np.take_along_axis(values, order, axis=1)


@pytest.mark.parametrize("name", WITH_EIGEN)
def test_stability_eigenvalues(name: str) -> None:
    case = CASES[name]
    out = case["output"]
    grid = np.asarray(out["argvals"])
    result = _stability(name)
    compare(result.t, grid, 1e-14)
    expected = np.asarray(out["eig_re"]) + 1j * np.asarray(out["eig_im"])
    expected = expected.reshape(grid.size, -1)
    ours = _canonical(np.asarray(result.eigenvalues))
    theirs = _canonical(expected)
    compare(ours.real, theirs.real, case_rtol(MODULE, case))
    compare(ours.imag, theirs.imag, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", ["eigen_pda_lip_bspline11", "eigen_pda_handwriting_system"])
def test_unforced_limits_are_zero(name: str) -> None:
    case = CASES[name]
    limits = np.asarray(_stability(name).limits)
    np.testing.assert_array_equal(limits, np.zeros_like(limits))
    np.testing.assert_array_equal(np.asarray(case["output"]["limvals"]), 0.0)


# --------------------------------------------------------------------------- #
# equilibria (eigen.pda's limvals)
# --------------------------------------------------------------------------- #

FORCED_LIMITS = (
    "pda_forced_refinery_constant",
    "pda_forced_refinery_two_forcings",
    "pda_forced_order2_constant",
    "pda_forced_system_order1",
)


def _r_limits(name: str, shape: tuple[int, ...]) -> np.ndarray:
    return np.asarray(CASES[name]["output"]["limvals"], dtype=float).reshape(shape)


@pytest.mark.parametrize(
    "name", ["pda_forced_refinery_constant", "pda_forced_refinery_two_forcings"]
)
def test_first_order_limits_are_r_limvals_with_the_sign_flipped(name: str) -> None:
    """For one first-order equation R's ``limvals`` is ``-z*``.

    ``Dx = -β x + a u`` has the equilibrium ``z* = a u / β`` (the solution is
    drawn to it when ``β > 0``); R reports ``-a u / β``.  Measured: Fabel's
    limits equal minus R's to 4e-16 (constant weights) and 1.5e-11 relative
    (B-spline weights, R's eval.fd rounding).
    """
    case = CASES[name]
    ours = np.asarray(_stability(name).limits)
    compare(ours, -_r_limits(name, ours.shape), case_rtol(MODULE, case))


@pytest.mark.xfail(
    strict=True,
    reason=(
        "R 6.3.0 eigen.pda limvals is not the equilibrium -A(t)^-1 f(t): measured on "
        "refinery (one first-order equation) it is exactly -z* (sign flipped); for the "
        "second-order equation it is (0, -a/b0) where z* = (a/b0, 0); for the forced "
        "two-equation system it is (-0.5000000338, -1.9e-8) where z* = (0.4999999929, "
        "-0.3861262521) -- the second equation's forcing is lost. Fabel returns z*, "
        "checked against the ODE solution in tests/unit/test_dynamics.py."
    ),
)
@pytest.mark.parametrize("name", FORCED_LIMITS)
def test_limits_match_r_limvals(name: str) -> None:
    case = CASES[name]
    ours = np.asarray(_stability(name).limits)
    compare(ours, _r_limits(name, ours.shape), case_rtol(MODULE, case))


def test_every_golden_case_is_covered() -> None:
    assert set(CASES) == {*FORCED_SINGLE, *SYSTEMS, *WITH_EIGEN}


def test_lip_weights_and_residuals() -> None:
    case = CASES["eigen_pda_lip_bspline11"]
    pda = _fit(case["name"])
    for weight, expected in zip(pda.weights_, case["output"]["bwt_coefs"], strict=True):
        compare(weight.coefs, expected, case_rtol(MODULE, case))
    compare(pda.residuals_.coefs, case["output"]["resfd_coefs"][0], case_rtol(MODULE, case))
