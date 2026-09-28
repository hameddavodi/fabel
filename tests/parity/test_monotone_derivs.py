"""Parity of monotone and positive smooth derivatives against R ``fda`` 6.3.0.

``tests/golden/monotone_derivs.json`` (from ``tools/golden_r/monotone_derivs.R``)
records ``eval.monfd`` (``h = ∫ exp W``, derivatives 0-3), ``eval.posfd``
(``exp W``, derivatives 0-2) and ``predict`` of a ``smooth.monotone`` fit
(``β₀ + β₁ h`` and ``β₁ Dᵏ h``) for given latent functions ``W``.  The
:class:`~fdatools.smoothing.SmoothResult` below is built from R's ``W`` and ``β``
directly, so the comparison tests evaluation only, not the iterative fit.

Two R outputs are measurably wrong and are strict xfails (with a passing test
that pins down what R computes instead):

* ``eval.monfd(t, W, 0)`` integrates ``exp W`` numerically.  Against an
  adaptive-quadrature reference (``scipy.integrate.quad``, tolerance 1e-14)
  fdatools' Gauss-Legendre value is off by at most 4.3e-14 of the curve's range,
  R's by 2.4e-6 (growth) and 1.8e-6 (fixed W); 9e-5 relative at a point.
* ``eval.posfd(t, W, 2)`` returns ``exp(W) D²W`` (to 1.5e-16), which misses the
  ``exp(W) (DW)²`` term of ``D² exp W``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy.integrate import quad

from fdatools import FData
from fdatools.smoothing import SmoothResult

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "monotone_derivs"
CASES = {case["name"]: case for case in golden_cases(MODULE)}


def _result(case: dict[str, Any], constraint: str, beta: Any) -> SmoothResult:
    inp = case["input"]
    coefs = np.asarray(inp["coefs"], dtype=float)
    fd = FData(coefs, build_basis(inp["basis"]))
    n_curves = coefs.shape[1]
    return SmoothResult(
        fd=fd,
        df=0.0,
        gcv=np.zeros(n_curves),
        sse=0.0,
        penalty_matrix=np.zeros((coefs.shape[0], coefs.shape[0])),
        lam=0.0,
        y2c_map=None,
        beta=None if beta is None else np.asarray(beta, dtype=float).reshape(2, -1),
        constraint=constraint,
    )


def _unit_beta(case: dict[str, Any]) -> np.ndarray:
    n_curves = np.asarray(case["input"]["coefs"]).shape[1]
    return np.vstack([np.zeros(n_curves), np.ones(n_curves)])


def test_every_golden_case_is_covered() -> None:
    assert set(CASES) == {"monfd_posfd_fixed_w_two_curves", "monfd_growth_girl1"}


_MONFD_QUADRATURE = (
    "R's eval.monfd(t, W, 0) integrates exp W numerically: measured against "
    "scipy.integrate.quad (tol 1e-14) R is off by 2.4e-6 (growth) / 1.8e-6 (fixed W) "
    "of the range, 9e-5 relative pointwise; fdatools by at most 4.3e-14"
)


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize(
    "deriv",
    [pytest.param(0, marks=pytest.mark.xfail(strict=True, reason=_MONFD_QUADRATURE)), 1, 2, 3],
)
def test_monotone_derivatives_match_eval_monfd(name: str, deriv: int) -> None:
    case = CASES[name]
    result = _result(case, "monotone", _unit_beta(case))
    t = np.asarray(case["input"]["t"])
    expected = np.asarray(case["output"]["monfd"][deriv], dtype=float)
    compare(np.asarray(result(t, deriv)).reshape(expected.shape), expected, case_rtol(MODULE, case))


_POSFD_SECOND = (
    "R's eval.posfd(t, W, 2) returns exp(W) D2W (to 1.5e-16), not D2 exp W = "
    "exp(W) (D2W + (DW)^2); measured 0.87 relative difference"
)


@pytest.mark.parametrize(
    "deriv", [0, 1, pytest.param(2, marks=pytest.mark.xfail(strict=True, reason=_POSFD_SECOND))]
)
def test_positive_derivatives_match_eval_posfd(deriv: int) -> None:
    case = CASES["monfd_posfd_fixed_w_two_curves"]
    result = _result(case, "positive", None)
    t = np.asarray(case["input"]["t"])
    expected = np.asarray(case["output"]["posfd"][deriv], dtype=float)
    compare(result(t, deriv), expected, case_rtol(MODULE, case))


@pytest.mark.parametrize(
    "deriv",
    [pytest.param(0, marks=pytest.mark.xfail(strict=True, reason=_MONFD_QUADRATURE)), 1, 2, 3],
)
def test_scaled_monotone_derivatives_match_predict(deriv: int) -> None:
    case = CASES["monfd_growth_girl1"]
    result = _result(case, "monotone", case["input"]["beta"])
    t = np.asarray(case["input"]["t"])
    expected = np.asarray(case["output"]["predict"][deriv], dtype=float)
    compare(np.asarray(result(t, deriv))[:, 0], expected, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", sorted(CASES))
def test_monotone_values_match_adaptive_quadrature(name: str) -> None:
    """``∫ₐᵗ exp W`` agrees with an independent adaptive quadrature."""
    case = CASES[name]
    result = _result(case, "monotone", _unit_beta(case))
    t = np.asarray(case["input"]["t"])[::4]
    latent = result.fd
    lower = latent.domain[0]
    n_curves = np.asarray(case["input"]["coefs"]).shape[1]

    def integrand(s: float, k: int) -> float:
        return float(np.exp(np.asarray(latent(np.array([s])))[0, k]))

    reference = np.array(
        [
            [
                quad(integrand, lower, x, args=(k,), epsabs=1e-13, epsrel=1e-13, limit=200)[0]
                for k in range(n_curves)
            ]
            for x in t
        ]
    )
    ours = np.asarray(result(t, 0)).reshape(reference.shape)
    compare(ours, reference, 1e-12)


def test_r_second_derivative_of_a_positive_function_drops_the_square_term() -> None:
    case = CASES["monfd_posfd_fixed_w_two_curves"]
    latent = _result(case, "positive", None).fd
    t = np.asarray(case["input"]["t"])
    r_value = np.asarray(case["output"]["posfd"][2], dtype=float)
    compare(np.exp(latent(t)) * latent(t, 2), r_value, 1e-14)
