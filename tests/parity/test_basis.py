"""Parity of :mod:`fabel.basis` against golden output from R ``fda`` 6.3.0.

Every case in ``tests/golden/basis.json`` is replayed here.  A handful of cases
are marked ``xfail`` because R ``fda`` 6.3.0 is demonstrably wrong there; each
carries the mathematical reason in :data:`R_FDA_DEFECTS`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import LDO, Basis, BSpline, Constant, Exponential, Fourier, Monomial, Polygonal, Power

from .conftest import case_rtol, golden_cases

pytestmark = pytest.mark.parity

MODULE = "basis"

#: Golden cases where R ``fda`` 6.3.0 disagrees with the mathematics.  Fabel
#: implements the correct result, so parity is expected to fail.
R_FDA_DEFECTS: dict[str, str] = {
    "bspline_penalty_k4_n4_dom0_1_L0": (
        "With no interior knots R returns the monomial Gram (the Hilbert matrix) "
        "instead of the Bernstein-basis Gram; Fabel's [0, 0] entry is 1/7."
    ),
    "bspline_penalty_k4_n4_dom0_1_L1": (
        "Same zero-interior-knot path as L0: R's first rows are identically zero, "
        "which is impossible for a basis with nonzero first derivatives."
    ),
    "bspline_penalty_k4_n4_dom0_1_L2": (
        "R returns [[0,0,0,0],[0,0,0,0],[0,0,4,6],[0,0,6,12]].  D^2 B_0 = 6(1-t), "
        "so entry [0, 0] is the integral of 36 (1-t)^2 = 12, not 0."
    ),
    "fourier_periodmismatch_penalty_L0": (
        "Period 2 over domain (0, 1) leaves R on a numerical quadrature accurate "
        "only to ~1.4e-6 (its [7, 7] is 0.499999); the exact value is 1/2."
    ),
    "fourier_periodmismatch_penalty_L1": "Same partial-period quadrature error as L0.",
    "fourier_periodmismatch_penalty_L2": "Same partial-period quadrature error as L0.",
    "monomial_eval_n3_d2": (
        "R's second-derivative coefficient is prod_r (e - 2r) rather than "
        "prod_r (e - r), so D^2 t^e comes out as e (e - 2) t^(e-2)."
    ),
    "monomial_eval_n4_d2": "Same second-derivative coefficient defect as n3_d2.",
    "monomial_eval_n5_d2": "Same second-derivative coefficient defect as n3_d2.",
    "monomial_eval_n6_d2": "Same second-derivative coefficient defect as n3_d2.",
    "monomial_customexp_eval_d2": "Same second-derivative coefficient defect as n3_d2.",
    "power_cfg1_eval_d1": (
        "R zeroes the derivative of every power except e = 2, and returns 4t there instead of 2t."
    ),
    "power_cfg2_eval_d1": "Same first-derivative defect as cfg1.",
    "power_cfg2_penalty_L0": (
        "Where e_i + e_j = -1 the antiderivative is a logarithm; R evaluates 0/0 "
        "and returns NaN.  The exact entry is log(3 / 0.5)."
    ),
    "polygonal_n5_eval_d1": (
        "R's eval.basis ignores nderiv for a polygonal basis and returns the "
        "values themselves; its own eval.penalty uses the true derivatives."
    ),
    "polygonal_n11_eval_d1": "Same ignored-nderiv defect as n5.",
    "polygonal_n20_eval_d1": "Same ignored-nderiv defect as n5.",
}


def build_basis(name: str, inp: dict[str, Any]) -> Basis:
    """Construct the Fabel basis described by one golden case."""
    family = name.split("_")[0]
    domain: Any = tuple(inp["domain"]) if "domain" in inp else None
    if family == "bspline":
        return BSpline(
            domain=domain,
            n_basis=inp.get("n_basis"),
            order=inp["order"],
            breaks=inp.get("breaks"),
        )
    if family == "fourier":
        return Fourier(domain=domain, n_basis=inp["n_basis"], period=inp.get("period"))
    if family == "monomial":
        return Monomial(
            domain=domain, n_basis=inp.get("n_basis", 2), exponents=inp.get("exponents")
        )
    if family == "exponential":
        return Exponential(domain=domain, rates=inp["ratevec"])
    if family == "power":
        return Power(domain=domain, exponents=inp["exponents"])
    if family == "constant":
        return Constant(domain=domain)
    if family == "polygonal":
        return Polygonal(argvals=inp["argvals"])
    raise AssertionError(f"unknown golden basis family {family!r}")


def compare(actual: Any, expected: Any, rtol: float) -> None:
    """Compare against R, scaling ``atol`` by the magnitude of the expected values.

    A bare ``atol=1e-12`` on a matrix whose entries reach 1e5 would demand
    agreement to 5e-18 relative, far below double precision.  Scaling keeps the
    floor at roughly 4500 machine epsilons of the matrix norm.
    """
    want = np.asarray(expected, dtype=float)
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float),
        want,
        rtol=rtol,
        atol=1e-12 * max(1.0, scale),
    )


CASES = golden_cases(MODULE)


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_basis_matches_r_fda(case: dict[str, Any], request: pytest.FixtureRequest) -> None:
    """Replay one golden basis case and compare evaluation or penalty output."""
    name = case["name"]
    if name in R_FDA_DEFECTS:
        request.applymarker(pytest.mark.xfail(strict=True, reason=R_FDA_DEFECTS[name]))
    inp, out = case["input"], case["output"]
    basis = build_basis(name, inp)
    assert basis.n_basis == out["nbasis"]

    rtol = case_rtol(MODULE, case)
    if "values" in out:
        compare(basis(np.asarray(inp["t"], dtype=float), deriv=inp["deriv"]), out["values"], rtol)
    else:
        compare(basis.penalty(LDO(inp["lfd"])), out["penalty"], rtol)
