"""Parity of :mod:`fabel.core` against golden output from R ``fda`` 6.3.0.

Cases whose output is a function (a basis plus coefficients) are compared by
evaluating both representations on a dense grid: the two sides may legitimately
choose different bases for the same function, and only the function matters.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import LDO, BiFData, BSpline, Constant, Exponential, FData, Fourier, Monomial, inprod
from fabel.basis import Basis

from .conftest import case_rtol, golden_cases

pytestmark = pytest.mark.parity

MODULE = "core"
GRID = 401

#: R's ``inprod`` stops refining its Romberg iteration at a relative change of
#: 1e-4 -- its own help text promises only "four to five significant digits" --
#: so its values carry an error far above this module's 1e-6 tolerance.  For a
#: pair of B-spline bases the golden file calls ``inprod.bspline`` instead,
#: which is exact; no such exact entry point exists for the other bases.
_INPROD_REASON = (
    "R's inprod integrates by Romberg iteration with a convergence tolerance of "
    "1e-4 and returns values ~1.3e-4 away from the exact integral.  Its public "
    "API exposes no tolerance argument and no exact alternative outside "
    "inprod.bspline, which does not accept a Fourier basis.  Fabel uses "
    "Gauss-Legendre panels that are exact for these integrands -- its Gram "
    "matrices match R's own eval.penalty to 1e-8."
)

#: Golden cases where R ``fda`` 6.3.0 returns a coarse approximation of a
#: function Fabel represents exactly.  The relative errors are measured against
#: the exact function reconstructed from the same golden inputs.
R_FDA_DEFECTS: dict[str, str] = {
    "fd_mul_bspline": (
        "times.fd keeps the factors' break points, so the order-7 product space "
        "is C^5 where the true product is only C^2 and cannot be represented.  "
        "R's answer is 12.8% away from the exact product (max 0.254 on a curve "
        "of size 1.99), and is not even the L2 projection onto its own basis.  "
        "Fabel raises the interior knot multiplicities and is exact."
    ),
    "fd_power2": (
        "^.fd projects onto an arbitrary uniform refinement (80 intervals, "
        "order 7) whose knots miss the curve's own breaks at k/7, leaving a "
        "5.7e-5 relative error.  Fabel squares the curve exactly."
    ),
    "inprod_fourier_L0_0": _INPROD_REASON,
    "inprod_fourier_L1_1": _INPROD_REASON,
    "inprod_fourier_L2_2": _INPROD_REASON,
    "inprod_basis_bspline_x_fourier": _INPROD_REASON,
    "deriv_fd_bspline_order6_L2": (
        "deriv.fd re-expands D^2 x in the *original* order-6 basis, which holds "
        "only C^4 functions while D^2 x is C^2; the result is 1.3% off (max "
        "59.9 on a curve of size 4655).  Fabel returns the exact order-4 spline."
    ),
}


def build_basis(spec: dict[str, Any]) -> Basis:
    """Construct the Fabel basis described by a golden basis specification."""
    kind = spec["type"]
    domain = (float(spec["rangeval"][0]), float(spec["rangeval"][1]))
    params = spec.get("params")
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


def build_operator(case_input: dict[str, Any], key: str = "lfd") -> LDO:
    """Build the operator a golden case asks for.

    A case records its operator either as a weight vector or as the R call that
    produced it, e.g. the string ``"int2Lfd(2)"``.
    """
    if "lfd_weights" in case_input:
        return LDO(weights=[float(w) for w in case_input["lfd_weights"]])
    value = case_input[key]
    if isinstance(value, str):
        return LDO(int(value.removeprefix("int2Lfd(").removesuffix(")")))
    return LDO(int(value))


def grid_for(basis: Basis) -> Any:
    """Return a dense evaluation grid over a basis domain."""
    lower, upper = basis.domain
    return np.linspace(lower, upper, GRID)


def compare(actual: Any, expected: Any, rtol: float) -> None:
    """Compare two arrays, scaling ``atol`` by the magnitude of the expected values."""
    want = np.asarray(expected, dtype=float)
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), want, rtol=rtol, atol=1e-12 * max(1.0, scale)
    )


def compare_functions(actual: FData, output: dict[str, Any], rtol: float) -> None:
    """Compare a computed :class:`FData` with R's, pointwise on a dense grid."""
    expected = FData(np.asarray(output["coefs"], dtype=float), build_basis(output["basis"]))
    grid = grid_for(expected.basis)
    compare(actual(grid), expected(grid), rtol)


def run_case(name: str, case: dict[str, Any], rtol: float) -> None:
    """Replay one golden core case."""
    inp, out = case["input"], case["output"]
    family = name.split("_")[0]

    if name == "inprod_basis_bspline_x_fourier":
        compare(inprod(build_basis(inp["basis1"]), build_basis(inp["basis2"])), out["values"], rtol)
        return

    if name == "bifd_construct_4d_coef":
        sbasis, tbasis = build_basis(inp["sbasis"]), build_basis(inp["tbasis"])
        coefs = np.asarray(inp["coef"], dtype=float)
        bifd = BiFData(coefs, sbasis, tbasis)
        expected = BiFData(
            np.asarray(out["coefs"], dtype=float),
            build_basis(out["sbasis"]),
            build_basis(out["tbasis"]),
        )
        s, t = grid_for(sbasis)[::40], grid_for(tbasis)[::40]
        compare(bifd(s, t), expected(s, t), rtol)
        return

    if name == "var_fd_bspline_n10curves":
        basis = build_basis(inp["basis"])
        fd = FData(np.asarray(inp["coefs"], dtype=float), basis)
        expected = BiFData(
            np.asarray(out["coefs"], dtype=float),
            build_basis(out["sbasis"]),
            build_basis(out["tbasis"]),
        )
        s = grid_for(basis)[::20]
        compare(fd.cov()(s, s), expected(s, s), rtol)
        return

    if family == "evalfd":
        basis = build_basis(inp["basis"])
        fd = FData(np.asarray(inp["coefs"], dtype=float), basis)
        t = np.asarray(inp["t"], dtype=float)
        operator: int | LDO = inp["deriv"] if "deriv" in inp else build_operator(inp)
        compare(fd(t, operator), out["values"], rtol)
        return

    if family == "inprod":
        basis = build_basis(inp["basis"])
        fd = FData(np.asarray(inp["coefs"], dtype=float), basis)
        compare(inprod(fd, fd, lfd1=inp["lfd1"], lfd2=inp["lfd2"]), out["values"], rtol)
        return

    basis = build_basis(inp["basis"])
    if name.startswith("fd_add"):
        left = FData(np.asarray(inp["coefs1"], dtype=float), basis)
        right = FData(np.asarray(inp["coefs2"], dtype=float), basis)
        compare_functions(left + right, out, rtol)
        return
    if name.startswith("fd_mul"):
        left = FData(np.asarray(inp["coefs1"], dtype=float), basis)
        right = FData(np.asarray(inp["coefs2"], dtype=float), basis)
        compare_functions(left * right, out, rtol)
        return

    fd = FData(np.asarray(inp["coefs"], dtype=float), basis)
    actions = {
        "mean": lambda: fd.mean(),
        "sd": lambda: fd.std(),
        "center": lambda: fd.center(),
        "deriv": lambda: fd.derivative(int(inp["lfd"])),
    }
    if family in actions:
        compare_functions(actions[family](), out, rtol)
        return
    if name == "fd_scalar_mul":
        compare_functions(fd * 2.0, out, rtol)
        return
    if name == "fd_power2":
        compare_functions(fd**2, out, rtol)
        return
    raise AssertionError(f"golden core case {name!r} is not covered by the parity harness")


CASES = golden_cases(MODULE)


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_core_matches_r_fda(case: dict[str, Any], request: pytest.FixtureRequest) -> None:
    """Replay one golden core case against R ``fda`` 6.3.0."""
    name = case["name"]
    if name in R_FDA_DEFECTS:
        request.applymarker(pytest.mark.xfail(strict=True, reason=R_FDA_DEFECTS[name]))
    run_case(name, case, case_rtol(MODULE, case))
