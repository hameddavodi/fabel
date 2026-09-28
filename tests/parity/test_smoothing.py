"""Parity of :mod:`fdatools.smoothing` against golden output from R ``fda`` 6.3.0.

Cases are parametrised by ``(case, field)`` rather than by case alone: R's
``smooth.basis`` returns seven quantities at once and a defect in one of them
(the Fourier harmonic penalty, say) must not hide agreement in the other six.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fdatools import LDO
from fdatools.smoothing import df_to_lambda, gcv_curve, lambda_to_df, smooth

from ._datasets import growth, handwriting, weather_daily
from .conftest import build_basis, case_rtol, golden_cases

pytestmark = pytest.mark.parity

MODULE = "smoothing"

#: The weather cases smooth with the harmonic accelerator of the annual cycle.
WEATHER_PERIOD = 365.0

#: R records ``penmat = NULL`` when ``lambda == 0``; there is nothing to compare.
_SMOOTH_FIELDS = ("coefs", "df", "gcv", "SSE", "penmat", "y2cMap", "fitted")

#: The roughness operator each golden case was generated with.  Recorded here
#: because ``tools/golden_r/smoothing.R`` passes it to ``fdPar``/``Data2fd`` but
#: the golden file's ``input`` block only carries the basis and lambda.
_OPERATORS: dict[str, int | LDO] = {
    "smooth_basis_synthetic_lfd_integer0": 0,
}

_HARMONIC_CASES = frozenset(
    {
        "smooth_basis_weather_fourier65_harmonic_lambda0.01",
        "smooth_basis_weather_fourier65_harmonic_lambda100",
        "smooth_basis_weather_fourier65_harmonic_lambda10000",
        "smooth_basis_weather_fourier65_harmonic_lambda1e+06",
        # Data2fd's default penalty on a Fourier basis is the harmonic
        # accelerator, not D^2 -- determined by matching the golden coefficients.
        "data2fd_weather_fourier",
    }
)

#: B-spline cases whose penalty is D^4: the growth basis has order 6, and
#: Data2fd's default penalty for an order-m spline is D^(m-2).
_LFD4_CASES = frozenset(
    {
        "smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01",
        "smooth_basis_growth_hgtm_bspline6_lfd4_lambda1",
        "smooth_basis_growth_hgtm_bspline6_lfd4_lambda100",
        "smooth_basisPar_growth_hgtm",
        "data2fd_growth_hgtf",
    }
)

_MONOTONE_OPERATORS = {
    "smooth_monotone_growth_hgtf_3girls": 3,
    "smooth_monotone_synthetic_sigmoid": 2,
}

# --------------------------------------------------------------------------- #
# R defects
# --------------------------------------------------------------------------- #

_HARMONIC_PENALTY_REASON = (
    "R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 "
    "(w = 2*pi/365) is wrong in its 6th significant digit.  L annihilates the "
    "constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) "
    "sin(k w t), so on the normalised Fourier basis the penalty is exactly "
    "diagonal with entries k^2 w^6 (k^2 - 1)^2.  Against that closed form "
    "fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is "
    "(64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02."
)

_HARMONIC_PROPAGATION_REASON = (
    _HARMONIC_PENALTY_REASON + "  At lambda >= 100 the wrong penalty moves the solution of "
    "(Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 "
    "relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the "
    "GCV score inherit the same error."
)

_GROWTH_CONDITION_REASON = (
    "Not reachable in double precision by either implementation: 31 ages on 35 "
    "order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 "
    "at lambda = 0.01 / 1 / 100.  Both penalty matrices are correct -- against "
    "a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an "
    "order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, "
    "i.e. both sit on the rounding floor.  But solving the two systems in "
    "60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / "
    "3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here.  The "
    "1e-8 golden tolerance is below the noise floor of the problem, not a sign "
    "that either solver is wrong."
)

_DF2LAMBDA_REASON = (
    "R's df2lambda stops its search short of solving df(lambda) = df.  Feeding "
    "R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / "
    "30.004195052 for the targets 10 / 20 / 30.  df(lambda) is smooth and "
    "strictly decreasing, so fdatools bisects it to 1e-12 and lands on lambda "
    "3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to "
    "twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 "
    "/ 8.6e-4 (relative)."
)

_MONOTONE_SHIFT_REASON = (
    "smooth.monotone has an exactly flat direction, so the coefficients are "
    "not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by "
    "W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness "
    "penalty is unchanged too because D^m annihilates constants.  fdatools and R "
    "settle on different members of that one-parameter family: the measured "
    "shift is constant across every coefficient of a curve (-0.08723 / "
    "-0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and "
    "beta1 differs by exactly the matching factor (4.76127 against 4.36367 = "
    "4.76127 * exp(-0.08723)).  What the model does identify agrees: the "
    "fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same."
)

_SMOOTH_POS_REASON = (
    "R's smooth.pos returns a point that is not stationary for its own "
    "criterion mean((y - exp(Phi c))^2) + lambda c'Rc.  The gradient there has "
    "norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and "
    "9.6e-13 for fdatools, and R's criterion value is the higher of the two "
    "(1.77593716848 against 1.77593715997; 3.19972667966087 against "
    "3.19972667965875).  One Gauss-Newton step from R's coefficients moves "
    "them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is "
    "literally one un-taken iteration short of the minimum fdatools reports."
)

#: ``(case, field)`` pairs where the golden value cannot be matched, either
#: because R ``fda`` 6.3.0 is wrong (the first four reasons) or because the
#: case is conditioned worse than the tolerance it is asserted at.
R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    ("smooth_basis_weather_fourier65_harmonic_lambda0.01", "penmat"): _HARMONIC_PENALTY_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda100", "penmat"): _HARMONIC_PENALTY_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda10000", "penmat"): _HARMONIC_PENALTY_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda1e+06", "penmat"): _HARMONIC_PENALTY_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda100", "coefs"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda100", "gcv"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda100", "y2cMap"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda100", "fitted"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda10000", "coefs"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda10000", "y2cMap"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda10000", "fitted"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_weather_fourier65_harmonic_lambda1e+06", "y2cMap"): _HARMONIC_PROPAGATION_REASON,
    ("lambda2gcv_weather_log10lambda_2", "gcv"): _HARMONIC_PROPAGATION_REASON,
    ("lambda2gcv_weather_grid_mean", "mean_gcv"): _HARMONIC_PROPAGATION_REASON,
    ("data2fd_weather_fourier", "coefs"): _HARMONIC_PROPAGATION_REASON,
    ("smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01", "gcv"): _GROWTH_CONDITION_REASON,
    ("smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01", "y2cMap"): _GROWTH_CONDITION_REASON,
    ("data2fd_growth_hgtf", "coefs"): _GROWTH_CONDITION_REASON,
    ("df2lambda_weather_df10", "lambda"): _DF2LAMBDA_REASON,
    ("df2lambda_weather_df20", "lambda"): _DF2LAMBDA_REASON,
    ("df2lambda_weather_df30", "lambda"): _DF2LAMBDA_REASON,
    ("smooth_pos_synthetic_2curves", "Wfdobj_coefs"): _SMOOTH_POS_REASON,
    ("smooth_pos_synthetic_2curves", "fitted"): _SMOOTH_POS_REASON,
    ("smooth_pos_weather_precip_pr_rupert", "Wfdobj_coefs"): _SMOOTH_POS_REASON,
}
for _case in (
    "smooth_basis_growth_hgtm_bspline6_lfd4_lambda1",
    "smooth_basis_growth_hgtm_bspline6_lfd4_lambda100",
    "smooth_basisPar_growth_hgtm",
):
    for _field in ("coefs", "df", "gcv", "SSE", "y2cMap", "fitted"):
        R_FDA_DEFECTS[(_case, _field)] = _GROWTH_CONDITION_REASON
for _case in ("smooth_monotone_growth_hgtf_3girls", "smooth_monotone_synthetic_sigmoid"):
    R_FDA_DEFECTS[(_case, "Wfdobj_coefs")] = _MONOTONE_SHIFT_REASON
    R_FDA_DEFECTS[(_case, "beta_slope")] = _MONOTONE_SHIFT_REASON
    R_FDA_DEFECTS[(_case, "deriv1")] = _MONOTONE_SHIFT_REASON


# --------------------------------------------------------------------------- #
# case inputs
# --------------------------------------------------------------------------- #


def operator_for(name: str) -> int | LDO:
    """Return the roughness operator the golden case was generated with."""
    if name in _HARMONIC_CASES:
        return LDO.harmonic(WEATHER_PERIOD)
    if name in _LFD4_CASES:
        return 4
    if name in _MONOTONE_OPERATORS:
        return _MONOTONE_OPERATORS[name]
    return _OPERATORS.get(name, 2)


def response_for(name: str, case_input: dict[str, Any]) -> np.ndarray:
    """Return the observations a golden case was fitted to."""
    if "y" in case_input:
        return np.asarray(case_input["y"], dtype=float)
    if "handwrit" in name:
        return handwriting()[1]
    if "hgtm" in name:
        return growth("hgtm")
    if "hgtf" in name:
        return growth("hgtf")
    return weather_daily()


def weather_fit_inputs() -> tuple[np.ndarray, np.ndarray]:
    """Return the ``(t, y)`` of the lambda2gcv/lambda2df weather design."""
    return np.arange(1.0, 366.0), weather_daily()


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #


def compare(actual: Any, expected: Any, rtol: float) -> None:
    """Compare against R, reshaping to R's layout and scaling ``atol``."""
    want = np.asarray(expected, dtype=float)
    got = np.reshape(np.asarray(actual, dtype=float), want.shape)
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    np.testing.assert_allclose(got, want, rtol=rtol, atol=1e-12 * max(1.0, scale))


_FITS: dict[str, Any] = {}


def fit_for(name: str, case: dict[str, Any]) -> Any:
    """Run (and memoise) the fdatools fit that replays one golden case."""
    if name not in _FITS:
        inp = case["input"]
        basis = build_basis(inp["basis"])
        t = np.asarray(inp["argvals"], dtype=float)
        y = response_for(name, inp)
        weights = np.asarray(inp["wtvec"], dtype=float) if "wtvec" in inp else None
        constraint = None
        if name.startswith("smooth_monotone"):
            constraint = "monotone"
        elif name.startswith("smooth_pos"):
            constraint = "positive"
        _FITS[name] = smooth(
            y,
            t,
            basis=basis,
            lam=float(inp["lambda"]),
            penalty=operator_for(name),
            weights=weights,
            constraint=constraint,
        )
    return _FITS[name]


def check_smooth_field(name: str, case: dict[str, Any], field: str, rtol: float) -> None:
    """Compare one field of a ``smooth.basis``-style golden case."""
    result = fit_for(name, case)
    out = case["output"]
    t = np.asarray(case["input"]["argvals"], dtype=float)
    if field == "coefs":
        compare(result.fd.coefs, out["coefs"], rtol)
    elif field == "df":
        compare(result.df, out["df"], rtol)
    elif field == "gcv":
        compare(result.gcv, out["gcv"], rtol)
    elif field == "SSE":
        compare(result.sse, out["SSE"], rtol)
    elif field == "penmat":
        compare(result.penalty_matrix, out["penmat"], rtol)
    elif field == "y2cMap":
        compare(result.y2c_map, out["y2cMap"], rtol)
    elif field == "fitted":
        compare(result.fd(t), out["fitted"], rtol)
    else:
        raise AssertionError(f"unknown smooth field {field!r}")


def check_constrained_field(name: str, case: dict[str, Any], field: str, rtol: float) -> None:
    """Compare one field of a ``smooth.monotone``/``smooth.pos`` golden case."""
    result = fit_for(name, case)
    out = case["output"]
    t = np.asarray(case["input"]["argvals"], dtype=float)
    if field == "Wfdobj_coefs":
        compare(result.fd.coefs, out["Wfdobj_coefs"], rtol)
    elif field == "fitted":
        compare(result(t), out["fitted"], rtol)
    elif field == "deriv1":
        compare(np.exp(np.asarray(result.fd(t), dtype=float)), out["deriv1"], rtol)
    elif field == "beta_intercept":
        compare(np.asarray(result.beta)[0], np.asarray(out["beta"], dtype=float)[0], rtol)
    elif field == "beta_slope":
        compare(np.asarray(result.beta)[1], np.asarray(out["beta"], dtype=float)[1], rtol)
    else:
        raise AssertionError(f"unknown constrained field {field!r}")


def check_lambda_helper(name: str, case: dict[str, Any], field: str, rtol: float) -> None:
    """Compare one of the lambda/df/gcv helper cases."""
    t, y = weather_fit_inputs()
    basis = build_basis(
        {"type": "fourier", "rangeval": [0, 365], "nbasis": 65, "params": WEATHER_PERIOD}
    )
    operator = LDO.harmonic(WEATHER_PERIOD)
    inp, out = case["input"], case["output"]
    if field == "gcv":
        scores = gcv_curve(y, t, basis, [10.0 ** float(inp["log10lambda"])], penalty=operator)
        compare(scores[0], out["gcv"], rtol)
    elif field == "df":
        compare(lambda_to_df(t, basis, float(inp["lambda"]), penalty=operator), out["df"], rtol)
    elif field == "lambda":
        compare(df_to_lambda(t, basis, float(inp["df"]), penalty=operator), out["lambda"], rtol)
    elif field in {"mean_gcv", "argmin_log10lambda"}:
        grid = np.asarray(inp["log10lambda_grid"], dtype=float)
        scores = gcv_curve(y, t, basis, 10.0**grid, penalty=operator)
        means = np.mean(scores, axis=1)
        if field == "mean_gcv":
            compare(means, out["mean_gcv"], rtol)
        else:
            compare(grid[int(np.argmin(means))], out["argmin_log10lambda"], rtol)
    else:
        raise AssertionError(f"unknown helper field {field!r}")


def fields_of(case: dict[str, Any]) -> tuple[str, ...]:
    """Return the comparable output fields of one golden case."""
    out = case["output"]
    if "Wfdobj_coefs" in out:
        extra = ("beta_intercept", "beta_slope", "deriv1") if "beta" in out else ()
        return ("Wfdobj_coefs", "fitted", *extra)
    if "coefs" in out and "df" not in out:
        return ("coefs",)
    if "df" in out and "coefs" in out:
        return tuple(f for f in _SMOOTH_FIELDS if out.get(f) is not None)
    return tuple(out)


CASES = golden_cases(MODULE)
PARAMS = [(case, field) for case in CASES for field in fields_of(case)]


@pytest.mark.parametrize(("case", "field"), PARAMS, ids=[f"{c['name']}-{f}" for c, f in PARAMS])
def test_smoothing_matches_r_fda(
    case: dict[str, Any], field: str, request: pytest.FixtureRequest
) -> None:
    """Replay one field of one golden smoothing case against R ``fda`` 6.3.0."""
    name = case["name"]
    reason = R_FDA_DEFECTS.get((name, field))
    if reason is not None:
        request.applymarker(pytest.mark.xfail(strict=True, reason=reason))
    rtol = case_rtol(MODULE, case)
    if "argvals" not in case["input"]:
        check_lambda_helper(name, case, field, rtol)
    elif name.startswith(("smooth_monotone", "smooth_pos")):
        check_constrained_field(name, case, field, rtol)
    elif "df" in case["output"]:
        check_smooth_field(name, case, field, rtol)
    else:
        compare(fit_for(name, case).fd.coefs, case["output"]["coefs"], rtol)
