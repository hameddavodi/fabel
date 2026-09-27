"""Parity of :mod:`fabel.registration` against golden output from R ``fda`` 6.3.0.

Cases are parametrised by ``(case, field)`` so that a defect in one R output
(the inverse warp of ``landmarkreg``, say) does not hide agreement in the others.

Two kinds of comparison are made:

* the golden call itself (``register`` against ``register.fd`` /
  ``landmarkreg``, ``decompose`` against ``AmpPhaseDecomp``);
* R's post-processing in isolation: ``register(..., init=<R's Wfd>,
  max_iter=0)`` maps R's own latent functions to warps and registered curves,
  which separates "same optimum?" from "same warp for the same ``W``?".
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel.core import FData
from fabel.registration import RegistrationResult, register

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "registration"
CASES: dict[str, dict[str, Any]] = {case["name"]: case for case in golden_cases(MODULE)}

LANDMARK = "landmarkreg_growth_hgtf_pubertal_spurt"
GROWTH = "register_fd_growth_hgtf_to_mean"
AMP_PHASE = "ampphasedecomp_growth_hgtf_to_mean"
WEATHER = "register_fd_weather_periodic_crit1"

#: ``tools/golden_r/registration.R`` builds every ``WfdPar`` with ``Lfdobj = 2``
#: and ``lambda = 1``; the golden ``input`` block records only the bases.
WARP_LAMBDA = 1.0
WARP_PENALTY = 2

# --------------------------------------------------------------------------- #
# measured reasons for the fields R gets wrong
# --------------------------------------------------------------------------- #

_REGISTER_STOPS_SHORT = (
    "register.fd returns a point that is not a minimum of its own criterion. "
    "Evaluated with R's own discretisation (grid mean over 351 points, "
    "crit=2, lambda=1), R's Wfd scores 0.621731 / 0.304731 / 0.912346 / "
    "0.114607 / 0.310316 / 0.070639 / 0.815999 / 0.194697 / 0.282611 / "
    "0.646418 for the ten girls, with gradients of max-norm 0.021 to 0.371 "
    "there; Fabel's Newton iterate is stationary (gradient < 1e-10) at "
    "0.619192 / 0.304554 / 0.911651 / 0.114439 / 0.310197 / 0.070245 / "
    "0.808155 / 0.194462 / 0.282337 / 0.639359 -- lower on every curve.  R's "
    "own dbglev=2 trace agrees: it stops girls 1-4 with gradient lengths "
    "0.0516 / 0.1017 / 0.028 / 0.0081 after 'Criterion increased, terminating "
    "iterations' and 'Reset twice, terminating'.  The line search fails because "
    "R's criterion is not a function of the warp it returns: at a constant W "
    "(h(t) = t exactly) with x0 = x = 1000 t on [1, 18] it reports 1.1564 "
    "instead of 0 and a gradient of 3026; the moments of its criterion put its "
    "internal warp at about (1 - 1.0e-4) t there.  The optima differ by up to 0.089 in W."
)

_PERIODIC_STOPS_SHORT = (
    "register.fd returns a point that is not a minimum of its own criterion. "
    "Evaluated with R's discretisation (grid mean over 651 points, crit=1, "
    "lambda=1, periodic shift), R's (Wfd, shift) scores 0.634340 / 0.286287 / "
    "0.101681 / 1.122814 / 0.806865 for the five stations, with gradients of "
    "max-norm 0.0088 to 0.205 there; Fabel's Newton iterate is stationary "
    "(gradient < 1e-10) at 0.631612 / 0.285765 / 0.101659 / 1.120003 / "
    "0.798741 -- lower on every curve.  Shifts and warps trade off along a "
    "shallow valley, so the stationary point is far from R's: shifts 6.269 / "
    "-4.475 / 2.127 / 8.781 / -16.230 days against R's 6.820 / -4.428 / "
    "2.104 / 9.895 / -18.618, and W up to 4.1 apart.  R's internal warp is "
    "also not the one it returns (see the growth case)."
)

_R_INTERNAL_WARP = (
    "R's regfd is not the registered curve of the warp R returns.  Given R's "
    "own Wfd, x(h(t)) projected on the curve basis differs from R's regfd by "
    "3.6e-4 (growth) / 3.2e-4 (weather) relative.  With x(t) = t on a basis "
    "that also holds the warp, R's regfd differs from R's warpfd by 9e-5 "
    "although both are the same function x(h(t)) = h(t); and R's criterion at "
    "a constant W (h(t) = t exactly) is not zero -- its internal warp is about "
    "(1 - 1.0e-4) t there, not t."
)

_TRAPEZOID_WARP = (
    "R integrates exp(W) for warpfd with the trapezoidal rule on 1025 points "
    "(monfn) -- reproduced to 2e-12.  For the weather warps that rule is off "
    "the exact integral by up to 4.7e-4 days, which moves the least-squares "
    "warp coefficients by up to 8.1e-4; the smallest coefficient (-0.1143) is "
    "then 2.5e-3 out relative.  Fabel integrates exp(W) by Gauss-Legendre to "
    "rounding error."
)

_LANDMARK_WARP_QUADRATURE = (
    "landmarkreg's warpfd is not the exact warp of its own Wfd.  Integrating "
    "exp(W) for R's recorded Wfd to rounding error gives warps that differ "
    "from R's warpfd by up to 3.5e-5 years; the order-6 spline coefficients "
    "amplify that to 7.6e-5, so 6 of 350 coefficients miss 1e-5 relative "
    "(worst 6.9e-5).  No 1025-point trapezoid or grid choice reproduces R's "
    "values, so the gap is R's quadrature, not a different warp: Fabel's own "
    "Wfd agrees with R's to 5.9e-7."
)

_LANDMARK_INVERSE = (
    "landmarkreg's warpinvfd is not the inverse of its warpfd: "
    "warpinvfd(warpfd(t)) - t reaches 1.73 years (girl 3, t = 2.48), and "
    "warpinvfd rises from 1 to 3.53 over the first year for girl 1, where the "
    "true inverse reaches 2.17.  Fabel inverts its warps by Newton's method; "
    "its inverse composes with the warp to the identity within 2e-8."
)

_LANDMARK_REGISTERED = (
    "landmarkreg's regfd is not x(h(t)): for girl 3 it has a coefficient of "
    "3360.37 and a value of 421.27 cm at age 1.2 (the curve spans 67.6-183.2 "
    "cm), where x(h(1.2)) = 80.25 cm.  It is sampled through R's faulty inverse "
    "warp (above), whose steep start leaves the first years almost unsampled. "
    "Fabel projects x(h(t)) on the registration grid."
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    (LANDMARK, "warpfd_coefs"): _LANDMARK_WARP_QUADRATURE,
    (LANDMARK, "warpinvfd_coefs"): _LANDMARK_INVERSE,
    (LANDMARK, "regfd_coefs"): _LANDMARK_REGISTERED,
    (GROWTH, "regfd_coefs"): _REGISTER_STOPS_SHORT,
    (GROWTH, "warpfd_coefs"): _REGISTER_STOPS_SHORT,
    (GROWTH, "Wfd_coefs"): _REGISTER_STOPS_SHORT,
    (WEATHER, "regfd_coefs"): _PERIODIC_STOPS_SHORT,
    (WEATHER, "warpfd_coefs"): _PERIODIC_STOPS_SHORT,
    (WEATHER, "Wfd_coefs"): _PERIODIC_STOPS_SHORT,
    (WEATHER, "shift"): _PERIODIC_STOPS_SHORT,
}

#: Post-processing comparisons (R's Wfd in, R's warps out) that R gets wrong.
R_POSTPROCESS_DEFECTS: dict[tuple[str, str], str] = {
    (GROWTH, "regfd_coefs"): _R_INTERNAL_WARP,
    (WEATHER, "regfd_coefs"): _R_INTERNAL_WARP,
    (WEATHER, "warpfd_coefs"): _TRAPEZOID_WARP,
}


# --------------------------------------------------------------------------- #
# case inputs
# --------------------------------------------------------------------------- #


def _curves(case: dict[str, Any], key: str = "coefs") -> FData:
    """Return the curves recorded under ``input[key]`` in the case's curve basis."""
    return FData(np.asarray(case["input"][key], dtype=float), build_basis(case["input"]["basis"]))


def _golden_latent(case: dict[str, Any]) -> FData:
    """Return R's latent functions ``Wfd`` for a case."""
    basis = build_basis(case["input"]["wbasis"])
    return FData(np.asarray(case["output"]["Wfd_coefs"], dtype=float), basis)


def _continuous_options(case: dict[str, Any]) -> dict[str, Any]:
    """Translate the recorded ``register.fd`` arguments into ``register`` options."""
    inputs = case["input"]
    return {
        "warp_basis": build_basis(inputs["wbasis"]),
        "lam": WARP_LAMBDA,
        "penalty": WARP_PENALTY,
        "criterion": "eigen" if inputs["crit"] == 2 else "least_squares",
        "periodic": bool(inputs["periodic"]),
    }


@cache
def _fabel_result(name: str) -> RegistrationResult:
    """Run Fabel on the golden inputs of case ``name`` once per session."""
    case = CASES[name]
    if name == LANDMARK:
        inputs = case["input"]
        return register(
            _curves(case),
            landmarks=np.asarray(inputs["ximarks"], dtype=float)[:, None],
            target_landmarks=[inputs["x0marks"]],
            warp_basis=build_basis(inputs["wbasis"]),
            lam=WARP_LAMBDA,
            penalty=WARP_PENALTY,
        )
    return register(_curves(case), _curves(case, "y0fd_coefs"), **_continuous_options(case))


@cache
def _from_r_latent(name: str) -> RegistrationResult:
    """Apply R's own ``Wfd`` (and shift) without optimising."""
    case = CASES[name]
    return register(
        _curves(case),
        _curves(case, "y0fd_coefs"),
        init=_golden_latent(case),
        init_shift=np.asarray(case["output"]["shift"], dtype=float),
        max_iter=0,
        **_continuous_options(case),
    )


def _field(result: RegistrationResult, field: str) -> Any:
    """Return the Fabel counterpart of an R output field."""
    if field == "regfd_coefs":
        return result.registered.coefs
    if field == "warpfd_coefs":
        return result.warp.coefs
    if field == "Wfd_coefs":
        return result.latent.coefs
    if field == "shift":
        return result.shift
    if field == "warpinvfd_coefs":
        assert result.warp_inverse is not None
        return result.warp_inverse.coefs
    raise AssertionError(f"unknown field {field!r}")


def _mark(request: pytest.FixtureRequest, reason: str | None) -> None:
    if reason is not None:
        request.applymarker(pytest.mark.xfail(strict=True, reason=reason))


REGISTRATION_FIELDS = [
    (name, field) for name in (LANDMARK, GROWTH, WEATHER) for field in CASES[name]["output"]
]


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "field"), REGISTRATION_FIELDS)
def test_registration_matches_r(name: str, field: str, request: pytest.FixtureRequest) -> None:
    """``register`` reproduces ``landmarkreg`` / ``register.fd`` output fields."""
    _mark(request, R_FDA_DEFECTS.get((name, field)))
    case = CASES[name]
    compare(_field(_fabel_result(name), field), case["output"][field], case_rtol(MODULE, case))


@pytest.mark.parametrize(
    ("name", "field"),
    [(name, field) for name in (GROWTH, WEATHER) for field in ("warpfd_coefs", "regfd_coefs")],
)
def test_warps_from_r_latent_match_r(name: str, field: str, request: pytest.FixtureRequest) -> None:
    """R's own ``Wfd`` maps to R's ``warpfd`` / ``regfd`` through Fabel's warps."""
    _mark(request, R_POSTPROCESS_DEFECTS.get((name, field)))
    case = CASES[name]
    result = _from_r_latent(name)
    compare(_field(result, field), case["output"][field], case_rtol(MODULE, case))


def test_from_r_latent_keeps_r_latent_and_shift() -> None:
    """``max_iter=0`` leaves R's latent functions and shifts untouched."""
    for name in (GROWTH, WEATHER):
        case = CASES[name]
        result = _from_r_latent(name)
        compare(result.latent.coefs, case["output"]["Wfd_coefs"], 1e-14)
        compare(result.shift, case["output"]["shift"], 1e-14)


_AMP_PHASE_FIELDS = {"MS.amp": "amp_mse", "MS.pha": "phase_mse", "RSQR": "rsq", "C": "c"}


@pytest.mark.parametrize("field", list(_AMP_PHASE_FIELDS))
def test_decompose_matches_ampphasedecomp(field: str) -> None:
    """``decompose`` on R's own registration reproduces ``AmpPhaseDecomp``."""
    case = CASES[AMP_PHASE]
    source = CASES[case["input"]["source_case"]]
    curve_basis = build_basis(source["input"]["basis"])
    warp_basis = build_basis(source["input"]["wbasis"])
    output = source["output"]
    result = RegistrationResult(
        registered=FData(np.asarray(output["regfd_coefs"], dtype=float), curve_basis),
        warp=FData(np.asarray(output["warpfd_coefs"], dtype=float), warp_basis),
        unregistered=_curves(source),
        latent=_golden_latent(source),
        shift=np.asarray(output["shift"], dtype=float),
    )
    decomposition = result.decompose()
    compare(
        getattr(decomposition, _AMP_PHASE_FIELDS[field]),
        case["output"][field],
        case_rtol(MODULE, case),
    )
