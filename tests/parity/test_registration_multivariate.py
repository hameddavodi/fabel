"""Parity of multivariate registration against golden output from R ``fda`` 6.3.0.

The golden cases (``tools/golden_r/registration_multivariate.R``) register
curves with two variables: the hip and knee angles of the gait data and the
``x`` and ``y`` pen positions of the handwriting data.

R's ``register.fd`` fits the warps of multivariate curves to the first variable
alone (every case records R's univariate first-variable fit too, and the two
agree to rounding), so Fabel is run with ``var_weights=[1, 0]`` to compare with
R.  Fabel's default sums the criteria of all variables.  ``landmarkreg`` does
not accept multivariate curves; its warps depend only on the landmarks, so the
golden file runs it on each variable with the same landmarks.

As in ``test_registration.py`` three kinds of comparison are made: the golden
call itself, R's post-processing in isolation (R's ``Wfd`` in, R's warps and
registered curves out, ``max_iter=0``), and the measured facts behind every
strict xfail (R stops short of its optimum; R integrates its warps with a
1025-point trapezoid rule), each checked by a passing test.
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel.basis import Basis
from fabel.core import FData
from fabel.registration import RegistrationResult, register

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "registration_multivariate"
CASES: dict[str, dict[str, Any]] = {case["name"]: case for case in golden_cases(MODULE)}

GAIT = "register_fd_gait_hip_knee_crit2"
GAIT_PERIODIC = "register_fd_gait_hip_knee_periodic_crit1"
HANDWRITING = "register_fd_handwriting_xy_crit2"
LANDMARK = "landmarkreg_gait_hip_knee_per_variable"
CONTINUOUS = (GAIT, GAIT_PERIODIC, HANDWRITING)

#: R's register.fd fits multivariate warps to the first variable only.
R_VAR_WEIGHTS = (1.0, 0.0)

#: R's warp values: a trapezoid rule on this many equally spaced points,
#: linearly interpolated (measured, see ``test_r_warps_use_a_trapezoid_rule``).
R_TRAPEZOID_POINTS = 1025

# --------------------------------------------------------------------------- #
# measured reasons for the fields R gets wrong
# --------------------------------------------------------------------------- #

_STOPS_SHORT = {
    GAIT: (
        "register.fd returns a point that is not a minimum of its own criterion. "
        "With R's discretisation (grid mean over 211 points, crit=2, lambda=0.01, "
        "first variable), R's Wfd scores 4.005844 / 1.376484 / 1.724921 / 3.736952 "
        "/ 7.200670 for the five boys, with gradients of max-norm 0.19 to 2.09 "
        "there; Fabel's Newton iterate is stationary (gradient < 1e-12) at "
        "3.993467 / 1.369739 / 1.716418 / 3.722099 / 6.791805 -- lower on every "
        "curve (checked by test_fabel_optimum_is_below_r_point).  The optima "
        "differ by up to 2.01 in W."
    ),
    GAIT_PERIODIC: (
        "register.fd returns a point that is not a minimum of its own criterion. "
        "With R's discretisation (grid mean over 211 points, crit=1, lambda=0.01, "
        "periodic shift, first variable), R's (Wfd, shift) scores 4.772789 / "
        "4.817748 / 13.194021 / 18.376605 for the four boys, with gradients of "
        "max-norm 1.09 to 6.20 there; Fabel's Newton iterate is stationary "
        "(gradient < 1e-12) at 4.578179 / 3.473667 / 12.702899 / 13.294800 -- "
        "lower on every curve (checked by test_fabel_optimum_is_below_r_point). "
        "Shifts differ by up to 2.32 and W by up to 7.6."
    ),
    HANDWRITING: (
        "register.fd returns a point that is not a minimum of its own criterion. "
        "With R's discretisation (grid mean over 211 points, crit=2, lambda=1, "
        "first variable), R's Wfd scores 0.495565 / 0.367009 / 1.008572 / "
        "0.924378 / 1.371460 for the five samples, with gradients of max-norm "
        "0.064 to 0.26 there; Fabel's Newton iterate is stationary (gradient < "
        "1e-11) at 0.495236 / 0.366892 / 1.008538 / 0.924099 / 1.371380 -- lower "
        "on every curve (checked by test_fabel_optimum_is_below_r_point).  The "
        "optima differ by up to 0.0073 in W."
    ),
}

_TRAPEZOID = (
    "R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and "
    "linear interpolation; that rule reproduces R's warpfd and register.newfd "
    "output to 1e-14 (test_r_warps_use_a_trapezoid_rule).  It is off the exact "
    "integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 "
    "(handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 "
    "/ 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / "
    "5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the "
    "smallest coefficients.  Fabel integrates exp(W) by Gauss-Legendre to "
    "rounding error."
)

_R_REGFD_NOT_ITS_WARP = (
    "R's regfd is not x(h(t)) for the warp R returns: it differs from R's own "
    "register.newfd(yfd, Wfd) -- x at R's trapezoid warp, reproduced to 1e-14 "
    "-- by 4.66e-3 (gait) / 1.02e-4 (handwriting) normwise relative, and from "
    "x(h(t) + shift) at R's trapezoid warp by 3.99e-3 (periodic gait).  The "
    "univariate cases show the same internal warp defect (test_registration.py)."
)

_LANDMARK_STOPS_SHORT = (
    "landmarkreg's Wfd does not meet the landmarks it is fitted to: its warps "
    "miss the target landmark by 3.9e-6 / 3.9e-7 / 1.3e-7 / 1.6e-6 / 3.3e-7 / "
    "1.1e-5 for the six boys (smooth.morph objective 1.5e-11 to 1.3e-10), while "
    "Fabel's meet it to 1e-10 (objective < 1e-23, checked by "
    "test_fabel_landmark_warps_meet_the_landmarks).  The Wfd of the boy whose "
    "landmark is 0.01 from the target is 2.8e-3 in size, so its 1e-6 gap is "
    "3.9e-4 relative; the largest gap is 3.4e-6 in W."
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    **{
        (name, field): _STOPS_SHORT[name]
        for name in CONTINUOUS
        for field in ("Wfd_coefs", "warpfd_coefs", "regfd_coefs")
    },
    (GAIT_PERIODIC, "shift"): _STOPS_SHORT[GAIT_PERIODIC],
}

R_POSTPROCESS_DEFECTS: dict[tuple[str, str], str] = {
    **{(name, "warpfd_coefs"): _TRAPEZOID for name in CONTINUOUS},
    **{(name, "regfd_coefs"): _R_REGFD_NOT_ITS_WARP for name in CONTINUOUS},
    **{(name, "newfd_coefs"): _TRAPEZOID for name in CONTINUOUS},
}


# --------------------------------------------------------------------------- #
# case inputs
# --------------------------------------------------------------------------- #


def _basis(case: dict[str, Any], key: str = "basis") -> Basis:
    return build_basis(case["input"][key])


def _curves(case: dict[str, Any], key: str = "coefs") -> FData:
    """Return the curves recorded under ``input[key]`` in the case's curve basis."""
    return FData(np.asarray(case["input"][key], dtype=float), _basis(case))


def _options(case: dict[str, Any]) -> dict[str, Any]:
    """Translate the recorded ``register.fd`` arguments into ``register`` options."""
    inputs = case["input"]
    return {
        "warp_basis": _basis(case, "wbasis"),
        "lam": float(inputs["lambda"]),
        "penalty": int(inputs["penalty"]),
        "criterion": "eigen" if inputs["crit"] == 2 else "least_squares",
        "periodic": bool(inputs["periodic"]),
    }


def _r_latent(case: dict[str, Any]) -> FData:
    return FData(np.asarray(case["output"]["Wfd_coefs"], dtype=float), _basis(case, "wbasis"))


@cache
def _fabel(name: str) -> RegistrationResult:
    """Fabel on the golden inputs with R's first-variable weighting."""
    case = CASES[name]
    return register(
        _curves(case), _curves(case, "y0fd_coefs"), var_weights=R_VAR_WEIGHTS, **_options(case)
    )


@cache
def _from_r_latent(name: str, *, with_shift: bool = True) -> RegistrationResult:
    """Apply R's own ``Wfd`` (and shift) without optimising."""
    case = CASES[name]
    options = _options(case)
    shift = np.asarray(case["output"]["shift"], dtype=float)
    return register(
        _curves(case),
        _curves(case, "y0fd_coefs"),
        init=_r_latent(case),
        init_shift=shift if with_shift else None,
        max_iter=0,
        var_weights=R_VAR_WEIGHTS,
        **{**options, "periodic": options["periodic"] and with_shift},
    )


def _field(result: RegistrationResult, field: str) -> Any:
    if field == "regfd_coefs":
        return result.registered.coefs
    if field == "warpfd_coefs":
        return result.warp.coefs
    if field == "Wfd_coefs":
        return result.latent.coefs
    if field == "shift":
        return result.shift
    raise AssertionError(f"unknown field {field!r}")


def _mark(request: pytest.FixtureRequest, reason: str | None) -> None:
    if reason is not None:
        request.applymarker(pytest.mark.xfail(strict=True, reason=reason))


def _trapezoid_warp(case: dict[str, Any], grid: np.ndarray) -> np.ndarray:
    """R's warp values at ``grid``: trapezoid rule on 1025 points, linear interpolation."""
    wbasis = _basis(case, "wbasis")
    lower, upper = wbasis.domain
    latent = np.asarray(case["output"]["Wfd_coefs"], dtype=float)
    fine = np.linspace(lower, upper, R_TRAPEZOID_POINTS)
    height = np.exp(np.asarray(wbasis(fine)) @ latent)
    steps = 0.5 * (height[1:] + height[:-1]) * (fine[1] - fine[0])
    cumulative = np.concatenate([np.zeros((1, latent.shape[1])), np.cumsum(steps, axis=0)])
    warp = lower + (upper - lower) * cumulative / cumulative[-1]
    return np.stack([np.interp(grid, fine, warp[:, i]) for i in range(latent.shape[1])], axis=1)


def _grid(fd: FData) -> np.ndarray:
    return np.linspace(*fd.domain, max(201, 10 * fd.basis.n_basis + 1))


FIELDS = ("regfd_coefs", "warpfd_coefs", "Wfd_coefs", "shift")


# --------------------------------------------------------------------------- #
# continuous registration
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "field"), [(n, f) for n in CONTINUOUS for f in FIELDS])
def test_register_matches_r(name: str, field: str, request: pytest.FixtureRequest) -> None:
    """``register(..., var_weights=[1, 0])`` reproduces multivariate ``register.fd``."""
    _mark(request, R_FDA_DEFECTS.get((name, field)))
    case = CASES[name]
    compare(_field(_fabel(name), field), case["output"][field], case_rtol(MODULE, case))


@pytest.mark.parametrize("name", CONTINUOUS)
def test_r_fits_multivariate_warps_to_the_first_variable(name: str) -> None:
    """R's multivariate warps are its univariate warps of the first variable."""
    output = CASES[name]["output"]
    compare(output["Wfd_coefs"], output["Wfd_first_variable"], 1e-8)
    compare(output["shift"], output["shift_first_variable"], 1e-8)


@pytest.mark.parametrize("name", CONTINUOUS)
def test_first_variable_weights_are_the_univariate_registration(name: str) -> None:
    """``var_weights=[1, 0]`` is the registration of the first variable alone."""
    case = CASES[name]
    curves, target = _curves(case), _curves(case, "y0fd_coefs")
    first = register(
        FData(curves.coefs[:, :, 0], curves.basis),
        FData(target.coefs[:, :, 0], target.basis),
        **_options(case),
    )
    multi = _fabel(name)
    np.testing.assert_allclose(multi.latent.coefs, first.latent.coefs, rtol=0, atol=1e-13)
    np.testing.assert_allclose(multi.shift, first.shift, rtol=0, atol=1e-13)
    np.testing.assert_allclose(
        multi.registered.coefs[:, :, 0], first.registered.coefs, rtol=0, atol=1e-10
    )
    second = multi.apply(FData(curves.coefs[:, :, 1], curves.basis))
    np.testing.assert_allclose(multi.registered.coefs[:, :, 1], second.coefs, rtol=0, atol=1e-10)


@pytest.mark.parametrize("name", CONTINUOUS)
def test_fabel_optimum_is_below_r_point(name: str) -> None:
    """The measured fact behind the xfails: R stops above a stationary point."""
    at_r = _from_r_latent(name)
    fabel = _fabel(name)
    assert at_r.criterion is not None
    assert fabel.criterion is not None
    assert np.all(np.asarray(fabel.criterion) < np.asarray(at_r.criterion))
    restarted = register(
        _curves(CASES[name]),
        _curves(CASES[name], "y0fd_coefs"),
        init=_r_latent(CASES[name]),
        init_shift=np.asarray(CASES[name]["output"]["shift"], dtype=float),
        var_weights=R_VAR_WEIGHTS,
        **_options(CASES[name]),
    )
    assert np.all(np.asarray(restarted.n_iter) > 0)


@pytest.mark.parametrize(
    ("name", "field"), [(n, f) for n in CONTINUOUS for f in ("warpfd_coefs", "regfd_coefs")]
)
def test_warps_from_r_latent_match_r(name: str, field: str, request: pytest.FixtureRequest) -> None:
    """R's own ``Wfd`` maps to R's ``warpfd`` / ``regfd`` through Fabel's warps."""
    _mark(request, R_POSTPROCESS_DEFECTS.get((name, field)))
    case = CASES[name]
    compare(_field(_from_r_latent(name), field), case["output"][field], case_rtol(MODULE, case))


def test_from_r_latent_keeps_r_latent_and_shift() -> None:
    """``max_iter=0`` leaves R's multivariate latent functions and shifts untouched."""
    for name in CONTINUOUS:
        result = _from_r_latent(name)
        compare(result.latent.coefs, CASES[name]["output"]["Wfd_coefs"], 1e-14)
        compare(result.shift, CASES[name]["output"]["shift"], 1e-14)
        assert result.registered.coefs.shape == np.shape(CASES[name]["output"]["regfd_coefs"])


# --------------------------------------------------------------------------- #
# register.newfd
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", CONTINUOUS)
def test_apply_matches_register_newfd(name: str, request: pytest.FixtureRequest) -> None:
    """``RegistrationResult.apply`` reproduces ``register.newfd(yfd, Wfd, 'monotone')``.

    R's ``type='monotone'`` ignores the shifts of a periodic registration, so
    R's latent functions are applied without them.
    """
    _mark(request, R_POSTPROCESS_DEFECTS.get((name, "newfd_coefs")))
    case = CASES[name]
    result = _from_r_latent(name, with_shift=False)
    compare(
        result.apply(_curves(case)).coefs, case["output"]["newfd_coefs"], case_rtol(MODULE, case)
    )


@pytest.mark.parametrize("name", CONTINUOUS)
def test_r_warps_use_a_trapezoid_rule(name: str) -> None:
    """The measured fact behind the trapezoid xfails, reproduced to rounding."""
    case = CASES[name]
    curves = _curves(case)
    grid = _grid(curves)
    trapezoid = _trapezoid_warp(case, grid)
    shift = np.asarray(case["output"]["shift"], dtype=float)
    wbasis = _basis(case, "wbasis")
    warp = np.linalg.lstsq(np.asarray(wbasis(grid)), trapezoid + shift, rcond=None)[0]
    compare(warp, case["output"]["warpfd_coefs"], 1e-12)
    values = np.stack(
        [np.asarray(curves[i](trapezoid[:, i]))[:, 0, :] for i in range(curves.n_curves)], axis=1
    )
    n_grid, n_curves, n_vars = values.shape
    design = np.asarray(curves.basis(grid))
    newfd = np.linalg.lstsq(design, values.reshape(n_grid, n_curves * n_vars), rcond=None)[0]
    compare(newfd.reshape(-1, n_curves, n_vars), case["output"]["newfd_coefs"], 1e-12)
    exact = np.asarray(_from_r_latent(name, with_shift=False).warp_values(grid))
    assert 1e-6 < np.max(np.abs(exact - trapezoid)) < 1e-4


# --------------------------------------------------------------------------- #
# landmark registration
# --------------------------------------------------------------------------- #


@cache
def _landmark_result() -> RegistrationResult:
    case = CASES[LANDMARK]
    inputs = case["input"]
    return register(
        _curves(case),
        landmarks=np.asarray(inputs["ximarks"], dtype=float)[:, None],
        target_landmarks=[inputs["x0marks"]],
        warp_basis=_basis(case, "wbasis"),
        lam=float(inputs["lambda"]),
        penalty=int(inputs["penalty"]),
    )


@pytest.mark.parametrize("field", ["Wfd_coefs_hip", "Wfd_coefs_knee"])
@pytest.mark.xfail(strict=True, reason=_LANDMARK_STOPS_SHORT)
def test_landmark_latent_matches_r(field: str) -> None:
    """Multivariate landmark warps reproduce ``landmarkreg`` on each variable."""
    case = CASES[LANDMARK]
    compare(_landmark_result().latent.coefs, case["output"][field], case_rtol(MODULE, case))


def test_r_landmark_warps_do_not_depend_on_the_variable() -> None:
    """R gives the hip and the knee the same landmark warps."""
    output = CASES[LANDMARK]["output"]
    compare(output["Wfd_coefs_hip"], output["Wfd_coefs_knee"], 1e-14)


def test_fabel_landmark_warps_meet_the_landmarks() -> None:
    """The measured fact behind the landmark xfail: Fabel's warps hit the landmarks."""
    case = CASES[LANDMARK]
    inputs = case["input"]
    result = _landmark_result()
    hit = np.asarray(result.warp_values(np.array([inputs["x0marks"]])))[0]
    np.testing.assert_allclose(hit, inputs["ximarks"], rtol=0, atol=1e-10)
    r_result = RegistrationResult(
        registered=result.registered,
        warp=result.warp,
        unregistered=result.unregistered,
        latent=FData(np.asarray(case["output"]["Wfd_coefs_hip"]), _basis(case, "wbasis")),
        shift=result.shift,
    )
    r_hit = np.asarray(r_result.warp_values(np.array([inputs["x0marks"]])))[0]
    assert np.max(np.abs(r_hit - np.asarray(inputs["ximarks"]))) > 1e-5
    r_latent = np.asarray(case["output"]["Wfd_coefs_hip"], dtype=float)
    np.testing.assert_allclose(result.latent.coefs, r_latent, rtol=0, atol=5e-6)


def test_multivariate_landmarks_register_every_variable() -> None:
    """Each variable of the multivariate result is the univariate registration."""
    case = CASES[LANDMARK]
    inputs = case["input"]
    curves = _curves(case)
    result = _landmark_result()
    for var in range(curves.n_vars):
        single = register(
            FData(curves.coefs[:, :, var], curves.basis),
            landmarks=np.asarray(inputs["ximarks"], dtype=float)[:, None],
            target_landmarks=[inputs["x0marks"]],
            warp_basis=_basis(case, "wbasis"),
            lam=float(inputs["lambda"]),
            penalty=int(inputs["penalty"]),
        )
        np.testing.assert_allclose(result.latent.coefs, single.latent.coefs, rtol=0, atol=0)
        np.testing.assert_allclose(
            result.registered.coefs[:, :, var], single.registered.coefs, rtol=0, atol=1e-12
        )
