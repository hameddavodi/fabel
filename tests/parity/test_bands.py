"""Parity of :func:`fabel.stats.confidence_band` against R ``fda`` 6.3.0.

``tests/golden/bands.json`` (``tools/golden_r/bands.R``) records pointwise
standard errors on ``CanadianWeather``:

* of a ``smooth.basis`` fit of log10 precipitation, ``sqrt(diag(Φ S Σ Sᵀ Φᵀ))``
  with ``S`` R's ``y2cMap`` and ``Σ`` either ``σ² I`` or the diagonal of the
  per-day residual variance across stations;
* of ``fRegress`` coefficients, ``sqrt(diag(Θ V Θᵀ))`` with ``V`` the ``bvar``
  R returns when given ``y2cMap`` and ``SigmaE``.

Fabel refits everything itself from the raw data.  The weather arrays come from
:func:`fabel.datasets.load_canadian_weather`; the tests are skipped when the
dataset is not in the local cache (``$FABEL_DATA_DIR``).
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel import Constant, FData
from fabel.datasets import CanadianWeather, _cache_dir, load_canadian_weather
from fabel.regression import FRegressResult, fregress
from fabel.smoothing import SmoothResult, smooth
from fabel.stats import confidence_band

from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "bands"
CASES = {case["name"]: case for case in golden_cases(MODULE)}

SMOOTH_SCALAR = "band_smooth_logprec_bspline53_sigma_scalar"
SMOOTH_POINTWISE = "band_smooth_logprec_bspline53_sigma_pointwise"
SCALAR_DESIGN = "band_fregress_logannualprec_on_meantemp_latitude"
FUNCTIONAL_DESIGN = "band_fregress_logannualprec_on_temp_bspline7"

_FUNCTIONAL_QUADRATURE = (
    "R integrates int temp_i(t) theta_k(t) dt approximately; Fabel's integrals "
    "are exact.  R's own dense trapezoid rule (365001 points) gives "
    "int temp_1 theta_4 = 943.854155886, Fabel's value to 1e-13, while R's "
    "inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the "
    "same kind of error (see tests/parity/test_regression.py).  Downstream, "
    "measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 "
    "1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the "
    "pointwise standard error of the temperature coefficient 4.9e-5."
)


# --------------------------------------------------------------------------- #
# data and fits
# --------------------------------------------------------------------------- #


@cache
def _weather() -> CanadianWeather:
    if not (_cache_dir() / "CanadianWeather.npz").exists():
        pytest.skip("CanadianWeather is not in the local dataset cache ($FABEL_DATA_DIR)")
    return load_canadian_weather()


def _annual_log_precip() -> np.ndarray:
    return np.log10(np.sum(_weather().precip, axis=0))


@cache
def _smooth_fit() -> SmoothResult:
    inputs = CASES[SMOOTH_SCALAR]["input"]
    return smooth(
        _weather().log10precip,
        np.asarray(inputs["argvals"]),
        basis=build_basis(inputs["basis"]),
        lam=inputs["lambda"],
        penalty=inputs["lfd"],
    )


def _residual_variance() -> np.ndarray:
    """Per-day residual variance across the 35 stations, as the book computes it."""
    inputs = CASES[SMOOTH_POINTWISE]["input"]
    fit = _smooth_fit()
    resid = _weather().log10precip - fit.fd(np.asarray(inputs["argvals"]))
    variance: np.ndarray = np.sum(resid**2, axis=1) / (resid.shape[1] - 1)
    return variance


@cache
def _scalar_model() -> FRegressResult:
    weather = _weather()
    return fregress(
        _annual_log_precip(),
        [1.0, np.mean(weather.temp, axis=0), weather.coordinates[:, 0]],
    )


@cache
def _functional_model() -> FRegressResult:
    inputs = CASES[FUNCTIONAL_DESIGN]["input"]
    temp = FData(np.asarray(inputs["temp_coefs"]), build_basis(inputs["temp_basis"]))
    # R's intercept is the constant *function* 1 on the year: its coefficient is
    # the intercept divided by the length of the domain.
    const = FData(np.ones((1, temp.n_curves)), Constant(domain=temp.domain))
    beta_basis = build_basis(inputs["beta_basis"])
    return fregress(
        np.asarray(inputs["y"]),
        [const, temp],
        beta=[None, (beta_basis, inputs["lambda"], 2)],
    )


def _sigma2(model: FRegressResult) -> float:
    assert model.df is not None
    resid = np.asarray(model.y) - np.asarray(model.fitted)
    return float(np.sum(resid**2)) / (int(resid.shape[0]) - float(model.df))


# --------------------------------------------------------------------------- #
# accessors
# --------------------------------------------------------------------------- #


def _smooth_scalar(field: str) -> Any:
    inputs = CASES[SMOOTH_SCALAR]["input"]
    t = np.asarray(inputs["t"])
    fit = _smooth_fit()
    band = confidence_band(fit, t, level=inputs["level"])
    accessors: dict[str, Callable[[], Any]] = {
        "df": lambda: fit.df,
        "SSE": lambda: fit.sse,
        "sigma2": lambda: fit.sse / (35 * (365 - fit.df)),
        "stderr": lambda: band.stderr,
        "stderr_deriv1": lambda: confidence_band(fit, t, deriv=1).stderr,
        "z": lambda: (band.upper[0, 0] - band.estimate[0, 0]) / band.stderr[0],
        "lower_curve1": lambda: band.lower[:, 0],
        "upper_curve1": lambda: band.upper[:, 0],
    }
    return accessors[field]()


def _smooth_pointwise(field: str) -> Any:
    t = np.asarray(CASES[SMOOTH_POINTWISE]["input"]["t"])
    if field == "varvec":
        return _residual_variance()
    return confidence_band(_smooth_fit(), t, sigma_e=_residual_variance()).stderr


def _scalar_design(field: str) -> Any:
    model = _scalar_model()
    level = CASES[SCALAR_DESIGN]["input"]["level"]
    bands = confidence_band(model, np.array([0.5]), level=level)
    accessors: dict[str, Callable[[], Any]] = {
        "beta": lambda: [band.estimate[0] for band in bands],
        "sigma2": lambda: _sigma2(model),
        "bvar": lambda: model.stderr().cov,
        "stderr": lambda: [band.stderr[0] for band in bands],
        "lower": lambda: [band.lower[0] for band in bands],
        "upper": lambda: [band.upper[0] for band in bands],
    }
    return accessors[field]()


def _functional_design(field: str) -> Any:
    model = _functional_model()
    t = np.asarray(CASES[FUNCTIONAL_DESIGN]["input"]["t"])
    const, temp = confidence_band(model, t)
    accessors: dict[str, Callable[[], Any]] = {
        "beta1": lambda: temp.estimate,
        "df": lambda: model.df,
        "sigma2": lambda: _sigma2(model),
        "bvar": lambda: model.stderr().cov,
        "stderr0": lambda: const.stderr[0],
        "stderr1": lambda: temp.stderr,
    }
    return accessors[field]()


_ACCESSORS = {
    SMOOTH_SCALAR: _smooth_scalar,
    SMOOTH_POINTWISE: _smooth_pointwise,
    SCALAR_DESIGN: _scalar_design,
    FUNCTIONAL_DESIGN: _functional_design,
}

#: ``(case, field)`` pairs where R is the inexact side, with measured numbers.
R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    (FUNCTIONAL_DESIGN, field): _FUNCTIONAL_QUADRATURE
    for field in ("beta1", "df", "sigma2", "bvar", "stderr0", "stderr1")
}


def _params() -> list[Any]:
    params = []
    for name, case in CASES.items():
        for field in case["output"]:
            reason = R_FDA_DEFECTS.get((name, field))
            marks = [pytest.mark.xfail(reason=reason, strict=True)] if reason else []
            params.append(pytest.param(name, field, id=f"{name}-{field}", marks=marks))
    return params


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "field"), _params())
def test_golden_field(name: str, field: str) -> None:
    case = CASES[name]
    want = np.asarray(case["output"][field], dtype=float)
    actual = np.asarray(_ACCESSORS[name](field), dtype=float)
    compare(np.reshape(actual, want.shape), want, case_rtol(MODULE, case))


def test_every_defect_names_a_real_field() -> None:
    known = {(name, field) for name, case in CASES.items() for field in case["output"]}
    assert set(R_FDA_DEFECTS) <= known


def test_every_case_has_an_accessor() -> None:
    assert set(CASES) == set(_ACCESSORS)
