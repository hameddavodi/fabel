"""Parity of :mod:`fabel.regression` against golden output from R ``fda`` 6.3.0.

Two references are used.

``tests/golden/regression.json`` records ``fRegress`` and its satellites on the
three book case studies plus a synthetic design.  Cases are parametrised by
``(case, field)``: every quantity R returns is asserted on its own, so a defect
in one does not hide agreement in the others.  The weather cases are built on
R's *own* smooths of the data -- the golden ``smooth.basis`` data-to-coefficient
map of the same design (``tests/golden/smoothing.json``) applied to the raw
observations -- so that the known defect of R's harmonic-accelerator penalty
(see ``tests/parity/test_smoothing.py``) does not leak into the regression
comparison.

Every golden case puts a functional term into the design, and R's ``fRegress``
evaluates the integrals of such a term by an approximate quadrature; Fabel's are
exact.  Where the difference exceeds the 1e-8 tolerance the field is a strict
xfail carrying the measured error of R's integrals.  To still check the rest of
the pipeline -- solve, df, GCV, OCV, coefficient covariance, standard errors,
cross-validation and prediction -- at 1e-8, ``tests/fixtures/
regression_scalar_design.json`` records R on a design of scalar covariates
only, which R computes without any quadrature.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from fabel import LDO, Constant, FData, Fourier
from fabel.regression import FRegressResult, fregress

from ._datasets import weather_daily, weather_region
from .conftest import build_basis, case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "regression"
CASES = {case["name"]: case for case in golden_cases(MODULE)}
FIXTURE = Path(__file__).resolve().parent.parent / "fixtures" / "regression_scalar_design.json"

WEATHER_SMOOTH_CASE = "smooth_basis_weather_fourier65_harmonic_lambda100"
PERIOD = 365.0

# --------------------------------------------------------------------------- #
# R defects
# --------------------------------------------------------------------------- #

_SCALAR_QUADRATURE = (
    "R integrates int x_i(t) theta_k(t) dt approximately.  The response basis is "
    "an orthonormal Fourier basis, so the exact integral is simply the k-th "
    "coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out "
    "and its Cmat up to 8.2e-8.  Fabel's Cmat and Dmat are the exact values."
)
_REGION_QUADRATURE = (
    "R's Gram matrix of the Fourier(65) coefficient basis is not the identity: "
    "the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = "
    "0.9999986 for an orthonormal basis.  With it R's Dmat is 1.38e-6 and its "
    "betas 1.13e-6 out (relative to the largest entry); entry by entry the "
    "near-zero coefficients differ by up to 1.1e-2 relative."
)
_CONCURRENT_QUADRATURE = (
    "R's concurrent-model integrals are 0.1-0.3% off.  R's own dense periodic "
    "trapezoid rule (20000 points, exact for trigonometric polynomials) gives "
    "sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i "
    "logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) "
    "and Dmat[6] 4172.91107 (3.4e-3 low).  Fabel reproduces the dense values "
    "(5861.19665, 4187.25630).  The betas inherit a 1.2e-2 relative error."
)
_SYNTHETIC_QUADRATURE = (
    "R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to "
    "1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact "
    "inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, "
    "which is Fabel's value.  Every quantity downstream of the design inherits "
    "the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4)."
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    **{
        ("fregress_scalar_precip_on_temp", field): _SCALAR_QUADRATURE
        for field in ("beta1", "OCV", "Cmat", "Dmat")
    },
    **{
        ("fregress_cv_scalar_precip_on_temp", field): _SCALAR_QUADRATURE
        for field in ("SSE.CV", "errfd.cv")
    },
    **{
        ("fregress_functional_temp_on_region", field): _REGION_QUADRATURE
        for field in ("beta0", "beta1", "beta2", "beta3", "yhatfdobj_coefs", "Cmat", "Dmat")
    },
    **{
        ("fregress_concurrent_logprecip_on_temp", field): _CONCURRENT_QUADRATURE
        for field in ("beta0", "beta1", "yhatfdobj_coefs", "Cmat", "Dmat")
    },
    **{
        ("fregress_synthetic_2scalar_1functional", field): _SYNTHETIC_QUADRATURE
        for field in (
            "beta0",
            "beta1",
            "beta2",
            "beta3",
            "yhatfdobj",
            "df",
            "gcv",
            "OCV",
            "Cmat",
            "Dmat",
        )
    },
    **{
        # At lambda = 0 the scalar coefficients beta0..beta2 still agree.  R's
        # approximate integrals are a fixed linear map of the exact ones
        # (Z_R = X' G_R against Z = X' G, G_R an approximate cross-Gram), so the
        # design spans the same column space and the unpenalised least-squares
        # fit only changes the functional coefficient beta3.
        ("fregress_synthetic_2scalar_1functional_lambda0", field): _SYNTHETIC_QUADRATURE
        for field in ("beta3", "yhatfdobj", "gcv", "OCV", "Cmat", "Dmat")
    },
}

#: Golden cases whose recorded output carries no numbers: the generator read
#: ``b$fd$coefs`` from each ``betastderrlist`` entry, but those entries are
#: plain ``fd`` objects, so every value was recorded as ``null``.
NULL_OUTPUT_CASES = frozenset(
    {"fregress_stderr_scalar_precip_on_temp", "fregress_stderr_functional_temp_on_region"}
)

# --------------------------------------------------------------------------- #
# model construction
# --------------------------------------------------------------------------- #


@cache
def _weather_y2c() -> np.ndarray:
    """R's data-to-coefficient map of the Fourier(65) harmonic smooth, lambda 1e2."""
    case = next(c for c in golden_cases("smoothing") if c["name"] == WEATHER_SMOOTH_CASE)
    return np.asarray(case["output"]["y2cMap"], dtype=float)


def _weather_fd(values: np.ndarray) -> FData:
    return FData(_weather_y2c() @ values, Fourier(domain=(0.0, PERIOD), n_basis=65))


@cache
def _scalar_model() -> FRegressResult:
    case = CASES["fregress_scalar_precip_on_temp"]
    temp = _weather_fd(weather_daily("Temperature.C"))
    const = FData(np.ones((1, temp.n_curves)), build_basis(case["input"]["basis1"]))
    beta_basis = build_basis(case["input"]["basis2"])
    return fregress(
        np.asarray(case["input"]["y"]),
        [const, temp],
        beta=[None, (beta_basis, case["input"]["lambda"], LDO.harmonic(PERIOD))],
    )


@cache
def _region_model() -> FRegressResult:
    case = CASES["fregress_functional_temp_on_region"]
    temp = _weather_fd(weather_daily("Temperature.C"))
    region = weather_region()
    dummies = {name: (region == name).astype(float) for name in case["input"]["regions"]}
    return fregress(temp, dummies, lam=case["input"]["lambda"], penalty=LDO.harmonic(PERIOD))


@cache
def _concurrent_model() -> FRegressResult:
    case = CASES["fregress_concurrent_logprecip_on_temp"]
    temp = _weather_fd(weather_daily("Temperature.C"))
    logprecip = _weather_fd(np.log(np.maximum(weather_daily("Precipitation.mm"), 0.01)))
    basis = build_basis(case["input"]["basis"])
    spec = (basis, case["input"]["lambda"], 2)
    return fregress(logprecip, [1.0, temp], beta=[spec, spec])


@cache
def _synthetic_model(name: str) -> FRegressResult:
    inputs = CASES[name]["input"]
    basis = build_basis(inputs["basis"])
    x = FData(np.asarray(inputs["coefs"]), basis)
    return fregress(
        np.asarray(inputs["y"]),
        [1.0, np.asarray(inputs["z1"]), np.asarray(inputs["z2"]), x],
        beta=[None, None, None, (basis, inputs["lambda"], 2)],
    )


def _fit_fields(model: Callable[[], FRegressResult]) -> dict[str, Callable[[], Any]]:
    """Accessors for the fields ``fregress_output()`` records for one fit."""

    def beta(j: int) -> Callable[[], Any]:
        return lambda: model().beta[j].coefs

    def fitted() -> Any:
        result = model().fitted
        return result.coefs if isinstance(result, FData) else result

    fields: dict[str, Callable[[], Any]] = {
        "Cmat": lambda: model().cmat,
        "Dmat": lambda: model().dmat,
        "df": lambda: model().df,
        "gcv": lambda: model().gcv,
        "OCV": lambda: model().ocv,
        "yhatfdobj": fitted,
        "yhatfdobj_coefs": fitted,
    }
    fields.update({f"beta{j}": beta(j) for j in range(4)})
    return fields


_MODELS: dict[str, Callable[[], FRegressResult]] = {
    "fregress_scalar_precip_on_temp": _scalar_model,
    "fregress_functional_temp_on_region": _region_model,
    "fregress_concurrent_logprecip_on_temp": _concurrent_model,
    "fregress_synthetic_2scalar_1functional": lambda: _synthetic_model(
        "fregress_synthetic_2scalar_1functional"
    ),
    "fregress_synthetic_2scalar_1functional_lambda0": lambda: _synthetic_model(
        "fregress_synthetic_2scalar_1functional_lambda0"
    ),
}


def _expected_fields(case: dict[str, Any]) -> dict[str, Any]:
    """Flatten a golden case's output: ``betaestlist_coefs`` becomes ``beta0..``."""
    out = dict(case["output"])
    betas = out.pop("betaestlist_coefs", None)
    if betas is not None:
        out.update({f"beta{j}": value for j, value in enumerate(betas)})
    return {key: value for key, value in out.items() if value is not None}


def _accessors(name: str) -> dict[str, Callable[[], Any]]:
    if name in _MODELS:
        return _fit_fields(_MODELS[name])
    if name == "predict_fregress_scalar_precip_on_temp":
        return {"predicted": lambda: _scalar_model().predict()}
    if name == "fregress_cv_scalar_precip_on_temp":
        return {
            "SSE.CV": lambda: _scalar_model().cv().sse,
            "errfd.cv": lambda: _scalar_model().cv().errors,
        }
    raise AssertionError(f"no accessor for golden case {name!r}")


def _params() -> list[Any]:
    params = []
    for name, case in CASES.items():
        if name in NULL_OUTPUT_CASES:
            continue
        for field in _expected_fields(case):
            marks = []
            reason = R_FDA_DEFECTS.get((name, field))
            if reason is not None:
                marks.append(pytest.mark.xfail(reason=reason, strict=True))
            params.append(pytest.param(name, field, id=f"{name}-{field}", marks=marks))
    return params


# --------------------------------------------------------------------------- #
# golden parity
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "field"), _params())
def test_golden_field(name: str, field: str) -> None:
    case = CASES[name]
    expected = _expected_fields(case)[field]
    actual = _accessors(name)[field]()
    want = np.asarray(expected, dtype=float)
    compare(np.reshape(np.asarray(actual, dtype=float), want.shape), want, case_rtol(MODULE, case))


@pytest.mark.parametrize("name", sorted(NULL_OUTPUT_CASES))
def test_null_golden_cases_really_are_empty(name: str) -> None:
    """Pin the reason these cases are not compared: every recorded value is null.

    If the golden file is regenerated with real numbers this fails, and the case
    must then move into the field-by-field comparison above.
    """
    values = CASES[name]["output"]["betastderrlist_coefs"]
    assert values
    assert all(value is None for value in values)


def test_every_defect_names_a_real_field() -> None:
    known = {(name, field) for name, case in CASES.items() for field in _expected_fields(case)}
    assert set(R_FDA_DEFECTS) <= known


# --------------------------------------------------------------------------- #
# scalar-only design: R computes it without quadrature
# --------------------------------------------------------------------------- #


@cache
def _fixture() -> dict[str, Any]:
    with FIXTURE.open(encoding="utf-8") as fh:
        doc: dict[str, Any] = json.load(fh)
    return doc


@cache
def _fixture_model() -> FRegressResult:
    inputs = _fixture()["input"]
    return fregress(
        np.asarray(inputs["y"]), [1.0, np.asarray(inputs["z1"]), np.asarray(inputs["z2"])]
    )


def _fixture_actual(field: str) -> Any:
    model = _fixture_model()
    inputs = _fixture()["input"]
    accessors: dict[str, Callable[[], Any]] = {
        "beta": lambda: np.concatenate([b.coefs[:, 0] for b in model.beta]),
        "yhat": lambda: model.fitted,
        "df": lambda: model.df,
        "gcv": lambda: model.gcv,
        "ocv": lambda: model.ocv,
        "Cmat": lambda: model.cmat,
        "Dmat": lambda: model.dmat,
        "sigma2": lambda: float(np.sum((model.y - model.fitted) ** 2)) / (30 - model.df),
        "bvar": lambda: model.stderr().cov,
        "betastderr": lambda: np.concatenate([b.coefs[:, 0] for b in model.stderr().beta]),
        "sse_cv": lambda: model.cv().sse,
        "errfd_cv": lambda: model.cv().errors,
        "predicted": lambda: model.predict(
            [1.0, np.asarray(inputs["z1_new"]), np.asarray(inputs["z2_new"])]
        ),
    }
    return accessors[field]()


@pytest.mark.parametrize("field", sorted(json.loads(FIXTURE.read_text("utf-8"))["output"]))
def test_scalar_design_matches_r(field: str) -> None:
    doc = _fixture()
    compare(_fixture_actual(field), doc["output"][field], doc["meta"]["rtol"])


def test_scalar_design_constant_basis_is_what_r_uses() -> None:
    model = _fixture_model()
    assert all(isinstance(b.basis, Constant) for b in model.beta)
