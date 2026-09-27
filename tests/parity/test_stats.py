"""Parity of :mod:`fabel.stats` against golden output from R ``fda`` 6.3.0.

Every case is built on the CanadianWeather temperatures, smoothed exactly as
``tools/golden_r/stats.R`` does: a 25-function Fourier basis on ``[0, 365]``,
observations at days ``1..365``, a second-derivative penalty and
``lambda = 1e2`` (``1e4`` for precipitation).

The two permutation tests are seeded in R (``set.seed(42)`` / ``set.seed(11)``).
Fabel draws permutations from whatever ``random_state`` it is given, so the
tests pass :class:`tests.parity._rrng.RRandom`, which answers ``permutation(n)``
with the very permutation R's ``sample(n)`` returns.  That turns the null
distributions into deterministic, comparable outputs.
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pytest

from fabel import FData, Fourier
from fabel.regression import fregress
from fabel.smoothing import smooth
from fabel.stats import boxplot, cor, cov, depth, f_test, t_test

from ._datasets import weather_daily, weather_region
from ._rrng import RRandom
from .conftest import case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "stats"

CASES = {case["name"]: case for case in golden_cases(MODULE)}

_FPERM_REASON = (
    "R's fRegress assembles its normal equations with inprod(), whose Romberg "
    "quadrature is inexact on a Fourier basis: R's inprod(basis, basis) for the "
    "25-function Fourier basis on [0, 365] departs from the identity -- the exact "
    "Gram matrix of that orthonormal basis -- by up to 1.383e-6, and the "
    "right-hand side Dmat carries a 7.28e-6 relative error.  Feeding R's own "
    "Cmat and Dmat (read from fRegress()$Cmat / $Dmat) through Fabel's solver and "
    "F statistic reproduces R's betas to 1.5e-15 and R's Fobs to 3.0e-15, so the "
    "whole gap is R's quadrature; with the exact Gram matrices the maximal "
    "F statistic is 0.44257145697174 against R's 0.44261503693042 (9.85e-5 "
    "relative).  Every refit under a permutation inherits the same error: Fnull "
    "is off by 3.9e-5 (median) and 2.5e-4 (worst) relative, and its 5% quantile "
    "qval by 1.1e-4.  The permutations themselves are R's exact draws (see "
    "RRandom), and pval, which only compares Fobs with Fnull, matches."
)

#: ``(case, field)`` pairs where R ``fda`` 6.3.0 is the less accurate side.
R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    ("fperm_fd_weather_temp_atlantic_dummy", "Fobs"): _FPERM_REASON,
    ("fperm_fd_weather_temp_atlantic_dummy", "Fnull"): _FPERM_REASON,
    ("fperm_fd_weather_temp_atlantic_dummy", "qval"): _FPERM_REASON,
}

WEATHER_BASIS = Fourier(domain=(0.0, 365.0), n_basis=25)
DAYS = np.arange(1.0, 366.0)


@cache
def weather_fd(variable: str = "Temperature.C", lam: float = 1e2) -> FData:
    """Return the smoothed weather curves the golden cases were computed from."""
    return smooth(weather_daily(variable), DAYS, basis=WEATHER_BASIS, lam=lam, penalty=2).fd


def _params(names: dict[str, tuple[str, ...]]) -> list[Any]:
    out = []
    for name, fields in names.items():
        for field in fields:
            reason = R_FDA_DEFECTS.get((name, field))
            marks = [pytest.mark.xfail(reason=reason, strict=True)] if reason else []
            out.append(pytest.param(name, field, id=f"{name}-{field}", marks=marks))
    return out


def _rtol(name: str) -> float:
    return case_rtol(MODULE, CASES[name])


# --------------------------------------------------------------------------- #
# the R generator stand-in
# --------------------------------------------------------------------------- #


def test_r_generator_matches_r() -> None:
    # Printed by R 4.6.1: set.seed(42); runif(3) / sample(20); set.seed(11);
    # sample(35); set.seed(-7); runif(2); set.seed(3); sample(100000, 3).
    rng = RRandom(42)
    assert [rng.unif() for _ in range(3)] == [
        0.91480604349635541,
        0.93707541329786181,
        0.28613953478634357,
    ]
    np.testing.assert_array_equal(
        RRandom(42).permutation(20) + 1,
        [17, 5, 1, 10, 4, 2, 20, 18, 8, 7, 16, 9, 19, 6, 14, 15, 12, 3, 13, 11],
    )
    np.testing.assert_array_equal(
        RRandom(11).permutation(35) + 1,
        [
            34,
            25,
            16,
            17,
            5,
            28,
            22,
            12,
            21,
            29,
            7,
            27,
            3,
            11,
            35,
            13,
            2,
            32,
            8,
            15,
            24,
            6,
            23,
            20,
            26,
            19,
            10,
            30,
            9,
            33,
            18,
            14,
            1,
            4,
            31,
        ],
    )
    rng = RRandom(-7)
    assert [rng.unif() for _ in range(2)] == [0.32649485557340086, 0.53986111702397466]
    np.testing.assert_array_equal(RRandom(3).sample(100000, 3) + 1, [52922, 87015, 84843])


# --------------------------------------------------------------------------- #
# var.fd / cor.fd / mean.fd / sd.fd
# --------------------------------------------------------------------------- #


def test_var_fd() -> None:
    name = "var_fd_weather_temperature"
    case = CASES[name]
    grid = np.asarray(case["input"]["grid"], dtype=float)
    temp = weather_fd()
    compare(cov(temp)(grid, grid), case["output"]["varmat"], _rtol(name))
    compare(cov(temp, temp)(grid, grid), case["output"]["varmat"], _rtol(name))


def test_cor_fd() -> None:
    name = "cor_fd_weather_temp_vs_precip"
    case = CASES[name]
    grid = np.asarray(case["input"]["grid"], dtype=float)
    temp = weather_fd()
    precip = weather_fd("Precipitation.mm", 1e4)
    compare(cor(temp, precip, s=grid, t=grid), case["output"]["cormat"], _rtol(name))


@pytest.mark.parametrize("field", ["mean_coefs", "mean_fitted", "sd_coefs", "sd_fitted"])
def test_mean_sd_fd(field: str) -> None:
    name = "mean_sd_fd_weather_temperature"
    case = CASES[name]
    grid = np.asarray(case["input"]["grid"], dtype=float)
    temp = weather_fd()
    summary = temp.mean() if field.startswith("mean") else temp.std()
    actual = summary.coefs if field.endswith("coefs") else summary(grid)
    compare(actual, case["output"][field], _rtol(name))


# --------------------------------------------------------------------------- #
# fbplot / fdepth
# --------------------------------------------------------------------------- #


def test_fbplot_mbd() -> None:
    name = "fbplot_weather_temperature_mbd"
    case = CASES[name]
    values = np.asarray(case["input"]["Y"], dtype=float)
    result = boxplot(values, method="MBD")
    compare(result.depth, case["output"]["depth"], _rtol(name))
    np.testing.assert_array_equal(result.outliers + 1, case["output"]["outpoint"])
    assert result.median_index + 1 == case["output"]["medcurve"]


def test_fbplot_input_is_the_smoothed_weather() -> None:
    # The golden fbplot/fdepth input is eval.fd(grid, temp_fd): checking it pins
    # the smoothing the other cases rely on.
    case = CASES["fbplot_weather_temperature_mbd"]
    grid = np.asarray(case["input"]["grid"], dtype=float)
    compare(weather_fd()(grid), case["input"]["Y"], _rtol(case["name"]))


@pytest.mark.parametrize("field", ["prof", "median", "lmed", "mtrim", "ltrim"])
def test_fdepth_fm(field: str) -> None:
    name = "fdepth_weather_temperature_fm"
    case = CASES[name]
    values = np.asarray(case["input"]["Y"], dtype=float)
    result = depth(values, method="FM")
    expected = case["output"][field]
    if field == "prof":
        compare(result.depth, expected, _rtol(name))
    elif field == "median":
        compare(result.median, expected, _rtol(name))
    elif field == "mtrim":
        compare(result.trimmed_mean, expected, _rtol(name))
    elif field == "lmed":
        assert result.median_index + 1 == expected
    else:
        np.testing.assert_array_equal(result.trimmed_index + 1, expected)


def test_fdepth_accepts_fdata() -> None:
    case = CASES["fdepth_weather_temperature_fm"]
    grid = np.asarray(case["input"]["grid"], dtype=float)
    result = depth(weather_fd(), method="FM", t=grid)
    compare(result.depth, case["output"]["prof"], _rtol(case["name"]))


# --------------------------------------------------------------------------- #
# tperm.fd / Fperm.fd
# --------------------------------------------------------------------------- #


@cache
def _tperm() -> Any:
    case = CASES["tperm_fd_weather_atlantic_vs_pacific"]
    region = weather_region()
    temp = weather_fd()
    return t_test(
        temp[region == "Atlantic"],
        temp[region == "Pacific"],
        n_perm=case["input"]["nperm"],
        q=case["input"]["q"],
        random_state=RRandom(case["input"]["seed"]),
    )


@cache
def _fperm() -> Any:
    case = CASES["fperm_fd_weather_temp_atlantic_dummy"]
    atlantic = np.asarray(case["input"]["atlantic"], dtype=float)
    return f_test(
        weather_fd(),
        [np.ones(atlantic.shape[0]), atlantic],
        basis=WEATHER_BASIS,
        lam=1e2,
        penalty=2,
        n_perm=case["input"]["nperm"],
        q=case["input"]["q"],
        random_state=RRandom(case["input"]["seed"]),
    )


_PERM_FIELDS = {
    "statistic": ("Tobs", "Fobs"),
    "null": ("Tnull", "Fnull"),
    "pvalue": ("pval", "pval"),
    "critical_value": ("qval", "qval"),
    "t": ("argvals", "argvals"),
}


@pytest.mark.parametrize(
    ("name", "field"),
    _params(
        {
            "tperm_fd_weather_atlantic_vs_pacific": ("Tobs", "Tnull", "pval", "qval", "argvals"),
            "fperm_fd_weather_temp_atlantic_dummy": ("Fobs", "Fnull", "pval", "qval", "argvals"),
        }
    ),
)
def test_permutation_tests(name: str, field: str) -> None:
    result = _tperm() if name.startswith("tperm") else _fperm()
    column = 0 if name.startswith("tperm") else 1
    attribute = next(attr for attr, pair in _PERM_FIELDS.items() if pair[column] == field)
    compare(getattr(result, attribute), CASES[name]["output"][field], _rtol(name))


def test_fperm_through_a_fitted_model() -> None:
    # SPEC 4.3 form f_test(fregress(...)): the golden Fperm.fd design (intercept +
    # Atlantic dummy, Fourier(25) betas, lambda = 1e2) fitted by fregress first.
    # It must give exactly the raw-form result above -- so the same golden fields
    # agree and the same three R defects apply -- and the R-exact fields match.
    case = CASES["fperm_fd_weather_temp_atlantic_dummy"]
    atlantic = np.asarray(case["input"]["atlantic"], dtype=float)
    model = fregress(
        weather_fd(),
        {"const": 1.0, "atlantic": atlantic},
        [WEATHER_BASIS, WEATHER_BASIS],
        lam=1e2,
        penalty=2,
    )
    result = f_test(
        model,
        n_perm=case["input"]["nperm"],
        q=case["input"]["q"],
        random_state=RRandom(case["input"]["seed"]),
    )
    raw = _fperm()
    assert result.statistic == raw.statistic
    assert result.pvalue == raw.pvalue
    assert result.critical_value == raw.critical_value
    np.testing.assert_array_equal(result.null, raw.null)
    np.testing.assert_array_equal(result.pointwise, raw.pointwise)
    compare(result.pvalue, case["output"]["pval"], _rtol(case["name"]))
    compare(result.t, case["output"]["argvals"], _rtol(case["name"]))
