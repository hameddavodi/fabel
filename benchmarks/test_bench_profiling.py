"""Benchmarks for generalized profiling (``fabel.profiling``).

The CSTR problem mirrors the golden case: 97 temperature and concentration
observations on ``[0, 24]``, a 51-function cubic B-spline basis per state
(240 Simpson nodes), ``kref`` and ``EoverR`` estimated.  The FitzHugh-Nagumo
problem has 201 observations of both states and a 203-function basis.

Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import BSpline
from fabel.profiling import CSTR_PARAMETERS, ProfiledODE, cstr_model, fitzhugh_nagumo_model


@pytest.fixture(scope="module")
def cstr_problem() -> ProfiledODE:
    t = np.arange(0.0, 24.001, 0.25)
    x = cstr_model().simulate(
        t, [1.5965, 341.3754], list(CSTR_PARAMETERS.values()), breaks=range(4, 24, 4)
    )
    y = x + np.random.default_rng(0).standard_normal(x.shape) * [0.02, 0.5]
    wt = y.var(axis=0, ddof=1)
    return ProfiledODE(
        cstr_model(estimate=("kref", "EoverR")),
        t,
        y,
        BSpline(domain=(0.0, 24.0), breaks=np.arange(0.0, 24.001, 0.5).tolist()),
        lam=[100.0, 100.0],
        state_weights=list(1.0 / wt),
    )


@pytest.fixture(scope="module")
def fhn_problem() -> ProfiledODE:
    t = np.linspace(0.0, 20.0, 201)
    fhn = fitzhugh_nagumo_model()
    x = fhn.simulate(t, [-1.0, 1.0], [0.2, 0.2, 3.0])
    y = x + 0.05 * np.random.default_rng(1).standard_normal(x.shape)
    basis = BSpline(domain=(0.0, 20.0), breaks=np.linspace(0.0, 20.0, 201).tolist())
    return ProfiledODE(fhn, t, y, basis, lam=1e3)


def test_bench_cstr_inner(benchmark: Any, cstr_problem: ProfiledODE) -> None:
    inner = benchmark(cstr_problem.fit_states, [0.461, 0.83301])
    assert inner.converged


def test_bench_cstr_profile(benchmark: Any, cstr_problem: ProfiledODE) -> None:
    fit = benchmark(cstr_problem.fit, [0.4, 0.8])
    assert fit.converged


def test_bench_fhn_profile(benchmark: Any, fhn_problem: ProfiledODE) -> None:
    fit = benchmark(fhn_problem.fit, [0.3, 0.3, 2.5])
    assert fit.converged
