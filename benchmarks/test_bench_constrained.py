"""Micro-benchmarks for evaluating constrained (monotone, positive) smooths.

Growth-sized problem: a latent ``W`` on 35 order-6 B-splines for 54 curves,
evaluated with its third derivative (the acceleration's rate) on 1001 points.

Run with ``pytest benchmarks -k constrained --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fdatools import BSpline, FData
from fdatools.smoothing import SmoothResult

N_CURVES = 54
N_BASIS = 35
POINTS = np.linspace(1.0, 18.0, 1001)


def result(constraint: str) -> SmoothResult:
    rng = np.random.default_rng(20260928)
    latent = FData(
        0.3 * rng.standard_normal((N_BASIS, N_CURVES)),
        BSpline(domain=(1.0, 18.0), n_basis=N_BASIS, order=6),
    )
    return SmoothResult(
        fd=latent,
        df=0.0,
        gcv=np.zeros(N_CURVES),
        sse=0.0,
        penalty_matrix=np.zeros((N_BASIS, N_BASIS)),
        lam=0.0,
        y2c_map=None,
        beta=np.vstack([np.zeros(N_CURVES), np.ones(N_CURVES)]),
        constraint=constraint,
    )


@pytest.mark.benchmark(group="constrained")
@pytest.mark.parametrize("deriv", [0, 3])
def test_bench_monotone_evaluation(benchmark: Any, deriv: int) -> None:
    fit = result("monotone")
    values = benchmark(fit, POINTS, deriv)
    assert np.asarray(values).shape == (POINTS.size, N_CURVES)


@pytest.mark.benchmark(group="constrained")
def test_bench_positive_third_derivative(benchmark: Any) -> None:
    fit = result("positive")
    values = benchmark(fit, POINTS, 3)
    assert np.asarray(values).shape == (POINTS.size, N_CURVES)
