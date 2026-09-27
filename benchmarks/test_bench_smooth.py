"""Micro-benchmarks for the hot paths of :mod:`fabel.smoothing`.

The shapes follow SPEC §6: the canonical weather problem (365 daily points,
65 basis functions, 35 curves), the 1000-curve smoothing target, and a GCV
grid search, which one eigendecomposition serves for every lambda.

Run with ``pytest benchmarks -k smooth --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import BSpline, Fourier
from fabel.smoothing import Smoother, gcv_curve, smooth

N_POINTS = 365
N_BASIS = 65
DOMAIN = (0.0, 365.0)
GRID = np.linspace(DOMAIN[0], DOMAIN[1], N_POINTS)
RNG = np.random.default_rng(20260927)


def curves(n_curves: int) -> np.ndarray:
    """Return ``n_curves`` noisy annual cycles sampled on the daily grid."""
    phase = RNG.uniform(0.0, 2 * np.pi, n_curves)
    signal = 10.0 * np.sin(2 * np.pi * GRID[:, None] / 365.0 + phase[None, :])
    return signal + RNG.standard_normal((N_POINTS, n_curves))


@pytest.fixture(scope="module")
def spline() -> BSpline:
    return BSpline(domain=DOMAIN, n_basis=N_BASIS, order=4)


@pytest.fixture(scope="module")
def fourier() -> Fourier:
    return Fourier(domain=DOMAIN, n_basis=N_BASIS)


@pytest.mark.benchmark(group="smooth")
def test_bench_smooth_weather_fixed_lambda(benchmark: Any, fourier: Fourier) -> None:
    y = curves(35)
    result = benchmark(smooth, y, GRID, basis=fourier, lam=1e2)
    assert np.asarray(result.fd.coefs).shape == (N_BASIS, 35)


@pytest.mark.benchmark(group="smooth")
def test_bench_smooth_1000_curves(benchmark: Any, spline: BSpline) -> None:
    y = curves(1000)
    result = benchmark(smooth, y, GRID, basis=spline, lam=1e2)
    assert np.asarray(result.fd.coefs).shape == (N_BASIS, 1000)


@pytest.mark.benchmark(group="smooth")
def test_bench_smooth_gcv_search(benchmark: Any, spline: BSpline) -> None:
    y = curves(35)
    result = benchmark(smooth, y, GRID, basis=spline, lam="gcv")
    assert result.lam > 0.0


@pytest.mark.benchmark(group="smooth")
def test_bench_smooth_gcv_curve_grid(benchmark: Any, spline: BSpline) -> None:
    y = curves(35)
    lambdas = 10.0 ** np.arange(-4.0, 8.25, 0.25)
    scores = benchmark(gcv_curve, y, GRID, spline, lambdas)
    assert scores.shape == (lambdas.size, 35)


@pytest.mark.benchmark(group="smooth")
def test_bench_smooth_monotone(benchmark: Any) -> None:
    t = np.linspace(0.0, 1.0, 101)
    y = np.tanh(6.0 * (t - 0.5)) + 0.02 * RNG.standard_normal(t.size)
    basis = BSpline(domain=(0.0, 1.0), n_basis=12)
    result = benchmark(smooth, y, t, basis=basis, lam=1e-4, constraint="monotone")
    assert np.all(np.diff(np.asarray(result(t))[:, 0]) >= -1e-12)


@pytest.mark.benchmark(group="smooth")
def test_bench_smoother_transform(benchmark: Any, spline: BSpline) -> None:
    x = curves(200).T
    smoother = Smoother(spline, t=GRID, lam=1e2).fit(x)
    coefs = benchmark(smoother.transform, x)
    assert coefs.shape == (200, N_BASIS)
