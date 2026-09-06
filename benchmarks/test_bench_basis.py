"""Micro-benchmarks for the hot paths of the basis layer.

Evaluation and penalty construction dominate every smoothing call, so both are
measured on a spline of realistic size (365 daily points, 65 basis functions --
the shape of the canonical weather problem).

Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import numpy as np
import pytest

from fabel import BSpline, Fourier
from fabel._linalg import clear_gram_cache

N_POINTS = 365
N_BASIS = 65
DOMAIN = (0.0, 365.0)
GRID = np.linspace(DOMAIN[0], DOMAIN[1], N_POINTS)


@pytest.fixture(scope="module")
def spline() -> BSpline:
    return BSpline(domain=DOMAIN, n_basis=N_BASIS, order=4)


@pytest.fixture(scope="module")
def fourier() -> Fourier:
    return Fourier(domain=DOMAIN, n_basis=N_BASIS)


@pytest.mark.benchmark(group="eval")
@pytest.mark.parametrize("deriv", [0, 1, 2])
def test_bench_bspline_eval(benchmark, spline: BSpline, deriv: int) -> None:
    values = benchmark(spline, GRID, deriv)
    assert values.shape == (N_POINTS, N_BASIS)


@pytest.mark.benchmark(group="eval")
def test_bench_fourier_eval(benchmark, fourier: Fourier) -> None:
    values = benchmark(fourier, GRID, 0)
    assert values.shape == (N_POINTS, N_BASIS)


@pytest.mark.benchmark(group="penalty")
@pytest.mark.parametrize("order", [0, 1, 2])
def test_bench_bspline_penalty(benchmark, spline: BSpline, order: int) -> None:
    # Penalties are memoised on the basis value, so the cache must be dropped
    # every round or this would time a dictionary lookup, not the quadrature.
    def build() -> np.ndarray:
        clear_gram_cache()
        return spline.penalty(order)

    matrix = benchmark(build)
    assert matrix.shape == (N_BASIS, N_BASIS)


@pytest.mark.benchmark(group="penalty")
def test_bench_bspline_penalty_cached(benchmark, spline: BSpline) -> None:
    spline.penalty(2)
    matrix = benchmark(spline.penalty, 2)
    assert matrix.shape == (N_BASIS, N_BASIS)


@pytest.mark.benchmark(group="penalty")
def test_bench_fourier_penalty(benchmark, fourier: Fourier) -> None:
    def build() -> np.ndarray:
        clear_gram_cache()
        return fourier.penalty(2)

    matrix = benchmark(build)
    assert matrix.shape == (N_BASIS, N_BASIS)
