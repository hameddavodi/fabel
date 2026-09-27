"""Benchmarks for continuous and landmark registration.

The problem size mirrors the growth golden case: ten curves in a 35-function
order-6 spline basis on ``[1, 18]`` (a 351-point registration grid) and a
five-function cubic warp basis.

Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import numpy as np
import pytest

from fabel import BSpline, FData
from fabel.registration import RegistrationResult, register

DOMAIN = (1.0, 18.0)
N_CURVES = 10


@pytest.fixture(scope="module")
def curves() -> FData:
    basis = BSpline(domain=DOMAIN, n_basis=35, order=6)
    grid = np.linspace(*DOMAIN, 400)
    rng = np.random.default_rng(20260927)
    centres = 11.5 + rng.normal(scale=1.2, size=N_CURVES)
    values = np.stack(
        [80.0 + 4.5 * grid + 20.0 / (1.0 + np.exp(-(grid - c) / 0.8)) for c in centres], axis=1
    )
    return FData(np.linalg.lstsq(basis(grid), values, rcond=None)[0], basis)


@pytest.mark.benchmark(group="registration")
@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
def test_bench_register_continuous(benchmark, curves: FData, criterion: str) -> None:
    warp_basis = BSpline(domain=DOMAIN, n_basis=5)
    result: RegistrationResult = benchmark(
        register, curves, warp_basis=warp_basis, lam=1.0, criterion=criterion
    )
    assert result.registered.n_curves == N_CURVES


@pytest.mark.benchmark(group="registration")
def test_bench_register_landmarks(benchmark, curves: FData) -> None:
    marks = np.linspace(10.5, 12.5, N_CURVES)
    result: RegistrationResult = benchmark(register, curves, landmarks=marks)
    assert result.warp_inverse is not None
