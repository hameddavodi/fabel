"""Benchmarks for multivariate registration and for applying warps to new curves.

The problem mirrors the gait golden case: five curves with two variables (a
hip-like and a knee-like angle) in a 21-function Fourier basis on ``[0, 20]``
(a 211-point registration grid) and a five-function cubic warp basis.  The
multivariate criterion costs one criterion evaluation per variable on top of
the shared warp derivatives.

Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import BSpline, FData, Fourier
from fabel.registration import RegistrationResult, register

DOMAIN = (0.0, 20.0)
N_CURVES = 5


@pytest.fixture(scope="module")
def curves() -> FData:
    basis = Fourier(domain=DOMAIN, n_basis=21)
    grid = np.linspace(*DOMAIN, 400, endpoint=False)
    rng = np.random.default_rng(20260928)
    shifts = rng.normal(scale=0.6, size=N_CURVES)
    phase = [2 * np.pi * (grid - s) / 20.0 for s in shifts]
    hip = np.stack([30.0 * np.cos(p) + 5.0 * np.sin(2 * p) for p in phase], axis=1)
    knee = np.stack([35.0 - 30.0 * np.cos(2 * p) + 8.0 * np.sin(p) for p in phase], axis=1)
    both = np.concatenate([hip, knee], axis=1)
    coefs = np.linalg.lstsq(basis(grid), both, rcond=None)[0]
    return FData(np.stack([coefs[:, :N_CURVES], coefs[:, N_CURVES:]], axis=2), basis)


@pytest.mark.benchmark(group="registration-multivariate")
def test_bench_register_multivariate(benchmark: Any, curves: FData) -> None:
    wbasis = BSpline(domain=DOMAIN, n_basis=5)
    result: RegistrationResult = benchmark(register, curves, warp_basis=wbasis, lam=1e-2)
    assert result.registered.coefs.shape == curves.coefs.shape


@pytest.mark.benchmark(group="registration-multivariate")
def test_bench_apply_warps(benchmark: Any, curves: FData) -> None:
    result = register(curves, landmarks=np.linspace(9.5, 10.5, N_CURVES))
    velocity = curves.derivative()
    warped: FData = benchmark(result.apply, velocity)
    assert warped.coefs.shape == velocity.coefs.shape
