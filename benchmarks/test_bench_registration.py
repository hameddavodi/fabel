"""Benchmarks for continuous and landmark registration.

The problem size mirrors the growth golden case: ten curves in a 35-function
order-6 spline basis on ``[1, 18]`` (a 351-point registration grid) and a
five-function cubic warp basis.

The continuous benchmarks run twice: ``numpy`` (analytic gradient and Hessian)
and ``torch`` (the same curves as tensors, so :func:`fdatools.registration.register`
takes the autodiff path: criterion in PyTorch, gradient and Hessian by
autograd, the same Newton iteration).  The torch cases skip without the
``fdatools[torch]`` extra.

Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import importlib.util
from typing import Any

import numpy as np
import pytest

from fdatools import BSpline, FData
from fdatools.registration import RegistrationResult, register

HAS_TORCH = importlib.util.find_spec("torch") is not None

BACKENDS = [
    "numpy",
    pytest.param("torch", marks=pytest.mark.skipif(not HAS_TORCH, reason="torch not installed")),
]

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


def _on_backend(fd: FData, backend: str) -> FData:
    """Return ``fd`` with NumPy or torch coefficients."""
    if backend == "numpy":
        return fd
    import torch

    return FData(torch.tensor(np.asarray(fd.coefs), dtype=torch.float64), fd.basis)


@pytest.mark.benchmark(group="registration")
@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
def test_bench_register_continuous(
    benchmark: Any, curves: FData, criterion: str, backend: str
) -> None:
    warp_basis = BSpline(domain=DOMAIN, n_basis=5)
    data = _on_backend(curves, backend)
    result: RegistrationResult = benchmark(
        register, data, warp_basis=warp_basis, lam=1.0, criterion=criterion
    )
    assert result.registered.n_curves == N_CURVES


@pytest.mark.benchmark(group="registration")
def test_bench_register_landmarks(benchmark: Any, curves: FData) -> None:
    marks = np.linspace(10.5, 12.5, N_CURVES)
    result: RegistrationResult = benchmark(register, curves, landmarks=marks)
    assert result.warp_inverse is not None
