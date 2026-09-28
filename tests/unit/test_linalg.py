"""Unit tests for fdatools._linalg (banded solvers, quadrature, Gram cache)."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fdatools import _linalg as la
from fdatools._backend import to_numpy

torch = pytest.importorskip("torch", reason="torch extra not installed")


def _spd_banded(n: int, bandwidth: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = np.zeros((n, n))
    for k in range(bandwidth + 1):
        vals = rng.uniform(0.1, 0.5, size=n - k)
        a += np.diag(vals, k)
        if k:
            a += np.diag(vals, -k)
    a += np.eye(n) * (bandwidth + 2.0)
    return a


# --------------------------------------------------------------------------- #
# banded storage
# --------------------------------------------------------------------------- #


def test_to_banded_upper_roundtrip() -> None:
    a = _spd_banded(8, 3)
    ab = la.to_banded(a, 3)
    assert ab.shape == (4, 8)
    for i in range(8):
        for j in range(8):
            if 0 <= j - i <= 3:
                assert ab[3 + i - j, j] == pytest.approx(a[i, j])


def test_bandwidth_of_detects_band() -> None:
    a = _spd_banded(10, 2)
    assert la.bandwidth_of(a) == 2
    assert la.bandwidth_of(np.eye(5)) == 0


# --------------------------------------------------------------------------- #
# banded / symmetric solves
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("bandwidth", [0, 1, 3])
def test_solve_spd_banded_numpy(bandwidth: int) -> None:
    a = _spd_banded(12, bandwidth)
    b = np.arange(12.0).reshape(12, 1) + 1.0
    x = la.solve_spd(a, b, bandwidth=bandwidth)
    np.testing.assert_allclose(a @ x, b, rtol=1e-12, atol=1e-12)


def test_solve_spd_1d_rhs() -> None:
    a = _spd_banded(6, 2)
    b = np.ones(6)
    x = la.solve_spd(a, b)
    assert x.shape == (6,)
    np.testing.assert_allclose(a @ x, b, rtol=1e-12, atol=1e-12)


def test_solve_spd_falls_back_when_not_positive_definite() -> None:
    a = np.array([[1.0, 2.0], [2.0, 1.0]])  # symmetric, indefinite
    b = np.array([[1.0], [0.0]])
    x = la.solve_spd(a, b)
    np.testing.assert_allclose(a @ x, b, rtol=1e-12, atol=1e-12)


def test_solve_spd_torch_matches_numpy() -> None:
    a = _spd_banded(9, 2)
    b = np.linspace(0.0, 1.0, 9).reshape(9, 1)
    want = la.solve_spd(a, b)
    got = la.solve_spd(torch.as_tensor(a), torch.as_tensor(b))
    assert isinstance(got, torch.Tensor)
    np.testing.assert_allclose(to_numpy(got), want, rtol=1e-10, atol=1e-12)


def test_solve_spd_torch_gradients_flow() -> None:
    a = torch.as_tensor(_spd_banded(5, 1))
    b = torch.ones(5, 1, dtype=torch.float64, requires_grad=True)
    x = la.solve_spd(a, b)
    x.sum().backward()
    assert b.grad is not None
    assert torch.all(torch.isfinite(b.grad))


def test_solve_spd_rejects_shape_mismatch() -> None:
    with pytest.raises(ValueError, match="shape"):
        la.solve_spd(np.eye(3), np.ones((4, 1)))


# --------------------------------------------------------------------------- #
# least squares
# --------------------------------------------------------------------------- #


def test_lstsq_recovers_exact_solution() -> None:
    rng = np.random.default_rng(3)
    a = rng.normal(size=(30, 5))
    coef = rng.normal(size=(5, 2))
    b = a @ coef
    got = la.lstsq(a, b)
    np.testing.assert_allclose(got, coef, rtol=1e-9, atol=1e-10)


def test_lstsq_torch_matches_numpy() -> None:
    rng = np.random.default_rng(4)
    a = rng.normal(size=(20, 4))
    b = rng.normal(size=(20, 3))
    want = la.lstsq(a, b)
    got = la.lstsq(torch.as_tensor(a), torch.as_tensor(b))
    np.testing.assert_allclose(to_numpy(got), want, rtol=1e-10, atol=1e-12)


def test_lstsq_torch_gradients_flow() -> None:
    a = torch.as_tensor(np.random.default_rng(5).normal(size=(12, 3)))
    b = torch.ones(12, 1, dtype=torch.float64, requires_grad=True)
    la.lstsq(a, b).sum().backward()
    assert b.grad is not None


# --------------------------------------------------------------------------- #
# sparse penalty assembly
# --------------------------------------------------------------------------- #


def test_sparse_penalty_matches_dense() -> None:
    a = _spd_banded(11, 2)
    sp = la.sparse_penalty(a)
    np.testing.assert_allclose(sp.toarray(), a, rtol=0.0, atol=1e-15)


def test_sparse_penalty_zeroes_outside_band() -> None:
    a = _spd_banded(9, 1)
    a[0, 8] = 1e-30  # numerically negligible off-band entry
    sp = la.sparse_penalty(a, bandwidth=1)
    assert sp.toarray()[0, 8] == 0.0


# --------------------------------------------------------------------------- #
# Gauss-Legendre quadrature
# --------------------------------------------------------------------------- #


def test_gauss_legendre_is_exact_for_polynomials() -> None:
    nodes, weights = la.gauss_legendre(4, 0.0, 2.0)
    # integral of t^7 over [0, 2] = 2^8 / 8 = 32
    assert float(np.sum(weights * nodes**7)) == pytest.approx(32.0, rel=1e-13)


def test_gauss_legendre_cached_nodes_are_readonly() -> None:
    x, w = la.gauss_legendre_reference(3)
    with pytest.raises(ValueError):
        x[0] = 0.0
    with pytest.raises(ValueError):
        w[0] = 0.0


def test_composite_gauss_legendre_exact_per_panel() -> None:
    breaks = np.array([0.0, 0.5, 1.0, 3.0])
    nodes, weights = la.composite_gauss_legendre(breaks, 3)
    assert float(np.sum(weights * nodes**5)) == pytest.approx(3.0**6 / 6.0, rel=1e-13)
    assert float(np.sum(weights)) == pytest.approx(3.0, rel=1e-14)


def test_composite_gauss_legendre_skips_degenerate_panels() -> None:
    breaks = np.array([0.0, 0.0, 1.0, 1.0])
    nodes, weights = la.composite_gauss_legendre(breaks, 2)
    assert nodes.shape == (2,)
    assert float(np.sum(weights)) == pytest.approx(1.0, rel=1e-14)


@given(deg=st.integers(min_value=1, max_value=12))
@settings(max_examples=12, deadline=None)
def test_gauss_legendre_weights_sum_to_interval_length(deg: int) -> None:
    _, w = la.gauss_legendre(deg, -1.5, 2.5)
    assert float(np.sum(w)) == pytest.approx(4.0, rel=1e-13)


# --------------------------------------------------------------------------- #
# Gram cache
# --------------------------------------------------------------------------- #


def test_gram_cache_returns_same_object_for_same_key() -> None:
    la.clear_gram_cache()
    calls = []

    def compute() -> np.ndarray:
        calls.append(1)
        return np.eye(3)

    first = la.cached_gram(("dummy", 1, 2.0), 0, compute)
    second = la.cached_gram(("dummy", 1, 2.0), 0, compute)
    assert first is second
    assert len(calls) == 1


def test_gram_cache_distinguishes_deriv() -> None:
    la.clear_gram_cache()
    a = la.cached_gram(("dummy",), 0, lambda: np.eye(2))
    b = la.cached_gram(("dummy",), 1, lambda: np.zeros((2, 2)))
    assert not np.array_equal(a, b)


def test_gram_cache_result_is_readonly() -> None:
    la.clear_gram_cache()
    g = la.cached_gram(("ro",), 0, lambda: np.eye(2))
    with pytest.raises(ValueError):
        g[0, 0] = 5.0
