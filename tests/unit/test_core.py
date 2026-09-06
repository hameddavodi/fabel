"""Unit tests for :mod:`fabel.core` (FData, BiFData, inner products)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fabel import LDO, BiFData, BSpline, Constant, Exponential, FData, Fourier, Monomial, inprod
from fabel import _linalg as la
from fabel._backend import to_numpy

torch = pytest.importorskip("torch", reason="torch extra not installed")

RNG = np.random.default_rng(20260906)


@pytest.fixture(autouse=True)
def _clear_cache() -> None:
    la.clear_gram_cache()


def make_fd(n_curves: int = 3, n_basis: int = 8, seed: int = 0) -> FData:
    """Return a random B-spline FData for tests."""
    basis = BSpline(domain=(0.0, 1.0), n_basis=n_basis)
    coefs = np.random.default_rng(seed).normal(size=(n_basis, n_curves))
    return FData(coefs, basis)


# --------------------------------------------------------------------------- #
# construction
# --------------------------------------------------------------------------- #


def test_one_dimensional_coefs_become_a_single_curve() -> None:
    fd = FData(np.arange(8.0), BSpline(domain=(0.0, 1.0), n_basis=8))
    assert fd.coefs.shape == (8, 1)
    assert fd.n_curves == 1
    assert len(fd) == 1


def test_rejects_coefs_that_do_not_match_the_basis() -> None:
    with pytest.raises(ValueError, match="8 rows"):
        FData(np.ones((5, 2)), BSpline(domain=(0.0, 1.0), n_basis=8))


def test_rejects_more_than_three_axes() -> None:
    with pytest.raises(ValueError, match="at most three"):
        FData(np.ones((8, 2, 2, 2)), BSpline(domain=(0.0, 1.0), n_basis=8))


def test_domain_and_repr_come_from_the_basis() -> None:
    fd = make_fd()
    assert fd.domain == (0.0, 1.0)
    assert "FData" in repr(fd)
    assert "n_curves=3" in repr(fd)


# --------------------------------------------------------------------------- #
# evaluation
# --------------------------------------------------------------------------- #


def test_evaluation_is_the_basis_matrix_times_the_coefficients() -> None:
    fd = make_fd()
    t = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose(fd(t), fd.basis(t) @ fd.coefs, rtol=0.0, atol=1e-15)


def test_evaluation_keeps_the_variable_axis() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    fd = FData(RNG.normal(size=(6, 4, 2)), basis)
    values = fd(np.linspace(0.0, 1.0, 9))
    assert values.shape == (9, 4, 2)
    np.testing.assert_allclose(values[..., 0], FData(fd.coefs[..., 0], basis)(np.linspace(0, 1, 9)))


def test_evaluation_accepts_an_ldo() -> None:
    fd = make_fd()
    t = np.linspace(0.0, 1.0, 7)
    want = 2.0 * fd(t) - 3.0 * fd(t, deriv=1) + fd(t, deriv=2)
    np.testing.assert_allclose(fd(t, LDO(weights=[2.0, -3.0])), want, rtol=1e-12, atol=1e-14)


def test_evaluation_in_torch_matches_numpy_and_flows_gradients() -> None:
    fd = make_fd()
    t = torch.linspace(0.05, 0.95, 13, dtype=torch.float64, requires_grad=True)
    values = fd(t)
    assert isinstance(values, torch.Tensor)
    np.testing.assert_allclose(to_numpy(values), fd(to_numpy(t)), rtol=1e-10, atol=1e-12)
    values.sum().backward()
    assert t.grad is not None
    assert bool(torch.all(torch.isfinite(t.grad)))


def test_to_numpy_and_to_torch_round_trip() -> None:
    fd = make_fd()
    t = np.linspace(0.0, 1.0, 5)
    np.testing.assert_allclose(to_numpy(fd.to_torch(t)), fd.to_numpy(t), rtol=1e-12, atol=1e-14)


# --------------------------------------------------------------------------- #
# derivatives
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "basis",
    [
        BSpline(domain=(0.0, 1.0), n_basis=8, order=5),
        Fourier(domain=(0.0, 1.0), n_basis=7),
        Monomial(domain=(0.0, 2.0), n_basis=5),
        Exponential(domain=(0.0, 1.0), rates=[0.0, 0.5, 1.5]),
        Constant(domain=(0.0, 1.0)),
    ],
    ids=lambda b: type(b).__name__,
)
@pytest.mark.parametrize("order", [1, 2])
def test_derivative_reproduces_pointwise_derivative(basis: Any, order: int) -> None:
    fd = FData(RNG.normal(size=(basis.n_basis, 2)), basis)
    t = np.linspace(0.05, 0.95, 41) * (basis.domain[1] - basis.domain[0]) + basis.domain[0]
    np.testing.assert_allclose(fd.derivative(order)(t), fd(t, deriv=order), rtol=1e-9, atol=1e-10)


def test_derivative_of_a_spline_lowers_the_order() -> None:
    fd = make_fd(n_basis=9)
    d = fd.derivative()
    assert isinstance(d.basis, BSpline)
    assert d.basis.order == fd.basis.order - 1
    assert d.basis.breaks == fd.basis.breaks


def test_derivative_beyond_the_order_is_zero() -> None:
    fd = FData(np.ones((4, 1)), BSpline(domain=(0.0, 1.0), n_basis=4, order=4))
    d = fd.derivative(4)
    t = np.linspace(0.0, 1.0, 7)
    np.testing.assert_allclose(d(t), np.zeros((7, 1)), rtol=0.0, atol=1e-12)


def test_derivative_rejects_a_negative_order() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        make_fd().derivative(-1)


def test_derivative_of_zero_returns_the_same_function() -> None:
    fd = make_fd()
    t = np.linspace(0.0, 1.0, 5)
    np.testing.assert_allclose(fd.derivative(0)(t), fd(t), rtol=0.0, atol=1e-15)


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #


def test_mean_is_the_pointwise_average() -> None:
    fd = make_fd(n_curves=6)
    t = np.linspace(0.0, 1.0, 9)
    np.testing.assert_allclose(fd.mean()(t), fd(t).mean(axis=1, keepdims=True), atol=1e-14)


def test_center_removes_the_mean() -> None:
    fd = make_fd(n_curves=6)
    t = np.linspace(0.0, 1.0, 9)
    np.testing.assert_allclose(fd.center()(t).mean(axis=1), np.zeros(9), atol=1e-14)


def test_std_is_exact_when_the_deviation_lies_in_the_basis() -> None:
    # x_j = mean + s_j * g with g > 0, so sd(t) = sd(s) * g(t) is in the span
    # and the least-squares projection reproduces it exactly.
    basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    shape = np.linspace(1.0, 3.0, 6)[:, None]
    scales = RNG.normal(size=(1, 9))
    fd = FData(np.ones((6, 9)) + shape * scales, basis)
    t = np.linspace(0.0, 1.0, 60)
    want = fd(t).std(axis=1, ddof=1, keepdims=True)
    np.testing.assert_allclose(fd.std()(t), want, rtol=1e-9, atol=1e-11)


def test_std_tracks_the_pointwise_sample_standard_deviation() -> None:
    fd = make_fd(n_curves=25, n_basis=6)
    t = np.linspace(0.02, 0.98, 60)
    want = fd(t).std(axis=1, ddof=1, keepdims=True)
    np.testing.assert_allclose(fd.std()(t), want, rtol=0.2, atol=0.02)


def test_std_needs_at_least_two_curves() -> None:
    with pytest.raises(ValueError, match="at least two curves"):
        make_fd(n_curves=1).std()


def test_cov_is_the_pointwise_covariance_surface() -> None:
    fd = make_fd(n_curves=12, n_basis=6)
    s = np.linspace(0.0, 1.0, 7)
    surface = fd.cov()(s, s)
    values = fd(s)
    want = np.cov(values, ddof=1)
    assert isinstance(fd.cov(), BiFData)
    np.testing.assert_allclose(surface, want, rtol=1e-9, atol=1e-12)


# --------------------------------------------------------------------------- #
# arithmetic
# --------------------------------------------------------------------------- #


def test_addition_and_subtraction_are_pointwise() -> None:
    a, b = make_fd(seed=1), make_fd(seed=2)
    t = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose((a + b)(t), a(t) + b(t), atol=1e-14)
    np.testing.assert_allclose((a - b)(t), a(t) - b(t), atol=1e-14)
    np.testing.assert_allclose((-a)(t), -a(t), atol=1e-14)


def test_scalar_arithmetic_is_pointwise() -> None:
    a = make_fd(seed=3)
    t = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose((a * 2.5)(t), 2.5 * a(t), atol=1e-14)
    np.testing.assert_allclose((2.5 * a)(t), 2.5 * a(t), atol=1e-14)
    np.testing.assert_allclose((a / 4.0)(t), a(t) / 4.0, atol=1e-14)
    np.testing.assert_allclose((a + 1.5)(t), a(t) + 1.5, atol=1e-12)
    np.testing.assert_allclose((1.5 - a)(t), 1.5 - a(t), atol=1e-12)


def test_product_of_two_functions_is_exact() -> None:
    a, b = make_fd(n_curves=2, seed=4), make_fd(n_curves=2, seed=5)
    t = np.linspace(0.0, 1.0, 101)
    np.testing.assert_allclose((a * b)(t), a(t) * b(t), rtol=1e-9, atol=1e-11)


def test_integer_power_is_exact() -> None:
    a = make_fd(n_curves=2, n_basis=6, seed=6)
    t = np.linspace(0.0, 1.0, 101)
    np.testing.assert_allclose((a**3)(t), a(t) ** 3, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose((a**0)(t), np.ones_like(a(t)), rtol=1e-10, atol=1e-12)


def test_fractional_power_is_approximated_on_a_positive_function() -> None:
    a = make_fd(n_curves=1, n_basis=5, seed=7) ** 2 + 1.0
    t = np.linspace(0.0, 1.0, 101)
    np.testing.assert_allclose((a**0.5)(t), a(t) ** 0.5, rtol=1e-4, atol=1e-5)


def test_division_by_a_function_is_rejected() -> None:
    a = make_fd()
    with pytest.raises(TypeError, match="not closed"):
        _ = a / a


def test_arithmetic_between_different_bases_is_rejected() -> None:
    a = make_fd()
    b = FData(np.ones((6, 3)), BSpline(domain=(0.0, 1.0), n_basis=6))
    with pytest.raises(ValueError, match="same basis"):
        _ = a + b


def test_arithmetic_between_different_curve_counts_is_rejected() -> None:
    with pytest.raises(ValueError, match="curves"):
        _ = make_fd(n_curves=3) + make_fd(n_curves=2)


# --------------------------------------------------------------------------- #
# inner products
# --------------------------------------------------------------------------- #


def test_inner_product_matrix_matches_quadrature() -> None:
    a, b = make_fd(n_curves=3, seed=8), make_fd(n_curves=4, seed=9)
    nodes, weights = la.composite_gauss_legendre(np.linspace(0.0, 1.0, 41), 12)
    want = a(nodes).T @ (weights[:, None] * b(nodes))
    np.testing.assert_allclose(a @ b, want, rtol=1e-9, atol=1e-12)
    assert (a @ b).shape == (3, 4)


def test_inner_product_honours_the_operators() -> None:
    a, b = make_fd(seed=10), make_fd(seed=11)
    nodes, weights = la.composite_gauss_legendre(np.linspace(0.0, 1.0, 41), 12)
    want = a(nodes, deriv=1).T @ (weights[:, None] * b(nodes, deriv=2))
    np.testing.assert_allclose(inprod(a, b, lfd1=1, lfd2=2), want, rtol=1e-8, atol=1e-11)


def test_inner_product_of_two_bases_is_the_cross_gram() -> None:
    left, right = BSpline(domain=(0.0, 1.0), n_basis=6), Fourier(domain=(0.0, 1.0), n_basis=5)
    got = inprod(left, right)
    assert got.shape == (6, 5)
    nodes, weights = la.composite_gauss_legendre(np.linspace(0.0, 1.0, 81), 12)
    want = left(nodes).T @ (weights[:, None] * right(nodes))
    np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-12)


def test_inner_product_of_a_basis_with_itself_is_its_gram() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=7)
    np.testing.assert_allclose(inprod(basis, basis), basis.gram(), rtol=1e-10, atol=1e-13)


def test_inner_product_rejects_different_domains() -> None:
    a = make_fd()
    b = FData(np.ones((8, 2)), BSpline(domain=(0.0, 2.0), n_basis=8))
    with pytest.raises(ValueError, match="domain"):
        _ = a @ b


# --------------------------------------------------------------------------- #
# indexing
# --------------------------------------------------------------------------- #


def test_indexing_selects_curves() -> None:
    fd = make_fd(n_curves=5)
    t = np.linspace(0.0, 1.0, 6)
    np.testing.assert_allclose(fd[2](t), fd(t)[:, [2]], atol=1e-15)
    np.testing.assert_allclose(fd[1:4](t), fd(t)[:, 1:4], atol=1e-15)
    np.testing.assert_allclose(fd[[0, 3]](t), fd(t)[:, [0, 3]], atol=1e-15)
    assert len(fd[1:4]) == 3


def test_iteration_yields_single_curves() -> None:
    fd = make_fd(n_curves=3)
    curves = list(fd)
    assert len(curves) == 3
    assert all(len(c) == 1 for c in curves)


# --------------------------------------------------------------------------- #
# BiFData
# --------------------------------------------------------------------------- #


def test_bifd_evaluates_the_tensor_product() -> None:
    sbasis = BSpline(domain=(0.0, 1.0), n_basis=5)
    tbasis = Fourier(domain=(0.0, 1.0), n_basis=3)
    coefs = RNG.normal(size=(5, 3))
    bifd = BiFData(coefs, sbasis, tbasis)
    s, t = np.linspace(0.0, 1.0, 4), np.linspace(0.0, 1.0, 6)
    want = sbasis(s) @ coefs @ tbasis(t).T
    np.testing.assert_allclose(bifd(s, t), want, rtol=1e-12, atol=1e-14)


def test_bifd_keeps_trailing_axes() -> None:
    sbasis = BSpline(domain=(0.0, 1.0), n_basis=4)
    tbasis = BSpline(domain=(0.0, 1.0), n_basis=3, order=3)
    bifd = BiFData(RNG.normal(size=(4, 3, 2)), sbasis, tbasis)
    assert bifd(np.linspace(0, 1, 5), np.linspace(0, 1, 7)).shape == (5, 7, 2)


def test_bifd_rejects_mismatched_coefficients() -> None:
    with pytest.raises(ValueError, match="shape"):
        BiFData(np.ones((3, 3)), BSpline(domain=(0.0, 1.0), n_basis=4), Constant())


def test_bifd_transpose_swaps_the_arguments() -> None:
    sbasis = BSpline(domain=(0.0, 1.0), n_basis=5)
    tbasis = BSpline(domain=(0.0, 1.0), n_basis=4, order=3)
    bifd = BiFData(RNG.normal(size=(5, 4)), sbasis, tbasis)
    s, t = np.linspace(0.0, 1.0, 3), np.linspace(0.0, 1.0, 6)
    np.testing.assert_allclose(bifd.transpose()(t, s), bifd(s, t).T, rtol=1e-12, atol=1e-14)


# --------------------------------------------------------------------------- #
# property tests
# --------------------------------------------------------------------------- #


@given(scale=st.floats(min_value=-5.0, max_value=5.0), shift=st.floats(-3.0, 3.0))
@settings(max_examples=25, deadline=None)
def test_linearity_of_evaluation(scale: float, shift: float) -> None:
    a, b = make_fd(seed=12), make_fd(seed=13)
    t = np.linspace(0.0, 1.0, 9)
    np.testing.assert_allclose(
        (a * scale + b * shift)(t), scale * a(t) + shift * b(t), rtol=1e-9, atol=1e-12
    )


@given(n_curves=st.integers(min_value=2, max_value=6))
@settings(max_examples=8, deadline=None)
def test_inner_product_is_symmetric(n_curves: int) -> None:
    fd = make_fd(n_curves=n_curves, seed=14)
    gram = fd @ fd
    np.testing.assert_allclose(gram, gram.T, rtol=1e-12, atol=1e-14)
    assert float(np.linalg.eigvalsh(gram).min()) > -1e-10
