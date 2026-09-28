"""Unit tests for fdatools.basis (all seven basis systems)."""

from __future__ import annotations

from math import pi, sqrt
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fdatools import (
    LDO,
    Basis,
    BSpline,
    Constant,
    Exponential,
    Fourier,
    Monomial,
    Polygonal,
    Power,
)
from fdatools import _linalg as la
from fdatools._backend import to_numpy

torch = pytest.importorskip("torch", reason="torch extra not installed")


def all_bases() -> list[tuple[Basis, float, float]]:
    """Return one instance of every basis with a safe evaluation interval."""
    return [
        (BSpline(domain=(0.0, 1.0), n_basis=7, order=4), 0.0, 1.0),
        (BSpline(domain=(-1.0, 2.0), breaks=[-1.0, 0.0, 0.5, 2.0], order=3), -1.0, 2.0),
        (Fourier(domain=(0.0, 1.0), n_basis=5), 0.0, 1.0),
        (Fourier(domain=(0.0, 365.0), n_basis=7, period=365.0), 0.0, 365.0),
        (Monomial(domain=(0.0, 1.0), n_basis=4), 0.0, 1.0),
        (Power(domain=(0.5, 2.0), exponents=[0.0, 0.5, 1.0, 2.0]), 0.5, 2.0),
        (Exponential(domain=(0.0, 1.0), rates=[0.0, 1.0, -2.0]), 0.0, 1.0),
        (Constant(domain=(0.0, 3.0)), 0.0, 3.0),
        (Polygonal([0.0, 0.25, 0.6, 1.0]), 0.0, 1.0),
    ]


@pytest.fixture(autouse=True)
def _clear_cache() -> None:
    la.clear_gram_cache()


# --------------------------------------------------------------------------- #
# shared behaviour
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_evaluation_shape_and_namespace(basis: Basis, lo: float, hi: float) -> None:
    t = np.linspace(lo, hi, 11)
    mat = basis(t)
    assert mat.shape == (11, basis.n_basis)
    assert isinstance(mat, np.ndarray)
    assert len(basis.names) == basis.n_basis


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_torch_matches_numpy(basis: Basis, lo: float, hi: float) -> None:
    grid = np.linspace(lo, hi, 9)
    for deriv in (0, 1, 2):
        want = basis(grid, deriv=deriv)
        got = basis(torch.as_tensor(grid), deriv=deriv)
        assert isinstance(got, torch.Tensor)
        np.testing.assert_allclose(to_numpy(got), want, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_gradients_flow_through_evaluation(basis: Basis, lo: float, hi: float) -> None:
    t = torch.linspace(lo, hi, 8, dtype=torch.float64, requires_grad=True)
    mat = basis(t)
    if not mat.requires_grad:  # the constant basis genuinely does not depend on t
        assert isinstance(basis, Constant)
        return
    mat.sum().backward()
    assert t.grad is not None
    assert torch.all(torch.isfinite(t.grad))


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_penalty_is_symmetric_and_positive_semidefinite(basis: Basis, lo: float, hi: float) -> None:
    for op in (0, 1, 2):
        pen = basis.penalty(op)
        assert pen.shape == (basis.n_basis, basis.n_basis)
        np.testing.assert_allclose(pen, pen.T, rtol=0.0, atol=1e-12)
        smallest = float(np.linalg.eigvalsh(pen).min())
        assert smallest > -1e-9 * max(1.0, float(np.abs(pen).max()))


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_gram_matches_numerical_quadrature(basis: Basis, lo: float, hi: float) -> None:
    panels = np.asarray(sorted(set(basis._natural_breaks()) | {lo, hi}))
    nodes, weights = la.composite_gauss_legendre(panels, 30)
    mat = basis(nodes)
    want = mat.T @ (weights[:, None] * mat)
    np.testing.assert_allclose(basis.gram(), want, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_frozen_hashable_and_equal_by_value(basis: Basis, lo: float, hi: float) -> None:
    twin = type(basis)(**_init_kwargs(basis))
    assert basis == twin
    assert hash(basis) == hash(twin)
    with pytest.raises(AttributeError):
        basis.domain = (0.0, 1.0)  # type: ignore[misc]


def _init_fields(basis: Basis) -> list[Any]:
    from dataclasses import fields

    skip = {"Polygonal": {"domain"}, "Fourier": {"n_harmonics"}}
    return [f for f in fields(basis) if f.name not in skip.get(type(basis).__name__, set())]


def _init_kwargs(basis: Basis) -> dict[str, Any]:
    kwargs = {f.name: getattr(basis, f.name) for f in _init_fields(basis)}
    if isinstance(basis, Fourier):
        kwargs["n_basis"] = basis.n_basis
    return kwargs


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_ldo_matches_the_weighted_derivative_sum(basis: Basis, lo: float, hi: float) -> None:
    t = np.linspace(lo, hi, 7)
    op = LDO(weights=[2.0, -3.0])
    want = 2.0 * basis(t, deriv=0) - 3.0 * basis(t, deriv=1) + basis(t, deriv=2)
    np.testing.assert_allclose(basis(t, op), want, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_negative_derivative_is_rejected(basis: Basis, lo: float, hi: float) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        basis(np.array([lo]), deriv=-1)


@pytest.mark.parametrize(("basis", "lo", "hi"), all_bases(), ids=lambda v: str(v)[:40])
def test_penalty_is_cached_and_read_only(basis: Basis, lo: float, hi: float) -> None:
    first = basis.penalty(2)
    assert basis.penalty(2) is first
    with pytest.raises(ValueError):
        first[0, 0] = 1.0


# --------------------------------------------------------------------------- #
# B-spline
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("order", [1, 2, 3, 4, 6])
def test_bspline_is_a_partition_of_unity(order: int) -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=order + 5, order=order)
    mat = basis(np.linspace(0.0, 1.0, 101))
    np.testing.assert_allclose(mat.sum(axis=1), np.ones(101), rtol=0.0, atol=1e-14)


@pytest.mark.parametrize("order", [2, 3, 4, 6])
def test_bspline_derivatives_sum_to_zero(order: int) -> None:
    basis = BSpline(domain=(0.0, 2.0), n_basis=order + 4, order=order)
    mat = basis(np.linspace(0.0, 2.0, 51), deriv=1)
    np.testing.assert_allclose(mat.sum(axis=1), np.zeros(51), rtol=0.0, atol=1e-11)


def test_bspline_is_nonnegative() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=4)
    assert float(basis(np.linspace(0.0, 1.0, 201)).min()) >= 0.0


def test_bspline_right_endpoint_is_included() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=6, order=4)
    row = basis(np.array([1.0]))[0]
    assert float(row[-1]) == pytest.approx(1.0)


def test_bspline_top_derivative_vanishes_at_the_right_endpoint() -> None:
    # D^(k-1) of an order-k spline is a step function with no right-hand limit
    # at the end of the domain, so it is reported as zero (matching R's fda).
    basis = BSpline(domain=(0.0, 1.0), n_basis=5, order=2)
    np.testing.assert_allclose(basis(np.array([1.0]), deriv=1), np.zeros((1, 5)))
    assert float(np.abs(basis(np.array([0.99]), deriv=1)).max()) > 0.0


def test_bspline_derivative_beyond_the_order_is_zero() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=6, order=4)
    np.testing.assert_allclose(basis(np.linspace(0.0, 1.0, 5), deriv=4), np.zeros((5, 6)))


def test_bspline_derivative_matches_finite_differences() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=8, order=5)
    t = np.linspace(0.15, 0.85, 20)
    h = 1e-6
    want = (basis(t + h) - basis(t - h)) / (2.0 * h)
    np.testing.assert_allclose(basis(t, deriv=1), want, rtol=1e-6, atol=1e-7)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_bspline_derivative_map_handles_repeated_interior_knots(n: int) -> None:
    """D^n of a spline with multiple knots is the a.e. derivative, off the breaks."""
    basis = BSpline(domain=(0.0, 1.0), order=7, breaks=[0.0, *([0.4] * 4), *([0.7] * 2), 1.0])
    derived, matrix = basis._derivative_map(n)
    t = np.linspace(0.02, 0.98, 61)
    t = t[np.min(np.abs(t[:, None] - np.array([0.4, 0.7])[None, :]), axis=1) > 1e-3]
    np.testing.assert_allclose(derived(t) @ matrix, basis(t, deriv=n), rtol=1e-9, atol=1e-9)


def test_bspline_derivative_map_caps_multiplicity_at_the_reduced_order() -> None:
    basis = BSpline(domain=(0.0, 1.0), order=7, breaks=[0.0, *([0.4] * 4), 1.0])

    # Multiplicity 4 survives while it fits the reduced order, then saturates.
    def multiplicity(n: int) -> int:
        derived = basis._derivative_map(n)[0]
        assert isinstance(derived, BSpline)
        return derived.breaks.count(0.4)

    assert multiplicity(2) == 4
    assert multiplicity(4) == 3
    assert multiplicity(6) == 1


def test_bspline_discontinuous_spline_can_be_differentiated() -> None:
    basis = BSpline(domain=(0.0, 1.0), order=2, breaks=[0.0, 0.5, 0.5, 1.0])
    derived, matrix = basis._derivative_map(1)
    t = np.array([0.1, 0.3, 0.7, 0.9])
    np.testing.assert_allclose(derived(t) @ matrix, basis(t, deriv=1), atol=1e-12)


def test_bspline_penalty_is_exactly_banded() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=12, order=4)
    pen = basis.penalty(2)
    index = np.arange(12)
    outside = np.abs(index[:, None] - index[None, :]) >= basis.order
    assert not np.any(pen[outside])


def test_bspline_breaks_and_n_basis_agree() -> None:
    basis = BSpline(domain=(0.0, 1.0), breaks=[0.0, 0.4, 1.0], order=4)
    assert basis.n_basis == 5
    assert len(basis.knots) == basis.n_basis + basis.order


def test_bspline_defaults_to_equally_spaced_breaks() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=6, order=3)
    np.testing.assert_allclose(basis.breaks, np.linspace(0.0, 1.0, 5))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"order": 0}, "order must be at least 1"),
        ({"n_basis": 2, "order": 4}, "n_basis must be at least order"),
        ({"breaks": [0.0], "order": 2}, "at least the two endpoints"),
        ({"breaks": [0.0, 0.6, 0.5, 1.0]}, "non-decreasing"),
        ({"breaks": [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0]}, "exceeds the order"),
        ({"breaks": [0.0, 0.0, 1.0]}, "end breaks must not repeat"),
        ({"breaks": [0.1, 0.5, 1.0]}, "must span the domain"),
        ({"breaks": [0.0, 0.5, 1.0], "n_basis": 99}, "conflicts with"),
        ({"domain": (1.0, 0.0)}, "must be increasing"),
    ],
)
def test_bspline_rejects_bad_arguments(kwargs: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        BSpline(**{"domain": (0.0, 1.0), **kwargs})


def test_bspline_rejects_points_outside_the_domain() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    with pytest.raises(ValueError, match="must lie in"):
        basis(np.array([0.5, 1.5]))


@given(
    order=st.integers(min_value=1, max_value=6),
    extra=st.integers(min_value=0, max_value=8),
)
@settings(max_examples=25, deadline=None)
def test_bspline_partition_of_unity_property(order: int, extra: int) -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=order + extra, order=order)
    mat = basis(np.linspace(0.0, 1.0, 37))
    np.testing.assert_allclose(mat.sum(axis=1), np.ones(37), rtol=0.0, atol=1e-13)


# --------------------------------------------------------------------------- #
# Fourier
# --------------------------------------------------------------------------- #


def test_fourier_rounds_even_sizes_up() -> None:
    assert Fourier(domain=(0.0, 1.0), n_basis=4).n_basis == 5
    assert Fourier(domain=(0.0, 1.0), n_basis=5).n_basis == 5


def test_fourier_period_defaults_to_the_domain_width() -> None:
    assert Fourier(domain=(0.0, 365.0), n_basis=3).period == 365.0
    assert Fourier(domain=(0.0, 1.0), n_basis=3, period=0.5).period == 0.5


def test_fourier_is_orthonormal_over_one_period() -> None:
    basis = Fourier(domain=(0.0, 365.0), n_basis=9)
    np.testing.assert_allclose(basis.gram(), np.eye(9), rtol=0.0, atol=1e-14)


def test_fourier_values_match_the_definition() -> None:
    basis = Fourier(domain=(0.0, 1.0), n_basis=5)
    t = np.linspace(0.0, 1.0, 13)
    omega = 2.0 * pi
    want = np.column_stack(
        [
            np.full_like(t, 1.0),
            sqrt(2.0) * np.sin(omega * t),
            sqrt(2.0) * np.cos(omega * t),
            sqrt(2.0) * np.sin(2.0 * omega * t),
            sqrt(2.0) * np.cos(2.0 * omega * t),
        ]
    )
    np.testing.assert_allclose(basis(t), want, rtol=1e-13, atol=1e-14)


def test_fourier_derivative_matches_finite_differences() -> None:
    basis = Fourier(domain=(0.0, 2.0), n_basis=7)
    t = np.linspace(0.1, 1.9, 15)
    h = 1e-6
    for deriv in (1, 2):
        want = (basis(t + h, deriv - 1) - basis(t - h, deriv - 1)) / (2.0 * h)
        np.testing.assert_allclose(basis(t, deriv), want, rtol=1e-6, atol=1e-6)


def test_fourier_constant_has_zero_derivative() -> None:
    basis = Fourier(domain=(0.0, 1.0), n_basis=5)
    np.testing.assert_allclose(basis(np.linspace(0.0, 1.0, 9), deriv=1)[:, 0], np.zeros(9))


def test_fourier_penalty_is_diagonal_over_whole_periods() -> None:
    basis = Fourier(domain=(0.0, 1.0), n_basis=7)
    pen = basis.penalty(2)
    np.testing.assert_allclose(pen - np.diag(np.diag(pen)), np.zeros((7, 7)))
    omega = 2.0 * pi
    assert float(pen[1, 1]) == pytest.approx(omega**4)


def test_fourier_penalty_handles_a_partial_period() -> None:
    basis = Fourier(domain=(0.0, 1.0), n_basis=5, period=2.0)
    nodes, weights = la.gauss_legendre(200, 0.0, 1.0)
    mat = basis(nodes, deriv=1)
    want = mat.T @ (weights[:, None] * mat)
    np.testing.assert_allclose(basis.penalty(1), want, rtol=1e-9, atol=1e-10)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [({"n_basis": 0}, "at least 1"), ({"n_basis": 3, "period": 0.0}, "must be positive")],
)
def test_fourier_rejects_bad_arguments(kwargs: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Fourier(domain=(0.0, 1.0), **kwargs)


# --------------------------------------------------------------------------- #
# Monomial / Power / Exponential
# --------------------------------------------------------------------------- #


def test_monomial_values_and_derivatives() -> None:
    basis = Monomial(domain=(0.0, 2.0), n_basis=4)
    np.testing.assert_allclose(basis(np.array([2.0])), [[1.0, 2.0, 4.0, 8.0]])
    np.testing.assert_allclose(basis(np.array([2.0]), 1), [[0.0, 1.0, 4.0, 12.0]])
    np.testing.assert_allclose(basis(np.array([2.0]), 2), [[0.0, 0.0, 2.0, 12.0]])
    np.testing.assert_allclose(basis(np.array([2.0]), 4), [[0.0, 0.0, 0.0, 0.0]])


def test_monomial_gram_is_the_hilbert_matrix() -> None:
    basis = Monomial(domain=(0.0, 1.0), n_basis=4)
    index = np.arange(4)
    want = 1.0 / (index[:, None] + index[None, :] + 1.0)
    np.testing.assert_allclose(basis.gram(), want, rtol=1e-14, atol=1e-15)


def test_monomial_custom_exponents() -> None:
    basis = Monomial(domain=(0.0, 1.0), exponents=[0, 2, 5])
    assert basis.exponents == (0, 2, 5)
    assert basis.n_basis == 3
    np.testing.assert_allclose(basis(np.array([2.0])), [[1.0, 4.0, 32.0]])


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"n_basis": 0}, "at least 1"),
        ({"exponents": [0, -1]}, "non-negative"),
        ({"exponents": [1, 1]}, "distinct"),
    ],
)
def test_monomial_rejects_bad_arguments(kwargs: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Monomial(domain=(0.0, 1.0), **kwargs)


def test_power_values_and_derivatives() -> None:
    basis = Power(domain=(0.25, 4.0), exponents=[0.0, 0.5, -1.0])
    np.testing.assert_allclose(basis(np.array([4.0])), [[1.0, 2.0, 0.25]])
    np.testing.assert_allclose(basis(np.array([4.0]), 1), [[0.0, 0.25, -1.0 / 16.0]])


def test_power_penalty_uses_the_logarithmic_branch() -> None:
    basis = Power(domain=(0.5, 3.0), exponents=[-1.0, 0.0, 1.0])
    pen = basis.penalty(0)
    assert float(pen[0, 1]) == pytest.approx(np.log(3.0 / 0.5))


def test_power_penalty_rejects_a_divergent_integral() -> None:
    basis = Power(domain=(0.0, 1.0), exponents=[0.0, 0.25])
    with pytest.raises(ValueError, match="diverges"):
        basis.penalty(1)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"domain": (-1.0, 1.0), "exponents": [0.0]}, "non-negative domain"),
        ({"domain": (0.0, 1.0), "exponents": []}, "at least one exponent"),
        ({"domain": (0.0, 1.0), "exponents": [1.0, 1.0]}, "distinct"),
        ({"domain": (0.0, 1.0), "exponents": [-1.0]}, "domain excluding zero"),
    ],
)
def test_power_rejects_bad_arguments(kwargs: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Power(**kwargs)


def test_exponential_values_and_derivatives() -> None:
    basis = Exponential(domain=(0.0, 1.0), rates=[0.0, 2.0])
    t = np.array([0.5])
    np.testing.assert_allclose(basis(t), [[1.0, np.exp(1.0)]])
    np.testing.assert_allclose(basis(t, 1), [[0.0, 2.0 * np.exp(1.0)]])
    np.testing.assert_allclose(basis(t, 3), [[0.0, 8.0 * np.exp(1.0)]])


def test_exponential_penalty_handles_cancelling_rates() -> None:
    basis = Exponential(domain=(0.0, 2.0), rates=[1.0, -1.0])
    # exp(t) * exp(-t) == 1, so the off-diagonal entry is the interval length.
    assert float(basis.gram()[0, 1]) == pytest.approx(2.0)


@pytest.mark.parametrize(("rates", "match"), [([], "at least one rate"), ([1.0, 1.0], "distinct")])
def test_exponential_rejects_bad_arguments(rates: list[float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Exponential(domain=(0.0, 1.0), rates=rates)


# --------------------------------------------------------------------------- #
# Constant / Polygonal
# --------------------------------------------------------------------------- #


def test_constant_is_one_with_zero_derivatives() -> None:
    basis = Constant(domain=(0.0, 2.0))
    assert basis.n_basis == 1
    np.testing.assert_allclose(basis(np.linspace(0.0, 2.0, 5)), np.ones((5, 1)))
    np.testing.assert_allclose(basis(np.linspace(0.0, 2.0, 5), 1), np.zeros((5, 1)))
    np.testing.assert_allclose(basis.gram(), [[2.0]])
    np.testing.assert_allclose(basis.penalty(1), [[0.0]])


def test_polygonal_interpolates_linearly_between_vertices() -> None:
    basis = Polygonal([0.0, 1.0, 2.0])
    np.testing.assert_allclose(
        basis(np.array([0.0, 0.5, 1.0, 1.25, 2.0])),
        [
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.75, 0.25],
            [0.0, 0.0, 1.0],
        ],
    )


def test_polygonal_is_a_partition_of_unity() -> None:
    basis = Polygonal([0.0, 0.3, 0.55, 1.0])
    mat = basis(np.linspace(0.0, 1.0, 41))
    np.testing.assert_allclose(mat.sum(axis=1), np.ones(41), rtol=0.0, atol=1e-15)


def test_polygonal_derivative_is_the_local_slope() -> None:
    basis = Polygonal([0.0, 0.5, 1.0])
    np.testing.assert_allclose(basis(np.array([0.25]), 1), [[-2.0, 2.0, 0.0]])


def test_polygonal_rejects_points_outside_the_domain() -> None:
    with pytest.raises(ValueError, match="must lie in"):
        Polygonal([0.0, 1.0])(np.array([-0.5]))


@pytest.mark.parametrize(
    ("argvals", "match"),
    [([0.0], "at least two arguments"), ([0.0, 0.0, 1.0], "strictly increasing")],
)
def test_polygonal_rejects_bad_arguments(argvals: list[float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Polygonal(argvals)


# --------------------------------------------------------------------------- #
# products
# --------------------------------------------------------------------------- #


def test_spline_product_merges_breaks_and_adds_orders() -> None:
    left = BSpline(domain=(0.0, 1.0), breaks=[0.0, 0.5, 1.0], order=4)
    right = BSpline(domain=(0.0, 1.0), breaks=[0.0, 0.25, 1.0], order=3)
    product = left * right
    assert isinstance(product, BSpline)
    assert product.order == 6
    # multiplicities restore the exact smoothness of the product: it is C^1 at
    # 0.25 (from the order-3 factor) and C^2 at 0.5 (from the order-4 factor).
    assert product.breaks == (0.0, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 1.0)
    assert product.n_basis == left.n_basis * right.n_basis - 7


def test_fourier_product_covers_every_harmonic() -> None:
    product = Fourier(domain=(0.0, 1.0), n_basis=5) * Fourier(domain=(0.0, 1.0), n_basis=7)
    assert isinstance(product, Fourier)
    assert product.n_basis == 13


def test_constant_is_the_product_identity() -> None:
    basis = Monomial(domain=(0.0, 1.0), n_basis=3)
    assert Constant(domain=(0.0, 1.0)) * basis == basis
    assert basis * Constant(domain=(0.0, 1.0)) == basis


def test_monomial_and_power_products_sum_exponents() -> None:
    mono = Monomial(domain=(0.0, 1.0), exponents=[0, 2]) * Monomial(
        domain=(0.0, 1.0), exponents=[1, 3]
    )
    assert isinstance(mono, Monomial)
    assert mono.exponents == (1, 3, 5)
    power = Power(domain=(1.0, 2.0), exponents=[0.5]) * Power(
        domain=(1.0, 2.0), exponents=[0.5, 1.0]
    )
    assert isinstance(power, Power)
    assert power.exponents == (1.0, 1.5)


def test_exponential_product_sums_rates() -> None:
    product = Exponential(domain=(0.0, 1.0), rates=[0.0, 1.0]) * Exponential(
        domain=(0.0, 1.0), rates=[2.0]
    )
    assert isinstance(product, Exponential)
    assert product.rates == (2.0, 3.0)


def test_evaluation_rejects_a_two_dimensional_argument() -> None:
    """A grid of points used to be silently flattened, hiding the caller's mistake."""
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    with pytest.raises(ValueError, match="one-dimensional"):
        basis(np.linspace(0.0, 1.0, 6).reshape(3, 2))


def test_evaluation_accepts_a_scalar_argument() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    np.testing.assert_allclose(basis(np.float64(0.5)), basis(np.array([0.5])), atol=1e-14)


def test_mixed_product_falls_back_to_a_rich_spline() -> None:
    left = Fourier(domain=(0.0, 1.0), n_basis=5)
    right = Monomial(domain=(0.0, 1.0), n_basis=3)
    product = left * right
    assert isinstance(product, BSpline)
    assert product.order >= 8
    assert product.n_basis >= 8


def test_fallback_product_reproduces_every_pairwise_product() -> None:
    """The fallback basis must hold each phi_i psi_j, not merely resemble it."""
    pairs = [
        (Fourier(domain=(0.0, 1.0), n_basis=5), Monomial(domain=(0.0, 1.0), n_basis=3)),
        (BSpline(domain=(0.0, 1.0), n_basis=10), Fourier(domain=(0.0, 1.0), n_basis=9)),
        (
            Fourier(domain=(0.0, 1.0), n_basis=7, period=1.0),
            Fourier(domain=(0.0, 1.0), n_basis=7, period=2.0),
        ),
    ]
    t = np.linspace(0.0, 1.0, 501)
    for left, right in pairs:
        product = left * right
        columns = (left(t)[:, :, None] * right(t)[:, None, :]).reshape(t.size, -1)
        residual = columns - product(t) @ np.linalg.lstsq(product(t), columns, rcond=None)[0]
        relative = np.max(np.abs(residual), axis=0) / np.max(np.abs(columns), axis=0)
        assert np.max(relative) <= 1e-8, (type(left).__name__, type(right).__name__)


def test_fallback_product_raises_when_the_product_is_not_a_spline() -> None:
    """sqrt(t) has an unbounded derivative at 0; no polynomial mesh resolves it."""
    left = Power(domain=(0.0, 1.0), exponents=[0.5])
    right = Fourier(domain=(0.0, 1.0), n_basis=3)
    with pytest.raises(ValueError, match="holds the product of Power and Fourier"):
        _ = left * right


def test_penalty_of_a_constant_non_derivative_operator_is_cached() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    la.clear_gram_cache()
    first = basis.penalty(LDO(weights=[1.0, 0.0]))
    assert len(la._GRAM_CACHE) == 1
    np.testing.assert_allclose(basis.penalty(LDO(weights=[1.0, 0.0])), first, atol=1e-14)
    assert len(la._GRAM_CACHE) == 1


def test_fallback_product_carries_the_spline_factor_knot_multiplicity() -> None:
    """A C^2 spline factor forces repeated knots, or the kink cannot be held."""
    spline = BSpline(domain=(0.0, 1.0), order=4, breaks=[0.0, 0.5, 1.0])
    product = spline * Fourier(domain=(0.0, 1.0), n_basis=5)
    assert isinstance(product, BSpline)
    assert product.breaks.count(0.5) == product.order - 3


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (BSpline(domain=(0.0, 1.0), n_basis=5), BSpline(domain=(0.0, 2.0), n_basis=5)),
        (Fourier(domain=(0.0, 1.0), n_basis=3), Monomial(domain=(0.0, 2.0), n_basis=3)),
    ],
)
def test_product_requires_a_shared_domain(left: Basis, right: Basis) -> None:
    with pytest.raises(ValueError, match="domains differ"):
        left * right


def test_product_of_a_basis_and_a_non_basis_is_not_implemented() -> None:
    with pytest.raises(TypeError):
        Constant(domain=(0.0, 1.0)) * 2.0  # type: ignore[operator]


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (BSpline(domain=(0.0, 1.0), n_basis=6), BSpline(domain=(0.0, 1.0), n_basis=5, order=3)),
        (Fourier(domain=(0.0, 1.0), n_basis=5), Fourier(domain=(0.0, 1.0), n_basis=3)),
        (Polygonal([0.0, 0.4, 1.0]), Polygonal([0.0, 0.7, 1.0])),
    ],
)
def test_product_basis_can_represent_the_pointwise_product(left: Basis, right: Basis) -> None:
    rng = np.random.default_rng(11)
    product = left * right
    t = np.linspace(0.0, 1.0, 400)
    target = (left(t) @ rng.normal(size=left.n_basis)) * (right(t) @ rng.normal(size=right.n_basis))
    fitted = product(t) @ la.lstsq(product(t), target[:, None])[:, 0]
    np.testing.assert_allclose(fitted, target, rtol=1e-8, atol=1e-9)
