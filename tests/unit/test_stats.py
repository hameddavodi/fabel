"""Unit and property tests for :mod:`fabel.stats`."""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fabel import LDO, BSpline, FData, Fourier, inprod
from fabel.regression import fregress
from fabel.stats import (
    BoxplotResult,
    DepthResult,
    PermutationTestResult,
    boxplot,
    cor,
    cov,
    depth,
    f_test,
    t_test,
)

# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #


def _random_fd(seed: int, n_curves: int = 9, n_basis: int = 7, shift: float = 0.0) -> FData:
    rng = np.random.default_rng(seed)
    basis = BSpline(domain=(0.0, 1.0), n_basis=n_basis)
    return FData(rng.normal(size=(n_basis, n_curves)) + shift, basis)


def _unique_values(n_points: int, n_curves: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n_points, n_curves)) + np.linspace(0.0, 1.0, n_points)[:, None]


# --------------------------------------------------------------------------- #
# brute-force references
# --------------------------------------------------------------------------- #


def _mbd_brute(values: np.ndarray) -> np.ndarray:
    """Modified band depth by counting every band (Lopez-Pintado and Romo)."""
    n = values.shape[1]
    out = np.zeros(n)
    pairs = list(itertools.combinations(range(n), 2))
    for k in range(n):
        inside = 0.0
        for i, j in pairs:
            low = np.minimum(values[:, i], values[:, j])
            high = np.maximum(values[:, i], values[:, j])
            inside += np.mean((low <= values[:, k]) & (values[:, k] <= high))
        out[k] = inside / len(pairs)
    return out


def _band_brute(values: np.ndarray) -> np.ndarray:
    """Band depth J = 2: the share of bands containing the whole curve."""
    n = values.shape[1]
    pairs = list(itertools.combinations(range(n), 2))
    out = np.zeros(n)
    for k in range(n):
        for i, j in pairs:
            low = np.minimum(values[:, i], values[:, j])
            high = np.maximum(values[:, i], values[:, j])
            out[k] += np.all((low <= values[:, k]) & (values[:, k] <= high))
    return out / len(pairs)


def _fm_brute(values: np.ndarray) -> np.ndarray:
    n = values.shape[1]
    ecdf = (values[:, None, :] <= values[:, :, None]).sum(axis=2) / n
    return np.mean(1.0 - np.abs(0.5 - ecdf), axis=0)


# --------------------------------------------------------------------------- #
# cov / cor
# --------------------------------------------------------------------------- #


def test_cov_without_second_argument_is_the_fdata_covariance() -> None:
    fd = _random_fd(0)
    np.testing.assert_allclose(cov(fd).coefs, fd.cov().coefs, rtol=1e-14)


def test_cross_covariance_matches_pointwise_sample_covariance() -> None:
    x = _random_fd(1)
    y = FData(np.random.default_rng(2).normal(size=(5, 9)), Fourier(domain=(0.0, 1.0), n_basis=5))
    s = np.linspace(0.0, 1.0, 7)
    t = np.linspace(0.0, 1.0, 5)
    surface = cov(x, y)
    assert surface.sbasis == x.basis
    assert surface.tbasis == y.basis
    xs, yt = x(s), y(t)
    expected = (xs - xs.mean(1, keepdims=True)) @ (yt - yt.mean(1, keepdims=True)).T / 8
    np.testing.assert_allclose(surface(s, t), expected, rtol=1e-12, atol=1e-14)


def test_cross_covariance_is_the_transpose_of_the_reverse() -> None:
    x, y = _random_fd(3), _random_fd(4)
    grid = np.linspace(0.0, 1.0, 6)
    np.testing.assert_allclose(cov(x, y)(grid, grid), cov(y, x)(grid, grid).T, rtol=1e-13)


def test_cov_rejects_mismatched_inputs() -> None:
    with pytest.raises(ValueError, match="same number of curves"):
        cov(_random_fd(0, n_curves=4), _random_fd(1, n_curves=5))
    with pytest.raises(ValueError, match="at least two curves"):
        cov(_random_fd(0, n_curves=1), _random_fd(1, n_curves=1))
    multi = FData(np.ones((7, 3, 2)), BSpline(domain=(0.0, 1.0), n_basis=7))
    with pytest.raises(ValueError, match="single variable"):
        cov(multi, multi)
    with pytest.raises(ValueError, match="domain"):
        cov(_random_fd(0), FData(np.ones((4, 9)), BSpline(domain=(0.0, 2.0), n_basis=4)))


def test_cor_of_a_function_with_itself_has_unit_diagonal() -> None:
    fd = _random_fd(5)
    grid = np.linspace(0.0, 1.0, 11)
    matrix = cor(fd, s=grid, t=grid)
    np.testing.assert_allclose(np.diag(matrix), 1.0, rtol=1e-12)
    np.testing.assert_allclose(matrix, matrix.T, rtol=1e-12)


def test_cor_default_grid_has_101_points() -> None:
    assert cor(_random_fd(6)).shape == (101, 101)
    assert cor(_random_fd(6), s=np.linspace(0.0, 1.0, 3)).shape == (3, 101)


def test_cor_matches_numpy_corrcoef() -> None:
    x, y = _random_fd(7), _random_fd(8)
    s, t = np.array([0.1, 0.4]), np.array([0.3, 0.9, 1.0])
    stacked = np.corrcoef(np.vstack([x(s), y(t)]))
    np.testing.assert_allclose(cor(x, y, s=s, t=t), stacked[:2, 2:], rtol=1e-12)


def test_cor_multivariate_stacks_variables() -> None:
    rng = np.random.default_rng(9)
    basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    fd = FData(rng.normal(size=(6, 8, 2)), basis)
    grid = np.linspace(0.0, 1.0, 4)
    matrix = cor(fd, s=grid, t=grid)
    assert matrix.shape == (4, 4, 2)
    for v in range(2):
        single = FData(fd.coefs[:, :, v], basis)
        np.testing.assert_allclose(matrix[:, :, v], cor(single, s=grid, t=grid), rtol=1e-12)
    with pytest.raises(ValueError, match="single variable"):
        cor(fd, fd)


@settings(max_examples=25, deadline=None)
@given(seed=st.integers(0, 10_000))
def test_cor_is_bounded(seed: int) -> None:
    x, y = _random_fd(seed), _random_fd(seed + 1)
    matrix = cor(x, y, s=np.linspace(0.0, 1.0, 9), t=np.linspace(0.0, 1.0, 9))
    assert np.all(np.abs(matrix) <= 1.0 + 1e-12)


# --------------------------------------------------------------------------- #
# depth
# --------------------------------------------------------------------------- #


@settings(max_examples=25, deadline=None)
@given(n_points=st.integers(1, 6), n_curves=st.integers(2, 8), seed=st.integers(0, 10_000))
def test_mbd_counts_every_band(n_points: int, n_curves: int, seed: int) -> None:
    values = _unique_values(n_points, n_curves, seed)
    np.testing.assert_allclose(depth(values).depth, _mbd_brute(values), rtol=1e-12)


def test_mbd_with_ties_uses_average_ranks() -> None:
    # R fbplot, method "MBD", on this tied matrix gives these depths.
    values = np.array([[1, 2, 3, 1, 3, 2], [1, 2, 3, 2, 1, 3]], dtype=float)
    ranks = np.array([[1.5, 3.5, 5.5, 1.5, 5.5, 3.5], [1.5, 3.5, 5.5, 3.5, 1.5, 5.5]])
    expected = ((ranks - 1) * (6 - ranks) + 5).mean(axis=0) / 15
    np.testing.assert_allclose(depth(values).depth, expected, rtol=1e-14)


def test_bd2_is_exact_for_curves_that_do_not_cross() -> None:
    offsets = np.array([0.0, 3.0, 1.0, 4.0, 2.0])
    values = np.sin(np.linspace(0.0, 3.0, 7))[:, None] + offsets
    np.testing.assert_allclose(depth(values, method="BD2").depth, _band_brute(values), rtol=1e-14)


def test_bd2_rank_formula() -> None:
    values = _unique_values(5, 7, 3)
    ranks = values.argsort(axis=1).argsort(axis=1) + 1.0
    expected = ((7 - ranks.max(0)) * (ranks.min(0) - 1) + 6) / 21
    np.testing.assert_allclose(depth(values, method="bd2").depth, expected, rtol=1e-14)


def test_both_breaks_band_depth_ties_with_mbd() -> None:
    values = _unique_values(6, 9, 11)
    bd2 = depth(values, method="BD2").depth
    mbd = depth(values, method="MBD").depth
    np.testing.assert_allclose(
        depth(values, method="both").depth, np.round(bd2 * 1e4) + mbd, rtol=1e-14
    )


@settings(max_examples=25, deadline=None)
@given(n_points=st.integers(1, 6), n_curves=st.integers(1, 8), seed=st.integers(0, 10_000))
def test_fm_depth_is_the_mean_ecdf_centrality(n_points: int, n_curves: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    values = rng.integers(0, 4, size=(n_points, n_curves)).astype(float)
    np.testing.assert_allclose(depth(values, method="FM").depth, _fm_brute(values), rtol=1e-14)


@settings(max_examples=20, deadline=None)
@given(seed=st.integers(0, 10_000), method=st.sampled_from(["MBD", "BD2", "FM", "both"]))
def test_depth_is_invariant_under_increasing_transformations(seed: int, method: str) -> None:
    values = _unique_values(5, 7, seed)
    np.testing.assert_allclose(
        depth(np.exp(values), method=method).depth, depth(values, method=method).depth
    )


def test_depth_trimming_and_median() -> None:
    values = _unique_values(8, 13, 21)
    result = depth(values, method="FM", trim=0.25)
    assert isinstance(result, DepthResult)
    assert result.median_index == int(np.argmax(result.depth))
    np.testing.assert_array_equal(result.median, values[:, result.median_index])
    threshold = np.quantile(result.depth, 0.25)
    np.testing.assert_array_equal(result.trimmed_index, np.flatnonzero(result.depth >= threshold))
    np.testing.assert_allclose(result.trimmed_mean, values[:, result.trimmed_index].mean(axis=1))
    assert result.t is None


def test_depth_trim_keeps_exact_ties() -> None:
    # Curves 0 and 4 tie exactly on the 25% quantile of the FM depths, reached
    # through different ranks at each point; both must be kept.
    values = np.array([[0.0, 1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0, 0.0]])
    result = depth(values, method="FM", trim=0.25)
    np.testing.assert_array_equal(result.trimmed_index, [0, 1, 2, 3, 4])
    assert depth(values, method="FM", trim=0.0).trimmed_index.shape == (5,)


def test_depth_of_fdata_uses_a_101_point_grid() -> None:
    fd = _random_fd(12)
    result = depth(fd)
    grid = np.linspace(0.0, 1.0, 101)
    np.testing.assert_allclose(result.t, grid)
    np.testing.assert_allclose(result.depth, depth(fd(grid)).depth)
    custom = depth(fd, t=np.array([0.2, 0.5]))
    assert custom.median.shape == (2,)


def test_depth_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="method"):
        depth(np.ones((3, 4)), method="mode")
    with pytest.raises(ValueError, match="trim"):
        depth(np.ones((3, 4)), trim=1.0)
    with pytest.raises(ValueError, match="two-dimensional"):
        depth(np.ones(4))
    with pytest.raises(ValueError, match="at least two curves"):
        depth(np.ones((4, 1)), method="MBD")
    with pytest.raises(ValueError, match="evaluation points"):
        depth(np.ones((4, 3)), t=np.arange(3.0))
    with pytest.raises(ValueError, match="single variable"):
        depth(FData(np.ones((4, 3, 2)), BSpline(n_basis=4)))
    with pytest.raises(ValueError, match="finite"):
        depth(np.array([[1.0, np.nan], [0.0, 1.0]]))


def test_depth_keeps_explicit_points_for_arrays() -> None:
    values = _unique_values(3, 5, 1)
    points = np.array([0.0, 0.5, 2.0])
    np.testing.assert_array_equal(depth(values, t=points).t, points)


# --------------------------------------------------------------------------- #
# boxplot
# --------------------------------------------------------------------------- #


def _with_outliers() -> np.ndarray:
    rng = np.random.default_rng(31)
    grid = np.linspace(0.0, 1.0, 40)
    values = np.sin(2 * np.pi * grid)[:, None] + 0.2 * rng.normal(size=(40, 20))
    values[:, 3] += 6.0  # magnitude outlier
    values[18:30, 11] = -8.0  # a shape outlier over part of the domain
    return values


def test_boxplot_flags_outliers() -> None:
    result = boxplot(_with_outliers())
    assert isinstance(result, BoxplotResult)
    np.testing.assert_array_equal(result.outliers, [3, 11])
    assert result.median_index == int(np.argmax(result.depth))
    np.testing.assert_array_equal(result.median, _with_outliers()[:, result.median_index])


def test_boxplot_envelopes() -> None:
    values = _with_outliers()
    result = boxplot(values, prob=0.5, factor=1.5)
    order = np.argsort(-result.depth, kind="stable")
    central = values[:, order[:10]]
    np.testing.assert_array_equal(result.central_lower, central.min(axis=1))
    np.testing.assert_array_equal(result.central_upper, central.max(axis=1))
    spread = result.central_upper - result.central_lower
    np.testing.assert_allclose(result.fence_lower, result.central_lower - 1.5 * spread)
    np.testing.assert_allclose(result.fence_upper, result.central_upper + 1.5 * spread)
    kept = np.setdiff1d(np.arange(20), result.outliers)
    np.testing.assert_array_equal(result.whisker_lower, values[:, kept].min(axis=1))
    np.testing.assert_array_equal(result.whisker_upper, values[:, kept].max(axis=1))
    np.testing.assert_array_equal(result.values, values)


def test_boxplot_central_region_rounds_up() -> None:
    values = _unique_values(6, 7, 5)
    result = boxplot(values, prob=0.5)
    order = np.argsort(-result.depth, kind="stable")
    np.testing.assert_array_equal(result.central_upper, values[:, order[:4]].max(axis=1))


def test_boxplot_of_fdata_and_errors() -> None:
    fd = _random_fd(13, n_curves=12)
    result = boxplot(fd, method="both")
    assert result.t is not None
    assert result.values.shape == (101, 12)
    with pytest.raises(ValueError, match="prob"):
        boxplot(fd, prob=0.0)
    with pytest.raises(ValueError, match="factor"):
        boxplot(fd, factor=-1.0)


def test_boxplot_plot() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values = _with_outliers()
    ax = boxplot(values).plot()
    # median, two whiskers, two outliers
    assert len(ax.lines) == 5
    assert len(ax.collections) == 1
    np.testing.assert_allclose(ax.lines[0].get_xydata()[:, 1], boxplot(values).median)
    _, other = plt.subplots()
    grid = np.linspace(0.0, 2.0, 40)
    assert boxplot(values, t=grid).plot(ax=other, color="k") is other
    np.testing.assert_allclose(other.lines[0].get_xydata()[:, 0], grid)
    plt.close("all")


# --------------------------------------------------------------------------- #
# t_test
# --------------------------------------------------------------------------- #


def _two_groups(effect: float) -> tuple[FData, FData]:
    first = _random_fd(40, n_curves=8)
    second = _random_fd(41, n_curves=6, shift=effect)
    return first, second


def test_t_test_statistic_is_the_maximal_welch_t() -> None:
    first, second = _two_groups(0.5)
    result = t_test(first, second, n_perm=20, random_state=0)
    assert isinstance(result, PermutationTestResult)
    grid = np.linspace(0.0, 1.0, 101)
    a, b = first(grid), second(grid)
    welch = np.abs(a.mean(1) - b.mean(1)) / np.sqrt(a.var(1, ddof=1) / 8 + b.var(1, ddof=1) / 6)
    np.testing.assert_allclose(result.pointwise, welch, rtol=1e-12)
    assert result.statistic == pytest.approx(welch.max(), rel=1e-12)
    np.testing.assert_allclose(result.t, grid)


def test_t_test_null_follows_the_permutations() -> None:
    first, second = _two_groups(0.0)
    result = t_test(first, second, n_perm=5, random_state=np.random.default_rng(3))
    rng = np.random.default_rng(3)
    grid = np.linspace(0.0, 1.0, 101)
    pooled = np.hstack([first(grid), second(grid)])
    for k in range(5):
        order = rng.permutation(14)
        a, b = pooled[:, order[:8]], pooled[:, order[8:]]
        stat = np.abs(a.mean(1) - b.mean(1)) / np.sqrt(a.var(1, ddof=1) / 8 + b.var(1, ddof=1) / 6)
        np.testing.assert_allclose(result.pointwise_null[k], stat, rtol=1e-12)
        assert result.null[k] == pytest.approx(stat.max(), rel=1e-12)


def test_t_test_summaries() -> None:
    first, second = _two_groups(3.0)
    result = t_test(first, second, n_perm=40, q=0.1, random_state=7)
    assert result.pvalue == np.mean(result.null >= result.statistic)
    assert result.critical_value == pytest.approx(np.quantile(result.null, 0.9), rel=1e-14)
    np.testing.assert_allclose(
        result.pointwise_pvalue, np.mean(result.pointwise_null >= result.pointwise, axis=0)
    )
    np.testing.assert_allclose(
        result.pointwise_critical_value, np.quantile(result.pointwise_null, 0.9, axis=0)
    )
    assert result.pvalue < 0.05


def test_t_test_is_reproducible_by_seed() -> None:
    first, second = _two_groups(0.2)
    one = t_test(first, second, n_perm=10, random_state=5)
    two = t_test(first, second, n_perm=10, random_state=5)
    np.testing.assert_array_equal(one.null, two.null)
    other = t_test(first, second, n_perm=10, random_state=6)
    assert not np.array_equal(one.null, other.null)
    assert t_test(first, second, n_perm=3).null.shape == (3,)


def test_t_test_custom_points_and_errors() -> None:
    first, second = _two_groups(0.0)
    points = np.array([0.1, 0.2, 0.9])
    assert t_test(first, second, n_perm=2, t=points, random_state=0).pointwise.shape == (3,)
    with pytest.raises(ValueError, match="n_perm"):
        t_test(first, second, n_perm=0)
    with pytest.raises(ValueError, match="q"):
        t_test(first, second, q=1.5)
    with pytest.raises(ValueError, match="two curves"):
        t_test(first[0], second)
    with pytest.raises(ValueError, match="domain"):
        t_test(first, FData(np.ones((4, 3)), BSpline(domain=(0.0, 2.0), n_basis=4)))
    multi = FData(np.ones((7, 3, 2)), BSpline(domain=(0.0, 1.0), n_basis=7))
    with pytest.raises(ValueError, match="single variable"):
        t_test(multi, multi)
    with pytest.raises(TypeError, match="random_state"):
        t_test(first, second, random_state="seed")  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# f_test
# --------------------------------------------------------------------------- #


def _regression_data(seed: int = 50) -> tuple[FData, np.ndarray]:
    rng = np.random.default_rng(seed)
    basis = Fourier(domain=(0.0, 1.0), n_basis=7)
    group = np.repeat([0.0, 1.0], 6)
    coefs = rng.normal(size=(7, 12)) + np.outer(rng.normal(size=7), group)
    return FData(coefs, basis), group


def test_f_test_unpenalised_scalar_covariates_is_pointwise_least_squares() -> None:
    y, group = _regression_data()
    result = f_test(y, [np.ones(12), group], n_perm=4, random_state=0)
    grid = np.linspace(0.0, 1.0, 101)
    values = y(grid)
    design = np.column_stack([np.ones(12), group])
    fitted = design @ np.linalg.lstsq(design, values.T, rcond=None)[0]
    stat = fitted.var(axis=0, ddof=1) / np.mean((values.T - fitted) ** 2, axis=0)
    np.testing.assert_allclose(result.pointwise, stat, rtol=1e-9)
    assert result.statistic == pytest.approx(stat.max(), rel=1e-9)


def test_f_test_null_permutes_the_response() -> None:
    y, group = _regression_data()
    covariates = [np.ones(12), group]
    result = f_test(y, covariates, lam=1e-3, n_perm=3, random_state=np.random.default_rng(8))
    rng = np.random.default_rng(8)
    for k in range(3):
        permuted = y[rng.permutation(12)]
        again = f_test(permuted, covariates, lam=1e-3, n_perm=1, random_state=0)
        np.testing.assert_allclose(result.null[k], again.statistic, rtol=1e-10)
        np.testing.assert_allclose(result.pointwise_null[k], again.pointwise, rtol=1e-10)


def test_f_test_penalty_shrinks_the_fit() -> None:
    y, group = _regression_data()
    loose = f_test(y, [np.ones(12), group], lam=0.0, n_perm=1, random_state=0)
    tight = f_test(y, [np.ones(12), group], lam=[0.0, 1e6], n_perm=1, random_state=0)
    assert tight.statistic < loose.statistic


def test_f_test_functional_covariate_matches_normal_equations() -> None:
    rng = np.random.default_rng(60)
    basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    x = FData(rng.normal(size=(6, 10)), basis)
    signal = x * FData(np.sin(np.linspace(0.0, 1.0, 6)), basis)
    y = signal + FData(0.1 * rng.normal(size=(signal.basis.n_basis, 10)), signal.basis)
    beta_basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    result = f_test(y, x, basis=beta_basis, lam=1e-4, penalty=LDO(2), n_perm=2, random_state=0)
    # Reference: Gauss-Legendre on 60 panels, whose edges include every break
    # (1/3, 1/2, 2/3), so each panel integrates polynomials of degree <= 23 exactly.
    nodes, weights = np.polynomial.legendre.leggauss(12)
    edges = np.linspace(0.0, 1.0, 61)
    grid = (edges[:-1, None] + (nodes[None, :] + 1.0) / 120.0).ravel()
    w = np.tile(weights / 120.0, 60)
    phi = beta_basis(grid)
    xv, yv = x(grid), y(grid)
    lhs = (phi * (w * (xv**2).sum(1))[:, None]).T @ phi + 1e-4 * beta_basis.penalty(2)
    rhs = phi.T @ (w * (xv * yv).sum(1))
    coef = np.linalg.solve(lhs, rhs)
    points = np.linspace(0.0, 1.0, 101)
    fitted = x(points) * (beta_basis(points) @ coef)[:, None]
    stat = fitted.var(axis=1, ddof=1) / np.mean((y(points) - fitted) ** 2, axis=1)
    np.testing.assert_allclose(result.pointwise, stat, rtol=1e-10)


def test_f_test_scalar_response() -> None:
    rng = np.random.default_rng(70)
    group = np.repeat([0.0, 1.0], 10)
    response = rng.normal(size=20) + 2.0 * group
    result = f_test(response, [np.ones(20), group], n_perm=30, q=0.05, random_state=1)
    design = np.column_stack([np.ones(20), group])
    fitted = design @ np.linalg.lstsq(design, response, rcond=None)[0]
    stat = fitted.var(ddof=1) / np.mean((response - fitted) ** 2)
    assert result.statistic == pytest.approx(stat, rel=1e-12)
    assert result.t is None
    assert result.pointwise.shape == (1,)
    assert result.pointwise_null.shape == (30, 1)
    assert result.pvalue < 0.05


def test_f_test_scalar_response_functional_covariate() -> None:
    rng = np.random.default_rng(71)
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    x = FData(rng.normal(size=(5, 15)), basis)
    response = np.asarray(inprod(x, FData(np.ones(5), basis)))[:, 0] + 0.1 * rng.normal(size=15)
    result = f_test(response, [np.ones(15), x], lam=[0.0, 0.0], n_perm=2, random_state=0)
    design = np.column_stack([np.ones(15), np.asarray(inprod(x, basis))])
    fitted = design @ np.linalg.lstsq(design, response, rcond=None)[0]
    stat = fitted.var(ddof=1) / np.mean((response - fitted) ** 2)
    assert result.statistic == pytest.approx(stat, rel=1e-8)


def test_f_test_rejects_bad_input() -> None:
    y, group = _regression_data()
    with pytest.raises(ValueError, match="covariate"):
        f_test(y, [])
    with pytest.raises(ValueError, match="12"):
        f_test(y, [np.ones(11)])
    with pytest.raises(ValueError, match="one entry per covariate"):
        f_test(y, [np.ones(12), group], lam=[1.0])
    with pytest.raises(ValueError, match="one entry per covariate"):
        f_test(y, [np.ones(12)], basis=[y.basis, y.basis])
    with pytest.raises(ValueError, match="one entry per covariate"):
        f_test(y, [np.ones(12)], penalty=[2, 2])
    with pytest.raises(ValueError, match="non-negative"):
        f_test(y, [np.ones(12)], lam=-1.0)
    with pytest.raises(ValueError, match="domain"):
        f_test(y, [np.ones(12)], basis=BSpline(domain=(0.0, 2.0), n_basis=5))
    with pytest.raises(ValueError, match="one-dimensional"):
        f_test(y, [np.ones((12, 2))])
    with pytest.raises(ValueError, match="one-dimensional"):
        f_test(np.ones((3, 3)), [np.ones(3)])
    with pytest.raises(ValueError, match="n_perm"):
        f_test(y, [np.ones(12)], n_perm=0)
    multi = FData(np.ones((7, 12, 2)), y.basis)
    with pytest.raises(ValueError, match="single variable"):
        f_test(multi, [np.ones(12)])
    with pytest.raises(ValueError, match="single variable"):
        f_test(y, [multi])
    with pytest.raises(ValueError, match="singular"):
        f_test(y, [np.ones(12), np.ones(12)])


def test_f_test_accepts_a_single_covariate() -> None:
    y, _ = _regression_data()
    one: Any = f_test(y, np.ones(12), n_perm=2, random_state=0)
    listed = f_test(y, [np.ones(12)], n_perm=2, random_state=0)
    np.testing.assert_allclose(one.pointwise, listed.pointwise)
    # An intercept alone explains nothing but the mean: yhat is constant in i.
    np.testing.assert_allclose(one.pointwise, 0.0, atol=1e-20)


# --------------------------------------------------------------------------- #
# f_test on a fitted fregress model
# --------------------------------------------------------------------------- #


def _same_result(left: PermutationTestResult, right: PermutationTestResult) -> None:
    """Assert two permutation test results are identical, bit for bit."""
    assert left.statistic == right.statistic
    assert left.pvalue == right.pvalue
    assert left.critical_value == right.critical_value
    np.testing.assert_array_equal(left.null, right.null)
    np.testing.assert_array_equal(left.pointwise, right.pointwise)
    np.testing.assert_array_equal(left.pointwise_null, right.pointwise_null)
    np.testing.assert_array_equal(left.pointwise_critical_value, right.pointwise_critical_value)
    if left.t is None or right.t is None:
        assert left.t is None
        assert right.t is None
    else:
        np.testing.assert_array_equal(left.t, right.t)


def test_f_test_on_a_functional_response_model_matches_the_raw_form() -> None:
    y, group = _regression_data()
    beta_basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    model = fregress(y, [np.ones(12), group], [beta_basis, beta_basis], lam=1e-3, penalty=2)
    from_model = f_test(model, n_perm=20, q=0.1, random_state=5)
    raw = f_test(
        y,
        [np.ones(12), group],
        basis=beta_basis,
        lam=1e-3,
        penalty=2,
        n_perm=20,
        q=0.1,
        random_state=5,
    )
    _same_result(from_model, raw)


def test_f_test_on_a_functional_covariate_model_matches_the_raw_form() -> None:
    rng = np.random.default_rng(61)
    basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    x = FData(rng.normal(size=(6, 10)), basis)
    y = x * FData(np.sin(np.linspace(0.0, 1.0, 6)), basis)
    beta_basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    model = fregress(y, {"x": x}, {"x": (beta_basis, 1e-4, LDO(2))})
    points = np.linspace(0.0, 1.0, 31)
    from_model = f_test(model, n_perm=4, t=points, random_state=2)
    raw = f_test(
        y, x, basis=beta_basis, lam=1e-4, penalty=LDO(2), n_perm=4, t=points, random_state=2
    )
    _same_result(from_model, raw)


def test_f_test_on_a_scalar_response_model_matches_the_raw_form() -> None:
    rng = np.random.default_rng(72)
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    x = FData(rng.normal(size=(5, 15)), basis)
    response = np.asarray(inprod(x, FData(np.ones(5), basis)))[:, 0] + 0.1 * rng.normal(size=15)
    model = fregress(response, [np.ones(15), x], lam=1e-3)
    from_model = f_test(model, n_perm=25, random_state=np.random.default_rng(9))
    raw = f_test(
        response,
        [np.ones(15), x],
        lam=[0.0, 1e-3],
        n_perm=25,
        random_state=np.random.default_rng(9),
    )
    _same_result(from_model, raw)
    assert from_model.t is None


def test_f_test_on_a_model_honours_the_intercept() -> None:
    y, group = _regression_data()
    model = fregress(y, {"const": 1.0, "group": group})
    from_model = f_test(model, n_perm=10, random_state=3)
    with_intercept = f_test(y, [np.ones(12), group], n_perm=10, random_state=3)
    without = f_test(y, group, n_perm=10, random_state=3)
    _same_result(from_model, with_intercept)
    assert not np.allclose(from_model.pointwise, without.pointwise)


def test_f_test_on_a_formula_model_matches_the_raw_form() -> None:
    y, group = _regression_data()
    labels = ["a" if g == 0.0 else "b" for g in group]
    model = fregress("y ~ g", {"y": y, "g": labels})
    assert model.names == ("const", "g.b")
    _same_result(
        f_test(model, n_perm=5, random_state=4),
        f_test(y, [np.ones(12), group], n_perm=5, random_state=4),
    )


def test_f_test_on_a_weighted_scalar_model_is_weighted_least_squares() -> None:
    rng = np.random.default_rng(73)
    group = np.repeat([0.0, 1.0], 10)
    response = rng.normal(size=20) + 2.0 * group
    weights = rng.uniform(0.5, 2.0, size=20)
    model = fregress(response, [1.0, group], weights=weights)
    result = f_test(model, n_perm=6, random_state=0)
    design = np.column_stack([np.ones(20), group])
    root = np.sqrt(weights)
    perm = np.random.default_rng(0)
    orders = [np.arange(20), *(perm.permutation(20) for _ in range(6))]
    expected = []
    for order in orders:
        coef = np.linalg.lstsq(design * root[:, None], response[order] * root, rcond=None)[0]
        fitted = design @ coef
        expected.append(fitted.var(ddof=1) / np.mean((response[order] - fitted) ** 2))
    assert result.statistic == pytest.approx(expected[0], rel=1e-10)
    np.testing.assert_allclose(result.null, expected[1:], rtol=1e-10)


def test_f_test_on_a_weighted_functional_model_is_weighted_least_squares() -> None:
    y, group = _regression_data()
    weights = np.linspace(0.5, 2.0, 12)
    model = fregress(y, [1.0, group], weights=weights)
    result = f_test(model, n_perm=3, random_state=1)
    # Scalar covariates with beta in the response basis: the concurrent model is
    # pointwise weighted least squares, the same hat matrix at every t.
    values = y(np.linspace(0.0, 1.0, 101))
    design = np.column_stack([np.ones(12), group])
    hat = design @ np.linalg.solve(design.T @ (weights[:, None] * design), design.T * weights)
    perm = np.random.default_rng(1)
    orders = [np.arange(12), *(perm.permutation(12) for _ in range(3))]
    expected = []
    for order in orders:
        fitted = values[:, order] @ hat.T
        expected.append(
            fitted.var(axis=1, ddof=1) / np.mean((values[:, order] - fitted) ** 2, axis=1)
        )
    np.testing.assert_allclose(result.pointwise, expected[0], rtol=1e-9)
    np.testing.assert_allclose(result.pointwise_null, np.array(expected[1:]), rtol=1e-9)
    unweighted = f_test(y, [np.ones(12), group], n_perm=3, random_state=1)
    assert not np.allclose(result.pointwise, unweighted.pointwise)


def test_f_test_on_a_model_rejects_raw_settings() -> None:
    y, group = _regression_data()
    model = fregress(y, [1.0, group])
    with pytest.raises(TypeError, match="x"):
        f_test(model, [np.ones(12), group])
    with pytest.raises(TypeError, match="basis"):
        f_test(model, basis=y.basis)  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="lam"):
        f_test(model, lam=1.0)  # type: ignore[call-overload]
    with pytest.raises(TypeError, match="penalty"):
        f_test(model, penalty=1)  # type: ignore[call-overload]
    with pytest.raises(ValueError, match="n_perm"):
        f_test(model, n_perm=0)
    with pytest.raises(TypeError, match="covariates"):
        f_test(y)  # type: ignore[call-overload]
