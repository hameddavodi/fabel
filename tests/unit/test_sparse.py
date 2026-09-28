"""Unit tests for :mod:`fabel.sparse` (PACE for sparse longitudinal data)."""

from __future__ import annotations

import pickle
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.base import clone
from sklearn.exceptions import NotFittedError

from fabel import LDO, BSpline, FData, Monomial
from fabel.sparse import (
    PACE,
    SparseCov,
    _curves,
    _harmonics,
    _whitened_sum_signs,
    sparse_cov,
    sparse_mean,
)

# Small simulated samples often give a non-positive sigma^2 estimate (the
# surface's sampling error on its diagonal exceeds the noise variance); the
# floor warning is tested on its own below.
pytestmark = pytest.mark.filterwarnings(
    "ignore:the estimated measurement-error variance:RuntimeWarning"
)

UNIT = (0.0, 1.0)


def _simulate(
    n_curves: int = 150,
    n_points: int = 6,
    noise: float = 0.2,
    seed: int = 0,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    """Rank-two model x_i(t) = 1 + t + a_i sqrt(2) sin(pi t) + b_i sqrt(2) cos(pi t)."""
    rng = np.random.default_rng(seed)
    times = [np.sort(rng.uniform(0.0, 1.0, n_points)) for _ in range(n_curves)]
    scores = np.column_stack([rng.normal(0.0, 2.0, n_curves), rng.normal(0.0, 0.7, n_curves)])
    values = [
        1.0
        + t
        + s[0] * np.sqrt(2.0) * np.sin(np.pi * t)
        + s[1] * np.sqrt(2.0) * np.cos(np.pi * t)
        + noise * rng.normal(size=t.shape)
        for t, s in zip(times, scores, strict=True)
    ]
    return times, values, scores


def _true_cov(s: np.ndarray, t: np.ndarray) -> np.ndarray:
    root2 = np.sqrt(2.0)
    first = 4.0 * np.outer(root2 * np.sin(np.pi * s), root2 * np.sin(np.pi * t))
    second = 0.49 * np.outer(root2 * np.cos(np.pi * s), root2 * np.cos(np.pi * t))
    return np.asarray(first + second)


# --------------------------------------------------------------------------- #
# input handling
# --------------------------------------------------------------------------- #


def test_curves_list_form_matches_pairs() -> None:
    t = [np.array([0.0, 0.5]), np.array([0.2, 0.4, 0.9])]
    y = [np.array([1.0, 2.0]), np.array([3.0, 4.0, 5.0])]
    pairs = [np.column_stack([a, b]) for a, b in zip(t, y, strict=True)]
    times, values = _curves(pairs, None)
    for got, want in zip(times + values, t + y, strict=True):
        np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize(
    ("y", "t", "error"),
    [
        ([], None, ValueError),
        ("abc", None, TypeError),
        (5.0, None, TypeError),
        ([np.array([1.0, 2.0])], None, ValueError),
        ([np.ones((2, 3))], None, ValueError),
        ([np.array([1.0])], [np.array([0.0]), np.array([1.0])], ValueError),
        ([np.array([1.0])], "t", ValueError),
        ([np.array([1.0, 2.0])], [np.array([0.0])], ValueError),
        ([np.array([])], [np.array([])], ValueError),
        ([np.array([np.nan])], [np.array([0.0])], ValueError),
        ([np.array([1.0])], [np.array([np.inf])], ValueError),
    ],
)
def test_curves_rejects_bad_input(y: object, t: object, error: type[Exception]) -> None:
    with pytest.raises(error):
        _curves(y, t)


# --------------------------------------------------------------------------- #
# mean
# --------------------------------------------------------------------------- #


def test_sparse_mean_reproduces_a_function_in_the_span() -> None:
    rng = np.random.default_rng(1)
    times = [np.sort(rng.uniform(0.0, 1.0, 3)) for _ in range(20)]
    values = [t**3 - 2.0 * t for t in times]
    mean = sparse_mean(values, times, BSpline(domain=UNIT, n_basis=6))
    grid = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose(mean(grid)[:, 0], grid**3 - 2.0 * grid, atol=1e-10)
    assert mean.coefs.shape == (6, 1)


def test_sparse_mean_penalty_shrinks_to_a_line() -> None:
    times, values, _ = _simulate(n_curves=40)
    basis = BSpline(domain=UNIT, n_basis=8)
    rough = sparse_mean(values, times, basis, lam=0.0)
    flat = sparse_mean(values, times, basis, lam=1e8, penalty=LDO(2))
    grid = np.linspace(0.0, 1.0, 21)
    assert np.max(np.abs(flat(grid, 2))) < 1e-3 < np.max(np.abs(rough(grid, 2)))


def test_sparse_mean_default_basis_covers_the_data() -> None:
    times, values, _ = _simulate(n_curves=10, n_points=3)
    mean = sparse_mean(values, times)
    pooled = np.concatenate(times)
    np.testing.assert_allclose(mean.domain, (pooled.min(), pooled.max()))
    assert mean.basis.n_basis == 10


def test_sparse_mean_default_basis_is_capped_by_distinct_times() -> None:
    t = [np.array([0.0, 1.0]), np.array([0.0, 1.0, 2.0])]
    y = [np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0])]
    assert sparse_mean(y, t).basis.n_basis == 3


@pytest.mark.parametrize("lam", [-1.0, float("inf"), float("nan")])
def test_sparse_mean_rejects_bad_lambda(lam: float) -> None:
    with pytest.raises(ValueError, match="lam"):
        sparse_mean([np.array([1.0, 2.0])], [np.array([0.0, 1.0])], lam=lam)


def test_sparse_mean_rejects_an_undetermined_basis() -> None:
    t = [np.array([0.1, 0.2])]
    with pytest.raises(ValueError, match="singular"):
        sparse_mean([np.array([1.0, 2.0])], t, BSpline(domain=UNIT, n_basis=8))


def test_sparse_mean_rejects_a_single_time() -> None:
    with pytest.raises(ValueError, match="non-degenerate"):
        sparse_mean([np.array([1.0, 2.0])], [np.array([0.5, 0.5])])


# --------------------------------------------------------------------------- #
# covariance
# --------------------------------------------------------------------------- #


def test_sparse_cov_recovers_the_surface_and_the_noise() -> None:
    times, values, _ = _simulate(n_curves=6000, n_points=8, noise=0.5, seed=2)
    est = sparse_cov(values, times, basis=BSpline(domain=UNIT, n_basis=6), lam=1e-4)
    assert isinstance(est, SparseCov)
    grid = np.linspace(0.1, 0.9, 9)
    # Measured: 0.24 (surface), 0.26 (diagonal), sigma^2 = 0.261 against 0.25.
    assert np.max(np.abs(est.cov(grid, grid) - _true_cov(grid, grid))) < 0.4
    assert est.sigma2 == pytest.approx(0.25, abs=0.05)
    np.testing.assert_allclose(est.cov.coefs, est.cov.coefs.T)
    # V(t) = G(t, t) + sigma^2 on the diagonal.
    diag = np.diag(_true_cov(grid, grid)) + 0.25
    assert np.max(np.abs(est.variance(grid)[:, 0] - diag)) < 0.4


def test_sparse_cov_excludes_the_diagonal() -> None:
    """Pure noise with no curve effect: off-diagonal pairs average to zero."""
    rng = np.random.default_rng(3)
    times = [np.sort(rng.uniform(0.0, 1.0, 5)) for _ in range(400)]
    values = [rng.normal(0.0, 1.0, 5) for _ in times]
    basis = BSpline(domain=UNIT, n_basis=1, order=1)
    mean = FData(np.zeros((1, 1)), basis)
    est = sparse_cov(values, times, mean=mean, basis=basis)
    assert abs(float(est.cov.coefs[0, 0])) < 0.05
    assert est.sigma2 == pytest.approx(1.0, abs=0.1)


def test_sparse_cov_matches_brute_force_pairs() -> None:
    rng = np.random.default_rng(4)
    times = [np.sort(rng.uniform(0.0, 1.0, int(rng.integers(1, 6)))) for _ in range(25)]
    values = [rng.normal(size=t.shape) for t in times]
    basis = BSpline(domain=UNIT, n_basis=4)
    mean = FData(np.full((4, 1), 0.1), basis)
    est = sparse_cov(values, times, mean=mean, basis=basis, lam=0.01)
    rows, target = [], []
    for t, y in zip(times, values, strict=True):
        psi = basis(t)
        r = y - 0.1
        for j in range(len(t)):
            for k in range(len(t)):
                if j != k:
                    rows.append(np.kron(psi[j], psi[k]))
                    target.append(r[j] * r[k])
    a = np.asarray(rows)
    w, p = np.asarray(basis.gram()), np.asarray(basis.penalty(2))
    pen = np.kron(w, p) + np.kron(p, w)
    want = np.linalg.solve(a.T @ a + 0.01 * pen, a.T @ np.asarray(target)).reshape(4, 4)
    np.testing.assert_allclose(est.cov.coefs, want, rtol=1e-9, atol=1e-12)


def test_sparse_cov_default_mean_and_basis() -> None:
    times, values, _ = _simulate(n_curves=80)
    est = sparse_cov(values, times)
    assert est.cov.coefs.shape == (6, 6)
    assert est.mean.basis.n_basis == 10


def test_sparse_cov_needs_pairs() -> None:
    t = [np.array([0.0]), np.array([1.0])]
    with pytest.raises(ValueError, match="two or more"):
        sparse_cov([np.array([1.0]), np.array([2.0])], t, basis=BSpline(n_basis=1, order=1))


def test_sparse_cov_rejects_a_multi_curve_mean() -> None:
    times, values, _ = _simulate(n_curves=10)
    basis = BSpline(domain=UNIT, n_basis=4)
    with pytest.raises(ValueError, match="single curve"):
        sparse_cov(values, times, mean=FData(np.zeros((4, 2)), basis), basis=basis)


def test_sparse_cov_warns_and_floors_a_non_positive_sigma2() -> None:
    # Noise-free constant curves: the pairs weight a curve by n(n-1), the
    # diagonal by n, so the long curve with the large level pulls G above V.
    times = [np.linspace(0.0, 1.0, 5), np.array([0.2, 0.8])]
    values = [np.full(5, 3.0), np.full(2, 0.1)]
    basis = BSpline(domain=UNIT, n_basis=1, order=1)
    mean = FData(np.zeros((1, 1)), basis)
    with pytest.warns(RuntimeWarning, match="not positive"):
        est = sparse_cov(values, times, mean=mean, basis=basis)
    assert 0.0 < est.sigma2 < 1e-3


def test_sparse_cov_rejects_an_undetermined_surface() -> None:
    times, values, _ = _simulate(n_curves=5, n_points=2)
    basis = BSpline(domain=UNIT, n_basis=8)
    with pytest.raises(ValueError, match="covariance surface"):
        sparse_cov(values, times, mean=FData(np.zeros((8, 1)), basis), basis=basis)


@settings(max_examples=25, deadline=None)
@given(scale=st.floats(0.1, 10.0), seed=st.integers(0, 10_000))
def test_sparse_cov_scales_quadratically(scale: float, seed: int) -> None:
    times, values, _ = _simulate(n_curves=30, n_points=5, seed=seed)
    basis = BSpline(domain=UNIT, n_basis=3, order=2)
    mean = FData(np.zeros((3, 1)), basis)
    base = sparse_cov(values, times, mean=mean, basis=basis, lam=0.1)
    scaled = sparse_cov([scale * v for v in values], times, mean=mean, basis=basis, lam=0.1)
    np.testing.assert_allclose(scaled.cov.coefs, scale**2 * base.cov.coefs, rtol=1e-8, atol=1e-10)


@settings(max_examples=25, deadline=None)
@given(seed=st.integers(0, 10_000))
def test_sparse_cov_ignores_curve_order(seed: int) -> None:
    times, values, _ = _simulate(n_curves=20, n_points=4, seed=seed)
    order = np.random.default_rng(seed).permutation(20)
    basis = BSpline(domain=UNIT, n_basis=3, order=2)
    mean = FData(np.zeros((3, 1)), basis)
    one = sparse_cov(values, times, mean=mean, basis=basis, lam=0.1).cov.coefs
    two = sparse_cov(
        [values[i] for i in order], [times[i] for i in order], mean=mean, basis=basis, lam=0.1
    ).cov.coefs
    np.testing.assert_allclose(one, two, rtol=1e-10, atol=1e-12)


# --------------------------------------------------------------------------- #
# eigenfunctions
# --------------------------------------------------------------------------- #


def test_harmonics_are_metric_orthonormal_and_signed() -> None:
    basis = BSpline(domain=UNIT, n_basis=6)
    rng = np.random.default_rng(6)
    root = rng.normal(size=(6, 6))
    cov = root @ root.T
    for lam in (0.0, 1e-3):
        values, coefs = _harmonics(cov, basis, basis, lam, LDO(2), 3)
        metric = np.asarray(basis.gram()) + lam * np.asarray(basis.penalty(2))
        np.testing.assert_allclose(coefs.T @ metric @ coefs, np.eye(3), atol=1e-10)
        assert np.all(np.diff(values) <= 0.0)
        upper = np.linalg.cholesky(metric).T
        assert np.all((upper @ coefs).sum(axis=0) > 0.0)


def test_harmonics_keep_at_most_the_basis_size() -> None:
    basis = Monomial(domain=UNIT, n_basis=2)
    values, coefs = _harmonics(np.eye(2), basis, basis, 0.0, LDO(2), 5)
    assert values.shape == (2,)
    assert coefs.shape == (2, 2)


def test_whitened_sum_signs_flips_negative_columns() -> None:
    metric = np.array([[2.0, 0.0], [0.0, 1.0]])
    harm = np.array([[-1.0, 0.5], [0.2, 0.5]])
    np.testing.assert_array_equal(_whitened_sum_signs(harm, metric), [-1.0, 1.0])


# --------------------------------------------------------------------------- #
# estimator
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def fitted() -> tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray]:
    times, values, scores = _simulate(n_curves=300, n_points=6, noise=0.3, seed=7)
    # sigma2 is fixed at the truth: with 300 curves the surface error on the
    # diagonal is larger than the noise variance (see the sparse_cov tests).
    model = PACE(n=2, basis=BSpline(domain=UNIT, n_basis=6), lam_cov=1e-6, sigma2=0.09)
    return model.fit(values, t=times), times, values, scores


def test_pace_recovers_eigenvalues_and_scores(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, _, _, truth = fitted
    np.testing.assert_allclose(model.values, [4.0, 0.49], rtol=0.35)
    np.testing.assert_allclose(model.varprop.sum(), 1.0)
    for k in range(2):
        corr = abs(np.corrcoef(model.scores[:, k], truth[:, k])[0, 1])
        assert corr > 0.8
    assert model.scores.shape == (300, 2)
    assert model.n_components_ == 2
    assert model.sigma2_ == 0.09
    assert model.cov_estimate_.sigma2 > 0.0


def test_pace_scores_are_the_conditional_expectation(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, times, values, _ = fitted
    lam = np.diag(np.maximum(model.values, 0.0))
    for i in (0, 17, 123):
        xi = model.harmonics(times[i])
        resid = values[i] - model.mean_fd(times[i])[:, 0]
        sigma = xi @ lam @ xi.T + model.sigma2_ * np.eye(len(times[i]))
        np.testing.assert_allclose(
            model.scores[i], lam @ xi.T @ np.linalg.solve(sigma, resid), rtol=1e-10
        )


def test_pace_transform_matches_fit(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, times, values, _ = fitted
    np.testing.assert_allclose(model.transform(values[:10], times[:10]), model.scores[:10])
    pairs = [np.column_stack([t, y]) for t, y in zip(times[:3], values[:3], strict=True)]
    np.testing.assert_allclose(model.transform(pairs), model.scores[:3])


def test_pace_accessors(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, _, _, _ = fitted
    assert model.cov.coefs.shape == (6, 6)
    assert model.mean_fd.n_curves == 1
    assert model.harmonics.n_curves == 2
    assert isinstance(model.cov_estimate_, SparseCov)


def test_pace_inverse_transform(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, _, _, _ = fitted
    grid = np.linspace(0.1, 0.9, 7)
    curves = model.inverse_transform(model.scores[:4], grid)
    want = model.mean_fd(grid)[:, 0] + model.scores[:4] @ model.harmonics(grid).T
    np.testing.assert_allclose(curves, want)
    np.testing.assert_allclose(model.inverse_transform(model.scores[0], grid), want[:1])
    with pytest.raises(ValueError, match="shape"):
        model.inverse_transform(np.ones((2, 3)), grid)


def test_pace_fit_transform_and_list_form() -> None:
    times, values, _ = _simulate(n_curves=60, seed=8)
    basis = BSpline(domain=UNIT, n_basis=4)
    scores = PACE(n=1, basis=basis).fit_transform(values, t=times)
    pairs = [np.column_stack([t, y]) for t, y in zip(times, values, strict=True)]
    np.testing.assert_allclose(PACE(n=1, basis=basis).fit(pairs).scores, scores)


def test_pace_separate_bases_and_penalties() -> None:
    times, values, _ = _simulate(n_curves=100, seed=9)
    model = PACE(
        n=2,
        basis=BSpline(domain=UNIT, n_basis=5),
        mean_basis=BSpline(domain=UNIT, n_basis=8),
        harmonic_basis=BSpline(domain=UNIT, n_basis=9),
        lam_mean=1e-4,
        lam_cov=1e-4,
        lam=1e-5,
        penalty=LDO(2),
    ).fit(values, t=times)
    assert model.harmonics.basis.n_basis == 9
    assert model.mean_fd.basis.n_basis == 8


def test_pace_sigma2_override() -> None:
    times, values, _ = _simulate(n_curves=50, seed=10)
    basis = BSpline(domain=UNIT, n_basis=4)
    small = PACE(n=1, basis=basis, sigma2=1e-6).fit(values, t=times)
    large = PACE(n=1, basis=basis, sigma2=1e3).fit(values, t=times)
    assert small.sigma2_ == 1e-6
    # A huge noise variance shrinks every score towards zero.
    assert np.abs(large.scores).max() < 0.1 * np.abs(small.scores).max()


def test_pace_negative_eigenvalues_do_not_enter_the_scores() -> None:
    times, values, _ = _simulate(n_curves=40, seed=11)
    model = PACE(n=1, basis=BSpline(domain=UNIT, n_basis=4)).fit(values, t=times)
    model.values_ = -np.abs(model.values_)
    np.testing.assert_allclose(model.transform(values, times), 0.0)


def test_pace_zero_spectrum_gives_zero_varprop() -> None:
    times = [np.array([0.0, 0.5, 1.0])] * 6
    values = [np.array([1.0, 2.0, 3.0])] * 6
    basis = BSpline(domain=UNIT, n_basis=2, order=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        model = PACE(n=1, basis=basis, mean_basis=basis).fit(values, t=times)
    np.testing.assert_allclose(model.varprop, 0.0)


@pytest.mark.parametrize(
    ("params", "match"),
    [
        ({"n": 0}, "n must be"),
        ({"n": 1.5}, "n must be"),
        ({"lam": -1.0}, "lam must"),
        ({"lam_cov": -1.0}, "lam_cov"),
        ({"lam_mean": -1.0}, "lam_mean"),
        ({"sigma2": 0.0}, "sigma2"),
    ],
)
def test_pace_rejects_bad_parameters(params: dict[str, object], match: str) -> None:
    times, values, _ = _simulate(n_curves=20)
    with pytest.raises(ValueError, match=match):
        PACE(**params).fit(values, t=times)  # type: ignore[arg-type]


def test_pace_unfitted_raises() -> None:
    model = PACE()
    for attr in ("mean_fd", "cov", "harmonics", "values", "varprop", "scores"):
        with pytest.raises(NotFittedError):
            getattr(model, attr)
    with pytest.raises(NotFittedError):
        model.transform([np.array([1.0])], [np.array([0.0])])


def test_pace_estimator_api(
    fitted: tuple[PACE, list[np.ndarray], list[np.ndarray], np.ndarray],
) -> None:
    model, times, values, _ = fitted
    params = model.get_params()
    assert params["n"] == 2
    assert params["lam_cov"] == 1e-6
    fresh = clone(model)
    assert not hasattr(fresh, "scores_")
    assert fresh.set_params(n=1).n == 1
    restored = pickle.loads(pickle.dumps(model))
    np.testing.assert_allclose(restored.transform(values[:5], times[:5]), model.scores[:5])


@settings(max_examples=15, deadline=None)
@given(scale=st.floats(0.2, 5.0))
def test_pace_scores_scale_linearly(scale: float) -> None:
    times, values, _ = _simulate(n_curves=40, seed=12)
    basis = BSpline(domain=UNIT, n_basis=4)
    base = PACE(n=2, basis=basis, sigma2=0.05).fit(values, t=times)
    scaled = PACE(n=2, basis=basis, sigma2=0.05 * scale**2).fit(
        [scale * v for v in values], t=times
    )
    np.testing.assert_allclose(scaled.scores, scale * base.scores, rtol=1e-7, atol=1e-9)
