"""Unit tests for :mod:`fabel.decomposition` (FPCA, FCCA)."""

from __future__ import annotations

from typing import Any, cast

import matplotlib
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

from fabel import LDO, BSpline, FData, Fourier
from fabel.decomposition import FCCA, FPCA, _block_diagonal, _varimax_rotation
from fabel.smoothing import Smoother

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def make_fd(n_curves: int = 30, n_basis: int = 8, seed: int = 0) -> FData:
    """Return random B-spline curves with a decaying covariance spectrum."""
    rng = np.random.default_rng(seed)
    scale = 1.0 / (1.0 + np.arange(n_basis))[:, None]
    coefs = 3.0 + scale * rng.standard_normal((n_basis, n_curves))
    return FData(coefs, BSpline(domain=(0.0, 1.0), n_basis=n_basis))


def make_multivariate(n_curves: int = 25, n_basis: int = 7, seed: int = 1) -> FData:
    rng = np.random.default_rng(seed)
    coefs = rng.standard_normal((n_basis, n_curves, 2))
    coefs[:, :, 1] += 0.5 * coefs[:, :, 0]
    return FData(coefs, Fourier(domain=(0.0, 1.0), n_basis=n_basis))


def varimax_criterion(loadings: np.ndarray) -> float:
    squared = loadings**2
    n = loadings.shape[0]
    return float(np.sum(np.sum(squared**2, axis=0) - np.sum(squared, axis=0) ** 2 / n))


# --------------------------------------------------------------------------- #
# FPCA: the decomposition itself
# --------------------------------------------------------------------------- #


def test_unpenalised_harmonics_are_l2_orthonormal() -> None:
    fd = make_fd()
    pca = FPCA(n=4).fit(fd)
    harm = np.asarray(pca.harmonics.coefs)
    gram = np.asarray(fd.basis.gram())
    np.testing.assert_allclose(harm.T @ gram @ harm, np.eye(4), atol=1e-10)


def test_penalised_harmonics_are_orthonormal_in_the_penalised_metric() -> None:
    fd = make_fd()
    lam = 1e-3
    pca = FPCA(n=3, lam=lam).fit(fd)
    harm = np.asarray(pca.harmonics.coefs)
    metric = np.asarray(fd.basis.gram()) + lam * np.asarray(fd.basis.penalty(2))
    np.testing.assert_allclose(harm.T @ metric @ harm, np.eye(3), atol=1e-10)
    assert pca.lam_ == lam


def test_values_are_the_full_descending_spectrum() -> None:
    fd = make_fd(n_basis=9)
    pca = FPCA(n=2).fit(fd)
    assert pca.values.shape == (9,)
    assert np.all(np.diff(pca.values) <= 1e-12)
    np.testing.assert_allclose(pca.varprop, pca.values[:2] / pca.values.sum())


def test_unpenalised_values_equal_score_variances() -> None:
    fd = make_fd()
    pca = FPCA(n=3).fit(fd)
    np.testing.assert_allclose((pca.scores**2).mean(axis=0), pca.values[:3], rtol=1e-10)


def test_total_variance_is_the_integrated_pointwise_variance() -> None:
    fd = make_fd()
    pca = FPCA(n=2).fit(fd)
    coefs = np.asarray(fd.coefs)
    centred = coefs - coefs.mean(axis=1, keepdims=True)
    gram = np.asarray(fd.basis.gram())
    total = np.trace(centred.T @ gram @ centred) / coefs.shape[1]
    assert pca.values.sum() == pytest.approx(total, rel=1e-10)


def test_harmonics_have_a_positive_coefficient_sum() -> None:
    pca = FPCA(n=5).fit(make_fd(seed=3))
    assert np.all(np.asarray(pca.harmonics.coefs).sum(axis=0) > 0.0)


def test_scores_equal_transform_of_the_training_curves() -> None:
    fd = make_fd()
    pca = FPCA(n=3).fit(fd)
    np.testing.assert_allclose(pca.transform(fd), pca.scores, atol=1e-12)
    np.testing.assert_allclose(pca.fit_transform(fd), pca.scores, atol=1e-12)


def test_full_rank_inverse_transform_reconstructs_the_curves() -> None:
    fd = make_fd(n_basis=6)
    pca = FPCA(n=6).fit(fd)
    rebuilt = pca.inverse_transform(pca.scores)
    np.testing.assert_allclose(rebuilt.coefs, fd.coefs, atol=1e-10)


def test_penalised_inverse_transform_is_the_projection_on_the_harmonics() -> None:
    fd = make_fd(n_basis=6)
    pca = FPCA(n=6, lam=1e-2).fit(fd)
    rebuilt = pca.inverse_transform(pca.scores)
    np.testing.assert_allclose(rebuilt.coefs, fd.coefs, atol=1e-9)


def test_inverse_transform_accepts_a_single_score_vector() -> None:
    pca = FPCA(n=2).fit(make_fd())
    rebuilt = pca.inverse_transform(pca.scores[0])
    assert rebuilt.n_curves == 1


def test_uncentred_fit_has_a_zero_mean() -> None:
    fd = make_fd()
    pca = FPCA(n=2, center=False).fit(fd)
    assert np.all(np.asarray(pca.mean_fd.coefs) == 0.0)
    # The uncentred leading harmonic follows the (large) mean level.
    assert pca.varprop[0] > 0.9


def test_mean_fd_is_the_sample_mean() -> None:
    fd = make_fd()
    pca = FPCA(n=2).fit(fd)
    np.testing.assert_allclose(
        np.asarray(pca.mean_fd.coefs)[:, 0], np.asarray(fd.coefs).mean(axis=1)
    )


def test_integer_like_n_is_accepted() -> None:
    pca = FPCA(n=cast(Any, np.int64(2))).fit(make_fd())
    assert pca.scores.shape == (30, 2)


def test_n_is_capped_at_the_basis_size() -> None:
    pca = FPCA(n=50).fit(make_fd(n_basis=5))
    assert pca.n_components_ == 5


def test_harmonic_operator_penalty() -> None:
    rng = np.random.default_rng(4)
    fd = FData(rng.standard_normal((7, 20)), Fourier(domain=(0.0, 1.0), n_basis=7))
    pca = FPCA(n=2, lam=1e-4, penalty=LDO.harmonic(1.0)).fit(fd)
    assert pca.harmonics.n_curves == 2


def test_constant_curves_give_zero_variance_proportions() -> None:
    coefs = np.ones((5, 4))
    pca = FPCA(n=2).fit(FData(coefs, BSpline(domain=(0.0, 1.0), n_basis=5)))
    np.testing.assert_array_equal(pca.varprop, np.zeros(2))


# --------------------------------------------------------------------------- #
# FPCA: parameters and errors
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("n", [0, -1, 2.5, True, "3"])
def test_invalid_n_raises(n: Any) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        FPCA(n=n).fit(make_fd())


def test_negative_lambda_raises() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        FPCA(lam=-1.0).fit(make_fd())


def test_unknown_lambda_string_raises() -> None:
    with pytest.raises(ValueError, match="'gcv'"):
        FPCA(lam="aic").fit(make_fd())


def test_gcv_needs_three_curves() -> None:
    with pytest.raises(ValueError, match="three curves"):
        FPCA(lam="gcv").fit(make_fd(n_curves=2))


def test_gcv_picks_a_lambda_from_the_grid() -> None:
    fd = make_fd(n_curves=12, n_basis=6)
    pca = FPCA(n=2, lam="gcv").fit(fd)
    exponent = np.log10(pca.lam_)
    assert exponent == pytest.approx(round(exponent))
    assert -8 <= round(exponent) <= 8


def test_gcv_minimises_the_leave_one_out_reconstruction_error() -> None:
    fd = make_fd(n_curves=8, n_basis=6, seed=13)
    pca = FPCA(n=2, lam="gcv").fit(fd)
    coefs = np.asarray(fd.coefs)
    errors = [pca._reconstruction_sse(coefs, fd.basis, 10.0**e) for e in range(-8, 9)]
    assert pca.lam_ == 10.0 ** (int(np.argmin(errors)) - 8)
    # A held-out curve is never reconstructed better than exactly.
    assert min(errors) > 0.0


def test_methods_before_fit_raise() -> None:
    pca = FPCA()
    with pytest.raises(NotFittedError):
        _ = pca.harmonics
    for name in ("values", "scores", "scores_by_var", "varprop", "mean_fd", "rotation"):
        with pytest.raises(NotFittedError):
            getattr(pca, name)
    with pytest.raises(NotFittedError):
        pca.transform(make_fd())


# --------------------------------------------------------------------------- #
# FPCA: plain arrays and pipelines
# --------------------------------------------------------------------------- #


def test_plain_array_uses_the_default_spline() -> None:
    rng = np.random.default_rng(6)
    x = rng.standard_normal((20, 6))
    pca = FPCA(n=2).fit(x)
    assert pca.harmonics.basis == BSpline(domain=(0.0, 1.0), n_basis=6, order=4)
    assert pca.n_features_in_ == 6
    np.testing.assert_allclose(pca.transform(x), pca.scores, atol=1e-12)


def test_plain_array_with_an_explicit_basis() -> None:
    basis = Fourier(domain=(0.0, 2.0), n_basis=5)
    x = np.random.default_rng(7).standard_normal((15, 5))
    pca = FPCA(n=2, basis=basis).fit(x)
    assert pca.harmonics.basis is basis


def test_plain_array_basis_size_mismatch_raises() -> None:
    with pytest.raises(ValueError, match="columns"):
        FPCA(basis=BSpline(n_basis=7)).fit(np.zeros((5, 6)))


def test_transform_rejects_a_different_basis_size() -> None:
    pca = FPCA(n=2).fit(make_fd(n_basis=8))
    with pytest.raises(ValueError, match="basis size"):
        pca.transform(make_fd(n_basis=6))


def test_pipeline_after_smoother() -> None:
    rng = np.random.default_rng(8)
    t = np.linspace(0.0, 1.0, 30)
    labels = np.repeat([0, 1], 20)
    clean = np.sin(2 * np.pi * t)[None, :] * (1.0 + labels[:, None])
    curves = clean + 0.05 * rng.standard_normal(clean.shape)
    basis = BSpline(domain=(0.0, 1.0), n_basis=10)
    pipe = Pipeline(
        [
            ("smooth", Smoother(basis, t=t, lam=1e-6)),
            ("fpca", FPCA(n=2, basis=basis)),
            ("clf", LogisticRegression()),
        ]
    )
    pipe.fit(curves, labels)
    assert pipe.score(curves, labels) == 1.0


def test_get_params_round_trip() -> None:
    pca = FPCA(n=3, lam=1e-2, center=False)
    assert pca.get_params()["n"] == 3
    assert FPCA(**pca.get_params()).get_params() == pca.get_params()


# --------------------------------------------------------------------------- #
# FPCA: multivariate curves
# --------------------------------------------------------------------------- #


def test_multivariate_shapes() -> None:
    fd = make_multivariate()
    pca = FPCA(n=3).fit(fd)
    assert np.asarray(pca.harmonics.coefs).shape == (7, 3, 2)
    assert pca.scores_by_var.shape == (25, 3, 2)
    assert pca.scores.shape == (25, 3)
    assert pca.values.shape == (14,)
    assert np.asarray(pca.mean_fd.coefs).shape == (7, 1, 2)
    np.testing.assert_allclose(pca.scores, pca.scores_by_var.sum(axis=2))


def test_multivariate_harmonics_are_jointly_orthonormal() -> None:
    fd = make_multivariate()
    pca = FPCA(n=4).fit(fd)
    harm = np.asarray(pca.harmonics.coefs)
    gram = np.asarray(fd.basis.gram())
    inner = sum(harm[:, :, k].T @ gram @ harm[:, :, k] for k in range(2))
    np.testing.assert_allclose(inner, np.eye(4), atol=1e-10)


def test_multivariate_transform_and_inverse() -> None:
    fd = make_multivariate(n_basis=5)
    pca = FPCA(n=10).fit(fd)
    np.testing.assert_allclose(pca.transform(fd), pca.scores, atol=1e-10)
    np.testing.assert_allclose(pca.inverse_transform(pca.scores).coefs, fd.coefs, atol=1e-9)


def test_multivariate_transform_rejects_univariate_input() -> None:
    pca = FPCA(n=2).fit(make_multivariate())
    with pytest.raises(ValueError, match="variables"):
        pca.transform(FData(np.zeros((7, 3)), Fourier(domain=(0.0, 1.0), n_basis=7)))


def test_multivariate_gcv_runs() -> None:
    pca = FPCA(n=1, lam="gcv").fit(make_multivariate(n_curves=6, n_basis=5))
    assert pca.lam_ > 0.0


def test_multivariate_rotation() -> None:
    fd = make_multivariate()
    pca = FPCA(n=3).fit(fd)
    rotated = pca.rotate()
    assert np.asarray(rotated.harmonics.coefs).shape == (7, 3, 2)
    np.testing.assert_allclose(rotated.scores, pca.scores @ rotated.rotation, atol=1e-10)


# --------------------------------------------------------------------------- #
# rotation
# --------------------------------------------------------------------------- #


def test_rotation_is_orthogonal_and_preserves_the_span() -> None:
    fd = make_fd()
    pca = FPCA(n=3).fit(fd)
    rotated = pca.rotate("varimax")
    rot = rotated.rotation
    assert rot is not None
    np.testing.assert_allclose(rot.T @ rot, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(
        rotated.harmonics.coefs, np.asarray(pca.harmonics.coefs) @ rot, atol=1e-12
    )
    np.testing.assert_allclose(rotated.scores, pca.scores @ rot, atol=1e-10)
    # Total retained variance is invariant under rotation.
    assert rotated.values.sum() == pytest.approx(pca.values[:3].sum(), rel=1e-10)
    np.testing.assert_allclose(rotated.varprop, rotated.values / pca.values.sum())
    assert pca.rotation is None


def test_rotation_increases_the_varimax_criterion() -> None:
    fd = make_fd(seed=9)
    pca = FPCA(n=3).fit(fd)
    rotated = pca.rotate()
    grid = np.linspace(0.0, 1.0, 501)
    before = varimax_criterion(np.asarray(pca.harmonics(grid)))
    after = varimax_criterion(np.asarray(rotated.harmonics(grid)))
    assert after >= before - 1e-12


def test_rotation_keeps_the_sign_rule() -> None:
    rotated = FPCA(n=4).fit(make_fd(seed=10)).rotate()
    assert np.all(np.asarray(rotated.harmonics.coefs).sum(axis=0) > 0.0)


def test_rotating_one_component_is_the_identity() -> None:
    pca = FPCA(n=1).fit(make_fd())
    rotated = pca.rotate()
    assert rotated.rotation is not None
    np.testing.assert_allclose(rotated.rotation, np.eye(1))


def test_unknown_rotation_raises() -> None:
    with pytest.raises(ValueError, match="varimax"):
        FPCA(n=2).fit(make_fd()).rotate("promax")


def test_uncentred_constant_rotation_has_zero_varprop() -> None:
    pca = FPCA(n=2, center=True).fit(FData(np.ones((5, 4)), BSpline(n_basis=5)))
    np.testing.assert_array_equal(pca.rotate().varprop, np.zeros(2))


@settings(max_examples=25, deadline=None)
@given(
    st.integers(min_value=2, max_value=4),
    st.integers(min_value=0, max_value=2**31 - 1),
)
def test_varimax_is_orthogonal_and_stationary(n_comp: int, seed: int) -> None:
    loadings = np.random.default_rng(seed).standard_normal((60, n_comp))
    rot = _varimax_rotation(loadings)
    np.testing.assert_allclose(rot.T @ rot, np.eye(n_comp), atol=1e-10)
    rotated = loadings @ rot
    assert varimax_criterion(rotated) >= varimax_criterion(loadings) - 1e-10
    # Stationarity on the rotation group: the skew part of L' dV/dL vanishes.
    squared = rotated**2
    grad = 4.0 * (rotated**3 - rotated * squared.sum(axis=0) / rotated.shape[0])
    moment = rotated.T @ grad
    assert np.max(np.abs(moment - moment.T)) <= 1e-7 * max(1.0, np.max(np.abs(moment)))


def test_block_diagonal() -> None:
    block = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(_block_diagonal(block, 1), block)
    expected = np.zeros((4, 4))
    expected[:2, :2] = block
    expected[2:, 2:] = block
    np.testing.assert_array_equal(_block_diagonal(block, 2), expected)


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #


def test_plot_creates_one_panel_per_harmonic() -> None:
    pca = FPCA(n=3).fit(make_fd())
    axes = pca.plot(n_points=51)
    assert len(axes) == 3
    lines = axes[0].get_lines()
    assert len(lines) == 3
    grid = np.linspace(0.0, 1.0, 51)
    mean = np.asarray(pca.mean_fd(grid)).ravel()
    np.testing.assert_allclose(lines[0].get_xydata()[:, 1], mean)
    spread = np.sqrt((pca.scores[:, 0] ** 2).mean())
    harm = np.asarray(pca.harmonics(grid))[:, 0]
    np.testing.assert_allclose(lines[1].get_xydata()[:, 1], mean + spread * harm)
    plt.close("all")


def test_plot_on_given_axes() -> None:
    pca = FPCA(n=2).fit(make_fd())
    _, grid_axes = plt.subplots(1, 2)
    assert len(pca.plot(grid_axes)) == 2
    assert len(pca.plot(list(grid_axes))) == 2
    single = FPCA(n=1).fit(make_fd())
    _, ax = plt.subplots()
    assert single.plot(ax) == [ax]
    plt.close("all")


def test_plot_multivariate() -> None:
    pca = FPCA(n=2).fit(make_multivariate())
    axes = pca.plot(n_points=21)
    assert len(axes[0].get_lines()) == 6
    plt.close("all")


def test_plot_with_too_few_axes_raises() -> None:
    pca = FPCA(n=3).fit(make_fd())
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="axes"):
        pca.plot(ax)
    plt.close("all")


# --------------------------------------------------------------------------- #
# FCCA
# --------------------------------------------------------------------------- #


def cca_pair(n_curves: int = 40, noise: float = 0.3, seed: int = 11) -> tuple[FData, FData]:
    rng = np.random.default_rng(seed)
    basis = Fourier(domain=(0.0, 1.0), n_basis=7)
    shared = rng.standard_normal((7, n_curves))
    x = FData(shared + noise * rng.standard_normal((7, n_curves)), basis)
    y = FData(-shared + noise * rng.standard_normal((7, n_curves)), basis)
    return x, y


def test_cca_correlations_are_descending_and_bounded() -> None:
    x, y = cca_pair()
    cca = FCCA(n=3, lam1=1e-3, lam2=1e-3).fit(x, y)
    corr = np.asarray(cca.correlations)
    assert corr.shape == (7,)
    assert np.all(np.diff(corr) <= 1e-12)
    assert np.all(corr <= 1.0 + 1e-12)
    assert corr[0] > 0.8


def test_cca_weights_have_unit_norm_and_pairs_correlate_positively() -> None:
    x, y = cca_pair()
    cca = FCCA(n=3, lam1=1e-3, lam2=1e-3).fit(x, y)
    a = np.asarray(cca.weights1.coefs)
    b = np.asarray(cca.weights2.coefs)
    np.testing.assert_allclose(np.einsum("ij,ik,kj->j", a, np.asarray(x.basis.gram()), a), 1.0)
    np.testing.assert_allclose(np.einsum("ij,ik,kj->j", b, np.asarray(y.basis.gram()), b), 1.0)
    assert np.all(a.sum(axis=0) > 0.0)
    for j in range(3):
        assert np.corrcoef(cca.scores1[:, j], cca.scores2[:, j])[0, 1] > 0.0


def test_cca_scores_are_the_penalised_canonical_correlations() -> None:
    # Unpenalised, the empirical correlation of each score pair equals the
    # canonical correlation (sample covariance uses 1/N throughout).
    x, y = cca_pair(n_curves=200)
    cca = FCCA(n=2).fit(x, y)
    for j in range(2):
        r = np.corrcoef(cca.scores1[:, j], cca.scores2[:, j])[0, 1]
        assert r == pytest.approx(float(cca.correlations[j]), rel=1e-10)


def test_cca_transform_matches_the_training_scores() -> None:
    x, y = cca_pair()
    cca = FCCA(n=2, lam1=1e-2, lam2=1e-2).fit(x, y)
    s1, s2 = cca.transform(x, y)
    np.testing.assert_allclose(s1, cca.scores1, atol=1e-12)
    np.testing.assert_allclose(s2, cca.scores2, atol=1e-12)


def test_cca_without_centering() -> None:
    x, y = cca_pair()
    cca = FCCA(n=1, lam1=1e-2, lam2=1e-2, center=False).fit(x, y)
    s1, _ = cca.transform(x, y)
    np.testing.assert_allclose(s1, cca.scores1, atol=1e-12)
    np.testing.assert_array_equal(cca.means_[0], np.zeros((7, 1)))


def test_cca_accepts_different_bases() -> None:
    x, _ = cca_pair()
    rng = np.random.default_rng(12)
    y = FData(rng.standard_normal((5, 40)), BSpline(domain=(0.0, 1.0), n_basis=5))
    cca = FCCA(n=4, lam1=1e-2, lam2=1e-2).fit(x, y)
    assert cca.correlations.shape == (5,)
    assert cca.n_components_ == 4


def test_cca_rejects_non_fdata() -> None:
    x, _ = cca_pair()
    with pytest.raises(TypeError, match="FData"):
        FCCA().fit(x, np.zeros((7, 40)))  # type: ignore[arg-type]


def test_cca_rejects_mismatched_curve_counts() -> None:
    x, _ = cca_pair(n_curves=40)
    _, y = cca_pair(n_curves=30)
    with pytest.raises(ValueError, match="must match"):
        FCCA().fit(x, y)


def test_cca_rejects_multivariate_curves() -> None:
    fd = make_multivariate(n_curves=40)
    x, _ = cca_pair()
    with pytest.raises(ValueError, match="univariate"):
        FCCA().fit(x, fd)


@pytest.mark.parametrize("n", [0, 1.5, False])
def test_cca_invalid_n_raises(n: Any) -> None:
    x, y = cca_pair()
    with pytest.raises(ValueError, match="positive integer"):
        FCCA(n=n).fit(x, y)


def test_cca_negative_lambda_raises() -> None:
    x, y = cca_pair()
    with pytest.raises(ValueError, match="non-negative"):
        FCCA(lam2=-1.0).fit(x, y)


def test_cca_transform_rejects_a_different_basis_size() -> None:
    x, y = cca_pair()
    cca = FCCA(n=1, lam1=1e-2, lam2=1e-2).fit(x, y)
    small = FData(np.zeros((5, 3)), Fourier(domain=(0.0, 1.0), n_basis=5))
    with pytest.raises(ValueError, match="basis sizes"):
        cca.transform(small, y)


def test_cca_accessors_before_fit_raise() -> None:
    cca = FCCA()
    for name in ("weights1", "weights2", "correlations", "scores1", "scores2"):
        with pytest.raises(NotFittedError):
            getattr(cca, name)


def test_cca_tags_require_a_target() -> None:
    assert FCCA().__sklearn_tags__().target_tags.required is True
    assert FPCA().__sklearn_tags__().target_tags.required is False


# --------------------------------------------------------------------------- #
# n accepts any integer-like value (typing.SupportsIndex)
# --------------------------------------------------------------------------- #


class _IndexLike:
    """A non-int object that still implements ``__index__``."""

    def __init__(self, value: int) -> None:
        self.value = value

    def __index__(self) -> int:
        return self.value


@pytest.mark.parametrize("n", [np.int64(3), np.int32(3), np.uint8(3), _IndexLike(3)])
def test_fpca_accepts_integer_like_n(n: Any) -> None:
    pca = FPCA(n=n).fit(make_fd())
    assert pca.n_components_ == 3
    assert pca.harmonics.n_curves == 3
    assert pca.get_params()["n"] is n


@pytest.mark.parametrize("n", [np.int64(2), _IndexLike(2)])
def test_cca_accepts_integer_like_n(n: Any) -> None:
    x, y = cca_pair()
    cca = FCCA(n=n, lam1=1e-4, lam2=1e-4).fit(x, y)
    assert cca.n_components_ == 2
    assert cca.get_params()["n"] is n


@pytest.mark.parametrize("n", [np.float64(2.0), 2.0, np.True_, np.int64(0), _IndexLike(-1)])
def test_integer_like_n_rejects_non_integers(n: Any) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        FPCA(n=n).fit(make_fd())
    x, y = cca_pair()
    with pytest.raises(ValueError, match="positive integer"):
        FCCA(n=n).fit(x, y)
