"""Unit tests for :mod:`fdatools.smoothing` (smooth, SmoothResult, Smoother, lambda helpers)."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.exceptions import NotFittedError
from sklearn.pipeline import Pipeline

from fdatools import LDO, BSpline, FData, Fourier, Monomial
from fdatools import smoothing as sm
from fdatools.smoothing import (
    Smoother,
    SmoothResult,
    df_to_lambda,
    gcv_curve,
    lambda_to_df,
    smooth,
)

DOMAIN = (0.0, 1.0)
T = np.linspace(0.0, 1.0, 41)
RNG = np.random.default_rng(20260927)


def spline(n_basis: int = 12) -> BSpline:
    return BSpline(domain=DOMAIN, n_basis=n_basis)


def noisy(n_curves: int = 3, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    signal = np.sin(2 * np.pi * T)[:, None] * (1.0 + np.arange(n_curves))[None, :]
    out: np.ndarray = signal + 0.1 * rng.standard_normal((T.size, n_curves))
    return out


def normal_equations(
    y: np.ndarray, basis: BSpline, lam: float, weights: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Return the reference coefficients and hat matrix of a penalised fit."""
    phi = np.asarray(basis(T))
    w = np.ones(T.size) if weights is None else weights
    lhs = phi.T @ (w[:, None] * phi) + lam * np.asarray(basis.penalty(2))
    y2c = np.linalg.solve(lhs, phi.T * w)
    return y2c @ y, phi @ y2c


# --------------------------------------------------------------------------- #
# the linear fit
# --------------------------------------------------------------------------- #


def test_fixed_lambda_solves_the_penalised_normal_equations() -> None:
    y = noisy()
    basis = spline()
    result = smooth(y, T, basis=basis, lam=1e-4)
    coefs, hat = normal_equations(y, basis, 1e-4)
    np.testing.assert_allclose(result.fd.coefs, coefs, rtol=1e-10, atol=1e-12)
    assert result.df == pytest.approx(np.trace(hat), rel=1e-10)
    resid = y - hat @ y
    assert result.sse == pytest.approx(float(np.sum(resid**2)), rel=1e-10)
    n = T.size
    gcv = (np.sum(resid**2, axis=0) / n) / (1.0 - np.trace(hat) / n) ** 2
    np.testing.assert_allclose(result.gcv, gcv, rtol=1e-10)
    np.testing.assert_allclose(result.penalty_matrix, basis.penalty(2))
    assert result.lam == 1e-4
    assert result.constraint is None
    assert result.beta is None


def test_weights_enter_the_normal_equations() -> None:
    y = noisy(n_curves=1)
    weights = np.linspace(0.5, 2.0, T.size)
    basis = spline()
    result = smooth(y, T, basis=basis, lam=1e-3, weights=weights)
    coefs, _ = normal_equations(y, basis, 1e-3, weights)
    np.testing.assert_allclose(result.fd.coefs, coefs, rtol=1e-10, atol=1e-12)
    unweighted = smooth(y, T, basis=basis, lam=1e-3)
    assert not np.allclose(result.fd.coefs, unweighted.fd.coefs)


def test_weight_shape_is_validated() -> None:
    with pytest.raises(ValueError, match="weights must have shape"):
        smooth(noisy(), T, basis=spline(), lam=1.0, weights=np.ones(3))


def test_result_call_evaluates_the_curves() -> None:
    result = smooth(noisy(), T, basis=spline(), lam=1e-6)
    np.testing.assert_allclose(result(T), result.fd(T))
    np.testing.assert_allclose(result(T, 1), result.fd(T, 1))


def test_single_curve_has_a_scalar_gcv() -> None:
    result = smooth(np.sin(T), T, basis=spline(), lam=1e-6)
    assert np.ndim(result.gcv) == 0
    assert result.y2c_map.shape == (12, T.size)


def test_multivariate_observations_keep_their_curve_axes() -> None:
    y = RNG.standard_normal((T.size, 4, 2))
    result = smooth(y, T, basis=spline(), lam=1e-2)
    assert np.asarray(result.fd.coefs).shape == (12, 4, 2)
    assert np.shape(result.gcv) == (4, 2)
    flat = smooth(y.reshape(T.size, 8), T, basis=spline(), lam=1e-2)
    np.testing.assert_allclose(np.asarray(result.fd.coefs).reshape(12, 8), flat.fd.coefs)


def test_interpolating_fit_has_infinite_gcv() -> None:
    t = np.linspace(0.0, 1.0, 8)
    result = smooth(np.cos(t), t, basis=BSpline(domain=DOMAIN, n_basis=8), lam=0.0)
    assert result.df == pytest.approx(8.0)
    assert np.isinf(result.gcv)


def test_operator_penalty() -> None:
    t = np.linspace(0.0, 1.0, 50)
    y = np.sin(2 * np.pi * t)
    basis = Fourier(domain=DOMAIN, n_basis=9)
    # The harmonic accelerator annihilates the fundamental, so a huge penalty
    # keeps a pure sine exactly.
    result = smooth(y, t, basis=basis, lam=1e8, penalty=LDO.harmonic(1.0))
    np.testing.assert_allclose(result.fd(t)[:, 0], y, atol=1e-6)


def test_default_basis_is_a_cubic_spline_on_the_data_range() -> None:
    t = np.linspace(2.0, 5.0, 30)
    result = smooth(np.exp(-t), t, lam=1e-6)
    basis = result.fd.basis
    assert isinstance(basis, BSpline)
    assert basis.domain == (2.0, 5.0)
    assert basis.n_basis == 32


def test_default_basis_is_capped() -> None:
    t = np.linspace(0.0, 1.0, 200)
    assert smooth(t, t, lam=1.0).fd.basis.n_basis == 42


def test_default_basis_needs_a_non_degenerate_range() -> None:
    with pytest.raises(ValueError, match="non-degenerate"):
        smooth(np.ones(4), np.zeros(4), lam=1.0)


def test_t_must_be_one_dimensional() -> None:
    with pytest.raises(ValueError, match="one-dimensional"):
        smooth(np.ones(4), np.zeros((2, 2)), basis=spline(), lam=1.0)


def test_y_rows_must_match_t() -> None:
    with pytest.raises(ValueError, match="rows"):
        smooth(np.ones(5), T, basis=spline(), lam=1.0)


# --------------------------------------------------------------------------- #
# choosing lambda
# --------------------------------------------------------------------------- #


def test_negative_lambda_raises() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        smooth(noisy(), T, basis=spline(), lam=-1.0)


def test_unknown_lambda_string_raises() -> None:
    with pytest.raises(ValueError, match="'gcv'"):
        smooth(noisy(), T, basis=spline(), lam="aic")


def test_gcv_lambda_minimises_the_mean_gcv() -> None:
    y = noisy(seed=4)
    basis = spline(15)
    result = smooth(y, T, basis=basis, lam="gcv")
    best = float(np.mean(result.gcv))
    for factor in (0.5, 2.0):
        other = smooth(y, T, basis=basis, lam=result.lam * factor)
        assert best <= float(np.mean(other.gcv)) + 1e-12
    lo, hi = sm.LOG10_LAMBDA_RANGE
    assert lo <= np.log10(result.lam) <= hi


def test_gcv_is_the_default() -> None:
    y = noisy(seed=5)
    assert smooth(y, T, basis=spline()).lam == smooth(y, T, basis=spline(), lam="gcv").lam


def test_df_string_and_keyword_agree() -> None:
    y = noisy()
    basis = spline()
    by_string = smooth(y, T, basis=basis, lam="df=6")
    by_keyword = smooth(y, T, basis=basis, df=6.0, lam=1e9)
    assert by_string.df == pytest.approx(6.0, abs=1e-9)
    assert by_keyword.lam == pytest.approx(by_string.lam, rel=1e-12)


def test_unreachable_df_raises() -> None:
    with pytest.raises(ValueError, match="exceeds"):
        smooth(noisy(), T, basis=spline(), df=50.0)
    with pytest.raises(ValueError, match="below"):
        smooth(noisy(), T, basis=spline(), df=0.5)


def test_lambda_to_df_matches_the_trace_of_the_hat_matrix() -> None:
    basis = spline()
    _, hat = normal_equations(np.zeros(T.size), basis, 1e-3)
    assert lambda_to_df(T, basis, 1e-3) == pytest.approx(np.trace(hat), rel=1e-10)


def test_lambda_to_df_with_weights_and_operator() -> None:
    basis = spline()
    weighted = lambda_to_df(T, basis, 1e-3, weights=np.linspace(1.0, 2.0, T.size))
    first = lambda_to_df(T, basis, 1e-3, penalty=LDO(1))
    assert 2.0 < weighted < 12.0
    assert 1.0 < first < 12.0


@settings(max_examples=20, deadline=None)
@given(st.floats(min_value=2.5, max_value=11.5))
def test_df_to_lambda_inverts_lambda_to_df(target: float) -> None:
    basis = spline()
    lam = df_to_lambda(T, basis, target)
    assert lambda_to_df(T, basis, lam) == pytest.approx(target, abs=1e-8)


@settings(max_examples=20, deadline=None)
@given(st.floats(min_value=-8.0, max_value=6.0), st.floats(min_value=0.1, max_value=2.0))
def test_df_decreases_with_lambda(log_lam: float, step: float) -> None:
    basis = spline()
    lam = 10.0**log_lam
    assert lambda_to_df(T, basis, lam * 10.0**step) <= lambda_to_df(T, basis, lam) + 1e-9


def test_gcv_curve_matches_individual_fits() -> None:
    y = noisy()
    basis = spline()
    lambdas = [1e-6, 1e-3, 1.0]
    scores = gcv_curve(y, T, basis, lambdas)
    assert scores.shape == (3, 3)
    for row, lam in zip(scores, lambdas, strict=True):
        np.testing.assert_allclose(row, smooth(y, T, basis=basis, lam=lam).gcv, rtol=1e-8)


def test_gcv_curve_is_infinite_when_the_fit_interpolates() -> None:
    t = np.linspace(0.0, 1.0, 8)
    scores = gcv_curve(np.cos(t), t, BSpline(domain=DOMAIN, n_basis=8), [0.0])
    assert np.isinf(scores[0])


def test_golden_section_finds_a_parabola_minimum() -> None:
    assert sm._golden_section(lambda x: (x - 0.3) ** 2, -1.0, 2.0) == pytest.approx(0.3, abs=1e-7)
    assert sm._golden_section(lambda x: x, 1.0, 1.0) == 1.0


def test_mean_gcv_maps_non_finite_values_to_infinity() -> None:
    class _Stub:
        _xp = np

        @staticmethod
        def gcv(lam: float) -> np.ndarray:
            return np.array([np.nan])

    assert sm._mean_gcv(_Stub(), 0.0) == np.inf  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# irregular designs
# --------------------------------------------------------------------------- #


def irregular() -> tuple[list[np.ndarray], list[np.ndarray]]:
    rng = np.random.default_rng(7)
    grids = [np.sort(rng.uniform(0.0, 1.0, size)) for size in (25, 31, 40)]
    values = [np.sin(2 * np.pi * g) + 0.05 * rng.standard_normal(g.size) for g in grids]
    return grids, values


def test_irregular_curves_are_fitted_independently() -> None:
    grids, values = irregular()
    basis = spline(10)
    result = smooth(values, grids, basis=basis, lam=1e-4)
    assert np.asarray(result.fd.coefs).shape == (10, 3)
    assert len(result.y2c_map) == 3
    total_df = 0.0
    for k, (g, v) in enumerate(zip(grids, values, strict=True)):
        single = smooth(v, g, basis=basis, lam=1e-4)
        np.testing.assert_allclose(np.asarray(result.fd.coefs)[:, k], single.fd.coefs[:, 0])
        assert result.gcv[k] == pytest.approx(float(single.gcv))
        total_df += single.df
    assert result.df == pytest.approx(total_df)


def test_irregular_weights_per_curve_or_shared() -> None:
    grids, values = irregular()
    grids = [grids[0], grids[0]]
    values = [values[0], 2.0 * values[0]]
    weights = np.linspace(1.0, 3.0, grids[0].size)
    shared = smooth(values, grids, basis=spline(10), lam=1e-3, weights=weights)
    listed = smooth(values, grids, basis=spline(10), lam=1e-3, weights=[weights, weights])
    np.testing.assert_allclose(shared.fd.coefs, listed.fd.coefs)


def test_irregular_single_curve_and_default_basis() -> None:
    grids, values = irregular()
    result = smooth(values[:1], grids[:1], lam="gcv")
    assert result.fd.n_curves == 1
    assert not isinstance(result.y2c_map, tuple)
    assert np.ndim(result.gcv) == 0
    assert result.fd.basis.domain == (float(grids[0].min()), float(grids[0].max()))


def test_irregular_input_validation() -> None:
    grids, values = irregular()
    with pytest.raises(ValueError, match="sequence matching"):
        smooth(np.zeros(3), grids, basis=spline(), lam=1.0)
    with pytest.raises(ValueError, match="sequence matching"):
        smooth(values[:2], grids, basis=spline(), lam=1.0)
    with pytest.raises(ValueError, match="as many points"):
        smooth([values[1], values[0], values[2]], grids, basis=spline(), lam=1.0)
    with pytest.raises(ValueError, match="shared argument vector"):
        smooth(values, grids, constraint="positive")


def test_empty_sequence_is_not_irregular() -> None:
    assert not sm._is_irregular([])
    assert not sm._is_irregular(np.zeros(3))
    assert sm._is_irregular([np.zeros(3)])


# --------------------------------------------------------------------------- #
# constrained fits
# --------------------------------------------------------------------------- #


def test_unknown_constraint_raises() -> None:
    with pytest.raises(ValueError, match="constraint must be one of"):
        smooth(noisy(), T, constraint="convex")


def test_positive_fit_is_positive_and_exact_on_an_exponential() -> None:
    y = np.exp(np.sin(2 * np.pi * T))
    result = smooth(y, T, basis=spline(15), lam=1e-10, constraint="positive")
    values = np.asarray(result(T))
    assert np.all(values > 0.0)
    np.testing.assert_allclose(values[:, 0], y, rtol=1e-3)
    assert result.constraint == "positive"
    assert result.y2c_map is None
    assert result.beta is None


def test_positive_derivative_is_the_chain_rule() -> None:
    y = np.exp(np.cos(3 * T))
    result = smooth(y, T, basis=spline(10), lam=1e-6, constraint="positive")
    h = 1e-6
    points = np.array([0.2, 0.5, 0.8])
    numeric = (np.asarray(result(points + h)) - np.asarray(result(points - h))) / (2 * h)
    np.testing.assert_allclose(result(points, 1), numeric, rtol=1e-5)


def test_monotone_fit_is_increasing() -> None:
    rng = np.random.default_rng(8)
    y = np.column_stack([T**2, np.tanh(4 * (T - 0.5))]) + 0.02 * rng.standard_normal((T.size, 2))
    result = smooth(y, T, basis=spline(8), lam=1e-4, penalty=2, constraint="monotone")
    fine = np.linspace(0.0, 1.0, 201)
    values = np.asarray(result(fine))
    assert values.shape == (201, 2)
    assert np.all(np.diff(values, axis=0) >= -1e-12)
    assert np.asarray(result.beta).shape == (2, 2)
    assert np.all(np.asarray(result(fine, 1)) > 0.0)
    assert np.asarray(result(fine, 2)).shape == (201, 2)
    with pytest.raises(ValueError, match="non-negative"):
        result(fine, -1)


def test_monotone_derivative_matches_finite_differences() -> None:
    y = np.log1p(5 * T)
    result = smooth(y, T, basis=spline(8), lam=1e-6, constraint="monotone")
    h = 1e-6
    points = np.array([0.25, 0.5, 0.75])
    numeric = (np.asarray(result(points + h)) - np.asarray(result(points - h))) / (2 * h)
    np.testing.assert_allclose(result(points, 1), numeric, rtol=1e-5)


def latent_result(latent: FData, constraint: str) -> SmoothResult:
    """A constrained result around a given latent ``W`` (monotone with ``β = (0, 1)``)."""
    n_curves = 1 if np.ndim(latent.coefs) == 1 else int(np.shape(latent.coefs)[1])
    beta = np.vstack([np.zeros(n_curves), np.ones(n_curves)])
    return SmoothResult(
        fd=latent,
        df=0.0,
        gcv=np.zeros(n_curves),
        sse=0.0,
        penalty_matrix=np.zeros((latent.basis.n_basis,) * 2),
        lam=0.0,
        y2c_map=None,
        beta=beta if constraint == "monotone" else None,
        constraint=constraint,
    )


def test_positive_derivatives_of_exp_t_squared_are_exact() -> None:
    # W = t², exp W and its derivatives in closed form (Hermite-type polynomials)
    latent = FData(np.array([[0.0], [0.0], [1.0]]), Monomial(domain=DOMAIN, n_basis=3))
    result = latent_result(latent, "positive")
    t = np.linspace(0.0, 1.0, 11)
    e = np.exp(t**2)
    closed = [
        e,
        2 * t * e,
        (2 + 4 * t**2) * e,
        (12 * t + 8 * t**3) * e,
        (12 + 48 * t**2 + 16 * t**4) * e,
    ]
    for order, expected in enumerate(closed):
        np.testing.assert_allclose(np.asarray(result(t, order))[:, 0], expected, rtol=1e-13)


def test_monotone_derivatives_of_a_linear_latent_are_exact() -> None:
    # W = 2t: x = (e^{2t} - 1) / 2 and Dⁿ x = 2^{n-1} e^{2t}
    latent = FData(np.array([[0.0], [2.0]]), Monomial(domain=DOMAIN, n_basis=2))
    result = latent_result(latent, "monotone")
    t = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose(
        np.asarray(result(t))[:, 0], (np.exp(2 * t) - 1) / 2, rtol=1e-13, atol=1e-15
    )
    for order in range(1, 7):
        np.testing.assert_allclose(
            np.asarray(result(t, order))[:, 0], 2.0 ** (order - 1) * np.exp(2 * t), rtol=1e-13
        )


@settings(max_examples=25, deadline=None)
@given(
    seed=st.integers(min_value=0, max_value=10_000),
    order=st.integers(min_value=1, max_value=4),
    constraint=st.sampled_from(["positive", "monotone"]),
)
def test_constrained_derivatives_match_finite_differences(
    seed: int, order: int, constraint: str
) -> None:
    rng = np.random.default_rng(seed)
    latent = FData(0.5 * rng.standard_normal((9, 2)), BSpline(domain=DOMAIN, n_basis=9, order=7))
    result = latent_result(latent, constraint)
    points = np.array([0.23, 0.51, 0.77])
    h = 1e-5
    numeric = (
        np.asarray(result(points + h, order - 1)) - np.asarray(result(points - h, order - 1))
    ) / (2 * h)
    exact = np.asarray(result(points, order))
    np.testing.assert_allclose(exact, numeric, rtol=1e-5, atol=1e-5 * np.max(np.abs(exact)))


def test_monotone_second_derivative_of_a_fit_matches_finite_differences() -> None:
    y = np.log1p(5 * T)
    result = smooth(y, T, basis=spline(8), lam=1e-6, constraint="monotone")
    h = 1e-5
    points = np.array([0.25, 0.5, 0.75])
    numeric = (np.asarray(result(points + h, 1)) - np.asarray(result(points - h, 1))) / (2 * h)
    np.testing.assert_allclose(result(points, 2), numeric, rtol=1e-6)


def test_constrained_derivatives_keep_torch_tensors() -> None:
    torch = pytest.importorskip("torch")
    latent = FData(
        torch.tensor([[0.1], [0.4], [-0.2], [0.3], [0.2]], dtype=torch.float64),
        BSpline(domain=DOMAIN, n_basis=5),
    )
    points = torch.linspace(0.0, 1.0, 5, dtype=torch.float64)
    for constraint in ("positive", "monotone"):
        result = latent_result(latent, constraint)
        if result.beta is not None:
            result = replace(result, beta=torch.from_numpy(result.beta))
        value = result(points, 3)
        assert isinstance(value, torch.Tensor)
        reference = latent_result(FData(latent.coefs.numpy(), latent.basis), constraint)
        np.testing.assert_allclose(value.numpy(), reference(points.numpy(), 3), rtol=1e-12)


def test_single_monotone_curve_evaluates_one_dimensional_points() -> None:
    result = smooth(T**3, T, basis=spline(8), lam=1e-6, constraint="monotone")
    assert np.asarray(result(0.5)).size == 1


def test_morph_maps_the_domain_onto_itself() -> None:
    y = T**2
    result = smooth(y, T, basis=spline(8), lam=1e-4, constraint="morph")
    values = np.asarray(result(np.array([0.0, 1.0])))[:, 0]
    np.testing.assert_allclose(values, [0.0, 1.0], atol=1e-12)
    np.testing.assert_allclose(np.asarray(result.beta)[0], 0.0)


def test_constrained_fit_with_gcv_lambda() -> None:
    y = np.exp(np.sin(2 * np.pi * T)) + 0.01 * np.cos(17 * T)
    result = smooth(y, T, basis=spline(10), constraint="positive")
    assert result.lam > 0.0
    assert result.df > 0.0
    assert np.isfinite(result.gcv)


def test_constrained_interpolation_has_infinite_gcv() -> None:
    t = np.linspace(0.0, 1.0, 6)
    result = smooth(
        np.exp(t), t, basis=BSpline(domain=DOMAIN, n_basis=6), lam=0.0, constraint="positive"
    )
    assert np.isinf(result.gcv)


def test_monotone_result_without_beta_raises() -> None:
    fit = smooth(T, T, basis=spline(6), lam=1e-4, constraint="monotone")
    broken = SmoothResult(
        fd=fit.fd,
        df=fit.df,
        gcv=fit.gcv,
        sse=fit.sse,
        penalty_matrix=fit.penalty_matrix,
        lam=fit.lam,
        y2c_map=None,
        beta=None,
        constraint="monotone",
    )
    with pytest.raises(ValueError, match="beta"):
        broken(T)


def test_gauss_newton_stops_when_no_step_improves() -> None:
    def evaluate(c: np.ndarray) -> tuple[float, None]:
        return float(c @ c), None

    def ascent(c: np.ndarray, parts: Any) -> np.ndarray:
        return np.ones_like(c)

    start = np.zeros(3)
    coefs, _ = sm._gauss_newton(start, evaluate, ascent, np)
    np.testing.assert_array_equal(coefs, start)


def test_morph_beta_with_a_degenerate_integral() -> None:
    beta = sm._monotone_beta(np.zeros(4), np.zeros(4), np.ones(4), (0.0, 2.0), np)
    np.testing.assert_allclose(beta, [0.0, 1.0])


@settings(max_examples=10, deadline=None)
@given(st.integers(min_value=0, max_value=2**31 - 1))
def test_monotone_fits_are_monotone_on_random_data(seed: int) -> None:
    rng = np.random.default_rng(seed)
    y = np.cumsum(np.abs(rng.standard_normal(T.size))) + 0.3 * rng.standard_normal(T.size)
    result = smooth(y, T, basis=spline(7), lam=1e-3, constraint="monotone")
    assert np.all(np.diff(np.asarray(result(T))[:, 0]) >= -1e-10)


# --------------------------------------------------------------------------- #
# Smoother
# --------------------------------------------------------------------------- #


def test_smoother_transform_returns_coefficients() -> None:
    curves = noisy(n_curves=5).T
    smoother = Smoother(spline(), t=T, lam=1e-4)
    coefs = smoother.fit_transform(curves)
    assert coefs.shape == (5, 12)
    np.testing.assert_allclose(coefs, np.asarray(smoother.fd_.coefs).T)
    np.testing.assert_allclose(smoother.transform(curves), coefs)
    assert smoother.n_features_in_ == T.size
    np.testing.assert_array_equal(smoother.t_, T)
    assert isinstance(smoother.result_, SmoothResult)


def test_smoother_reuses_the_fitted_lambda() -> None:
    curves = noisy(n_curves=4).T
    smoother = Smoother(spline(), t=T).fit(curves)
    new = noisy(n_curves=2, seed=9).T
    direct = smooth(new.T, T, basis=spline(), lam=smoother.result_.lam)
    np.testing.assert_allclose(smoother.transform(new), np.asarray(direct.fd.coefs).T)


def test_smoother_default_grid_and_fit_time_grid() -> None:
    curves = noisy(n_curves=3).T
    default = Smoother(lam=1e-2).fit(curves)
    np.testing.assert_array_equal(default.t_, np.arange(T.size, dtype=float))
    at_fit = Smoother(lam=1e-2).fit(curves, t=T)
    np.testing.assert_array_equal(at_fit.t_, T)


def test_smoother_grid_length_is_checked() -> None:
    with pytest.raises(ValueError, match="points"):
        Smoother(t=np.linspace(0.0, 1.0, 5)).fit(noisy().T)


def test_smoother_feature_names_come_from_the_basis() -> None:
    smoother = Smoother(spline(6), t=T, lam=1e-3).fit(noisy().T)
    names = smoother.get_feature_names_out()
    assert list(names) == list(spline(6).names)


def test_smoother_weights_and_penalty_are_used() -> None:
    weights = np.linspace(1.0, 2.0, T.size)
    curves = noisy(n_curves=2).T
    smoother = Smoother(spline(), t=T, lam=1e-3, penalty=3, weights=weights).fit(curves)
    direct = smooth(curves.T, T, basis=spline(), lam=1e-3, penalty=3, weights=weights)
    np.testing.assert_allclose(smoother.fd_.coefs, direct.fd.coefs)


def test_smoother_before_fit_raises() -> None:
    with pytest.raises(NotFittedError):
        Smoother().transform(noisy().T)


def test_smoother_tags_and_pipeline() -> None:
    tags = Smoother().__sklearn_tags__()
    assert tags.target_tags.required is False
    assert tags.input_tags.sparse is False
    pipe = Pipeline([("smooth", Smoother(spline(), t=T, lam=1e-3))])
    assert pipe.fit_transform(noisy().T).shape == (3, 12)


def test_smoother_get_params_round_trip() -> None:
    smoother = Smoother(spline(), t=T, lam=1e-2, penalty=3)
    assert Smoother(**smoother.get_params()).get_params()["penalty"] == 3


# --------------------------------------------------------------------------- #
# backends
# --------------------------------------------------------------------------- #


def test_torch_input_gives_torch_output_with_gradients() -> None:
    torch = pytest.importorskip("torch", reason="torch extra not installed")
    y = torch.tensor(noisy(n_curves=2), dtype=torch.float64, requires_grad=True)
    result = smooth(y, torch.tensor(T), basis=spline(), lam=1e-3)
    assert isinstance(result.fd.coefs, torch.Tensor)
    total = result.fd.coefs.sum()
    total.backward()
    assert y.grad is not None
    reference = smooth(noisy(n_curves=2), T, basis=spline(), lam=1e-3)
    np.testing.assert_allclose(result.fd.coefs.detach().numpy(), reference.fd.coefs, rtol=1e-10)
    # d(sum of coefficients)/dy is the column sum of the data-to-coefficient map.
    expected = np.asarray(reference.y2c_map).sum(axis=0)
    np.testing.assert_allclose(y.grad.numpy()[:, 0], expected, rtol=1e-10)


def test_fdata_result_type() -> None:
    assert isinstance(smooth(noisy(), T, basis=spline(), lam=1.0).fd, FData)
