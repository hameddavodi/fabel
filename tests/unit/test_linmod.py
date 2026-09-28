"""Unit tests for :func:`fabel.regression.linmod`, the bivariate-coefficient model."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fabel import LDO, BiFData, BSpline, FData, Fourier, Monomial, inprod
from fabel.regression import LinmodResult, linmod
from fabel.smoothing import smooth

SBASIS = BSpline(domain=(0.0, 1.0), n_basis=6)
TBASIS = BSpline(domain=(-1.0, 2.0), n_basis=5)
YBASIS = BSpline(domain=(-1.0, 2.0), n_basis=8)


def covariates(n: int, seed: int = 0, basis: BSpline = SBASIS) -> FData:
    return FData(np.random.default_rng(seed).standard_normal((basis.n_basis, n)), basis)


def exact_response(x: FData, surface: np.ndarray, intercept: float = 0.5) -> FData:
    """Curves ``intercept + ∫ x_i(s) β(s, t) ds`` held exactly in ``TBASIS``."""
    z = np.asarray(inprod(x, SBASIS))
    return FData(intercept + surface.T @ z.T, TBASIS)


def noisy_data(n: int = 20, seed: int = 1) -> tuple[FData, FData]:
    rng = np.random.default_rng(seed)
    return FData(rng.standard_normal((YBASIS.n_basis, n)), YBASIS), covariates(n, seed + 10)


def objective(model: LinmodResult, alpha: np.ndarray, surface: np.ndarray) -> float:
    """The penalised criterion of the fit, by brute-force quadrature."""
    t = np.linspace(-1.0, 2.0, 3001)
    s = np.linspace(0.0, 1.0, 1001)
    wt = np.full(t.size, t[1] - t[0])
    wt[[0, -1]] /= 2.0
    ws = np.full(s.size, s[1] - s[0])
    ws[[0, -1]] /= 2.0
    abasis, sbasis, tbasis = model.alpha.basis, model.beta.sbasis, model.beta.tbasis
    z = np.asarray(inprod(model.x, sbasis))
    fitted = (abasis(t) @ alpha)[:, None] + tbasis(t) @ surface.T @ z.T
    resid = model.y(t) - fitted
    fit = float(np.sum(np.asarray(model.weights) * (wt @ resid**2)))
    la, ls, lt = model.lam
    pa, ps, pt = model.penalty
    rough_alpha = float(wt @ (abasis(t, pa) @ alpha) ** 2)
    beta_s = sbasis(s, ps) @ surface @ tbasis(t).T
    beta_t = sbasis(s) @ surface @ tbasis(t, pt).T
    rough_s = float(ws @ beta_s**2 @ wt)
    rough_t = float(ws @ beta_t**2 @ wt)
    return fit + la * rough_alpha + ls * rough_s + lt * rough_t


# --------------------------------------------------------------------------- #
# the fit
# --------------------------------------------------------------------------- #


def test_recovers_an_exactly_determined_surface_on_different_continua() -> None:
    x = covariates(30)
    surface = np.random.default_rng(5).standard_normal((6, 5))
    y = exact_response(x, surface)
    model = linmod(y, x)
    assert isinstance(model.beta, BiFData)
    np.testing.assert_allclose(model.beta.coefs, surface, atol=1e-9)
    np.testing.assert_allclose(model.alpha(np.array([-1.0, 0.3, 2.0])), 0.5, atol=1e-10)
    np.testing.assert_allclose(model.fitted.coefs, y.coefs, atol=1e-9)
    np.testing.assert_allclose(model.residuals.coefs, 0.0, atol=1e-9)
    assert model.beta.domain == ((0.0, 1.0), (-1.0, 2.0))


def test_default_bases_follow_r() -> None:
    y, x = noisy_data()
    model = linmod(y, x, lam_s=1e-2, lam_t=1e-2)
    assert model.alpha.basis == YBASIS
    assert model.beta.sbasis == SBASIS
    assert model.beta.tbasis == YBASIS
    assert model.lam == (0.0, 1e-2, 1e-2)
    assert model.penalty == (LDO(2), LDO(2), LDO(2))
    size = YBASIS.n_basis + SBASIS.n_basis * YBASIS.n_basis
    assert model.cmat.shape == (size, size)
    np.testing.assert_allclose(model.cmat, model.cmat.T, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("lams", [(0.0, 0.0, 0.0), (1e-2, 1e-3, 0.0), (1.0, 1e-1, 1e-2)])
def test_solution_minimises_the_penalised_criterion(lams: tuple[float, float, float]) -> None:
    y, x = noisy_data(15)
    model = linmod(
        y, x, t_basis=TBASIS, alpha_basis=TBASIS, lam_alpha=lams[0], lam_s=lams[1], lam_t=lams[2]
    )
    alpha = np.asarray(model.alpha.coefs[:, 0])
    surface = np.asarray(model.beta.coefs)
    best = objective(model, alpha, surface)
    rng = np.random.default_rng(3)
    for _ in range(4):
        step = 1e-3
        worse_a = objective(model, alpha + step * rng.standard_normal(alpha.shape), surface)
        worse_b = objective(model, alpha, surface + step * rng.standard_normal(surface.shape))
        assert worse_a > best
        assert worse_b > best


def test_larger_lambdas_give_smoother_surfaces() -> None:
    y, x = noisy_data()
    op = LDO(2)
    rough = linmod(y, x, t_basis=TBASIS, lam_s=1e-6, lam_t=1e-6)
    smooth_s = linmod(y, x, t_basis=TBASIS, lam_s=1e2, lam_t=1e-6)
    smooth_t = linmod(y, x, t_basis=TBASIS, lam_s=1e-6, lam_t=1e2)

    def roughness(model: LinmodResult) -> tuple[float, float]:
        b = np.asarray(model.beta.coefs)
        gs, gt = SBASIS.gram(), TBASIS.gram()
        rs, rt = SBASIS.penalty(op), TBASIS.penalty(op)
        return float(np.trace(b.T @ rs @ b @ gt)), float(np.trace(b.T @ gs @ b @ rt))

    assert roughness(smooth_s)[0] < 1e-2 * roughness(rough)[0]
    assert roughness(smooth_t)[1] < 1e-2 * roughness(rough)[1]
    heavy = linmod(y, x, lam_alpha=1e6)
    assert float(np.sum(np.asarray(heavy.alpha.coefs) ** 2)) > 0.0
    second = heavy.alpha(np.linspace(-1.0, 2.0, 7), 2)
    np.testing.assert_allclose(second, 0.0, atol=1e-3)


def test_integer_weights_equal_replicated_curves() -> None:
    y, x = noisy_data(10)
    counts = np.array([1, 2, 3, 1, 1, 2, 1, 4, 1, 1])
    reps = np.repeat(np.arange(10), counts)
    weighted = linmod(y, x, lam_s=1e-2, lam_t=1e-3, lam_alpha=1e-2, weights=counts)
    replicated = linmod(
        FData(np.asarray(y.coefs)[:, reps], YBASIS),
        FData(np.asarray(x.coefs)[:, reps], SBASIS),
        lam_s=1e-2,
        lam_t=1e-3,
        lam_alpha=1e-2,
    )
    np.testing.assert_allclose(weighted.beta.coefs, replicated.beta.coefs, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(weighted.alpha.coefs, replicated.alpha.coefs, rtol=1e-9, atol=1e-12)
    unit = linmod(y, x, lam_s=1e-2, weights=np.ones(10))
    plain = linmod(y, x, lam_s=1e-2)
    np.testing.assert_array_equal(unit.beta.coefs, plain.beta.coefs)


@settings(max_examples=15, deadline=None)
@given(
    seed=st.integers(0, 10_000),
    scale=st.floats(-5.0, 5.0).filter(lambda v: abs(v) > 1e-3),
    lam=st.sampled_from([0.0, 1e-3, 1.0]),
)
def test_fit_is_linear_in_the_response(seed: int, scale: float, lam: float) -> None:
    rng = np.random.default_rng(seed)
    x = covariates(12, seed)
    y1 = FData(rng.standard_normal((5, 12)), TBASIS)
    y2 = FData(rng.standard_normal((5, 12)), TBASIS)
    both = FData(np.asarray(y1.coefs) + scale * np.asarray(y2.coefs), TBASIS)
    kw: dict[str, Any] = {"lam_alpha": lam, "lam_s": lam, "lam_t": lam}
    m1, m2, m12 = linmod(y1, x, **kw), linmod(y2, x, **kw), linmod(both, x, **kw)
    np.testing.assert_allclose(
        m12.beta.coefs, m1.beta.coefs + scale * m2.beta.coefs, rtol=1e-7, atol=1e-8
    )
    np.testing.assert_allclose(
        m12.fitted.coefs, m1.fitted.coefs + scale * m2.fitted.coefs, rtol=1e-7, atol=1e-8
    )


def test_fourier_and_monomial_bases_with_operator_penalties() -> None:
    period = 2.0
    fourier = Fourier(domain=(0.0, period), n_basis=7)
    rng = np.random.default_rng(8)
    x = FData(rng.standard_normal((7, 16)), fourier)
    y = FData(rng.standard_normal((7, 16)), fourier)
    harmonic = LDO.harmonic(period=period)
    model = linmod(
        y,
        x,
        alpha_basis=Fourier(domain=(0.0, period), n_basis=3),
        s_basis=Monomial(domain=(0.0, period), n_basis=3),
        lam_alpha=1e-3,
        lam_s=1e-3,
        lam_t=1e-3,
        penalty_alpha=harmonic,
        penalty_s=1,
        penalty_t=harmonic,
    )
    assert model.penalty == (harmonic, LDO(1), harmonic)
    assert model.beta.coefs.shape == (3, 7)
    assert np.all(np.isfinite(np.asarray(model.fitted.coefs)))


# --------------------------------------------------------------------------- #
# inputs and prediction
# --------------------------------------------------------------------------- #


def test_predict() -> None:
    y, x = noisy_data()
    model = linmod(y, x, t_basis=TBASIS, lam_s=1e-3, lam_t=1e-3)
    assert model.predict() is model.fitted
    np.testing.assert_allclose(model.predict(x).coefs, model.fitted.coefs, rtol=1e-12, atol=1e-13)
    # the same curves in a richer basis predict the same response
    richer = BSpline(domain=(0.0, 1.0), n_basis=12)
    t = np.linspace(0.0, 1.0, 400)
    refit = FData(np.linalg.lstsq(richer(t), x(t), rcond=None)[0], richer)
    np.testing.assert_allclose(model.predict(refit).coefs, model.fitted.coefs, atol=1e-8)
    new = model.predict(x[:3])
    assert new.n_curves == 3
    assert new.basis == YBASIS


def test_smooth_results_are_accepted() -> None:
    rng = np.random.default_rng(4)
    ts = np.linspace(0.0, 1.0, 30)
    tt = np.linspace(-1.0, 2.0, 40)
    xs = smooth(rng.standard_normal((30, 12)), ts, basis=SBASIS, lam=1e-4)
    ys = smooth(rng.standard_normal((40, 12)), tt, basis=TBASIS, lam=1e-4)
    from_smooth = linmod(ys, xs, lam_s=1e-2)
    from_fd = linmod(ys.fd, xs.fd, lam_s=1e-2)
    np.testing.assert_array_equal(from_smooth.beta.coefs, from_fd.beta.coefs)
    np.testing.assert_array_equal(from_smooth.predict(xs).coefs, from_fd.fitted.coefs)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"alpha_basis": SBASIS}, "alpha_basis lives on"),
        ({"t_basis": SBASIS}, "t_basis lives on"),
        ({"s_basis": TBASIS}, "s_basis lives on"),
        ({"lam_s": -1.0}, "lam_s must be finite"),
        ({"lam_t": float("inf")}, "lam_t must be finite"),
        ({"lam_alpha": float("nan")}, "lam_alpha must be finite"),
        ({"weights": np.r_[np.ones(19), 0.0]}, "weights must be positive"),
        ({"weights": np.ones(3)}, "must be a number or hold 20 values"),
    ],
)
def test_argument_validation(kwargs: dict[str, Any], match: str) -> None:
    y, x = noisy_data()
    with pytest.raises(ValueError, match=match):
        linmod(y, x, **kwargs)


def test_curve_validation() -> None:
    y, x = noisy_data()
    with pytest.raises(ValueError, match="x has 19 curves but y has 20"):
        linmod(y, x[:19], lam_s=1.0)
    with pytest.raises(TypeError, match="x must be an FData"):
        linmod(y, np.ones((20, 5)))
    with pytest.raises(TypeError, match="y must be an FData"):
        linmod([1.0, 2.0], x)
    multivariate = FData(np.ones((6, 20, 2)), SBASIS)
    with pytest.raises(ValueError, match="x must not carry a variable axis"):
        linmod(y, multivariate)
    model = linmod(y, x, lam_s=1.0, lam_t=1.0)
    with pytest.raises(ValueError, match="x lives on"):
        model.predict(FData(np.ones((5, 2)), TBASIS))
    with pytest.raises(TypeError, match="x must be an FData"):
        model.predict(np.ones(3))
