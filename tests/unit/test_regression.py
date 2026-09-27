"""Unit tests for :mod:`fabel.regression`."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fabel import LDO, BSpline, Constant, FData, Fourier, inprod
from fabel.regression import FRegress, FRegressCV, FRegressResult, FRegressStderr, fregress
from fabel.smoothing import smooth

RNG = np.random.default_rng(20260927)
BASIS = BSpline(domain=(0.0, 1.0), n_basis=7)


def curves(n: int, seed: int = 0, basis: BSpline = BASIS) -> FData:
    return FData(np.random.default_rng(seed).standard_normal((basis.n_basis, n)), basis)


def scalar_data(n: int = 25, seed: int = 1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    z1 = rng.standard_normal(n)
    z2 = rng.uniform(size=n)
    y = 1.0 + 2.0 * z1 - 0.5 * z2 + 0.2 * rng.standard_normal(n)
    return y, z1, z2


# --------------------------------------------------------------------------- #
# scalar response
# --------------------------------------------------------------------------- #


def test_scalar_on_scalar_is_ordinary_least_squares() -> None:
    y, z1, z2 = scalar_data()
    model = fregress(y, [1.0, z1, z2])
    design = np.column_stack([np.ones_like(z1), z1, z2])
    expected, *_ = np.linalg.lstsq(design, y, rcond=None)
    np.testing.assert_allclose(model.coefficients, expected, rtol=1e-12)
    np.testing.assert_allclose(model.fitted, design @ expected, rtol=1e-12)
    assert model.df == pytest.approx(3.0)
    sse = float(np.sum((y - design @ expected) ** 2))
    assert model.gcv == pytest.approx(sse / (len(y) - 3) ** 2)
    assert not model.functional_response
    assert model.names == ("x0", "x1", "x2")
    assert all(isinstance(b.basis, Constant) for b in model.beta)


def test_scalar_on_function_recovers_the_coefficient() -> None:
    x = curves(40, seed=2)
    truth = FData(np.linspace(-1.0, 1.0, 7), BASIS)
    y = 0.5 + np.asarray(inprod(x, truth))[:, 0]
    model = fregress(y, {"const": 1.0, "x": x})
    np.testing.assert_allclose(model.beta[1].coefs[:, 0], truth.coefs[:, 0], atol=1e-10)
    assert model.beta[0].coefs[0, 0] == pytest.approx(0.5)
    assert model.beta[1].basis == BASIS


def test_scalar_penalty_shrinks_roughness() -> None:
    x = curves(40, seed=3)
    y = np.random.default_rng(4).standard_normal(40)
    rough = fregress(y, [1.0, x])
    smoothed = fregress(y, [1.0, x], beta=[None, (BASIS, 1.0)])
    pen = BASIS.penalty(2)
    b_rough = rough.beta[1].coefs[:, 0]
    b_smooth = smoothed.beta[1].coefs[:, 0]
    assert b_smooth @ pen @ b_smooth < b_rough @ pen @ b_rough
    assert smoothed.df is not None
    assert rough.df is not None
    assert smoothed.df < rough.df


def test_ocv_is_the_leave_one_out_sum() -> None:
    y, z1, z2 = scalar_data(n=15)
    model = fregress(y, [1.0, z1, z2])
    loo = []
    for i in range(15):
        keep = np.arange(15) != i
        sub = fregress(y[keep], [1.0, z1[keep], z2[keep]])
        loo.append(y[i] - sub.predict([1.0, z1[i : i + 1], z2[i : i + 1]])[0])
    cv = model.cv()
    assert isinstance(cv, FRegressCV)
    np.testing.assert_allclose(cv.errors, loo, rtol=1e-10)
    assert cv.sse == pytest.approx(model.ocv, rel=1e-10)


def test_scalar_stderr_default_is_sigma2_times_inverse() -> None:
    y, z1, z2 = scalar_data()
    model = fregress(y, [1.0, z1, z2])
    se = model.stderr()
    assert isinstance(se, FRegressStderr)
    design = np.column_stack([np.ones_like(z1), z1, z2])
    sigma2 = float(np.sum((y - model.fitted) ** 2)) / (len(y) - 3)
    expected = sigma2 * np.linalg.inv(design.T @ design)
    np.testing.assert_allclose(se.cov, expected, rtol=1e-10)
    np.testing.assert_allclose(
        [b.coefs[0, 0] for b in se.beta], np.sqrt(np.diag(expected)), rtol=1e-10
    )
    np.testing.assert_allclose(
        se.fitted, np.sqrt(np.einsum("ij,jk,ik->i", design, expected, design)), rtol=1e-10
    )


def test_scalar_stderr_with_explicit_covariance_and_map() -> None:
    y, z1, _ = scalar_data(n=10)
    model = fregress(y, [1.0, z1])
    sigma = np.diag(np.linspace(0.5, 1.5, 10))
    se = model.stderr(sigma_e=sigma, y2c_map=np.eye(10))
    design = np.column_stack([np.ones(10), z1])
    bread = np.linalg.inv(design.T @ design) @ design.T
    np.testing.assert_allclose(se.cov, bread @ sigma @ bread.T, rtol=1e-10)
    with pytest.raises(ValueError, match="y2c_map"):
        model.stderr(y2c_map=np.eye(3))
    with pytest.raises(ValueError, match="sigma_e"):
        model.stderr(sigma_e=np.eye(3))


def test_scalar_stderr_of_functional_coefficient_is_projected() -> None:
    x = curves(30, seed=5)
    y = np.random.default_rng(6).standard_normal(30)
    model = fregress(y, [1.0, x], beta=[None, (BASIS, 1e-4)])
    se = model.stderr(sigma_e=1.0)
    # The stored curve is the L2 projection of sqrt(theta(t)' V theta(t)):
    # check it against an independent dense trapezoid projection.
    t = np.linspace(0.0, 1.0, 40001)
    theta = BASIS(t)
    pointwise = np.sqrt(np.einsum("ij,jk,ik->i", theta, se.cov[1:, 1:], theta))
    rhs = np.trapezoid(theta * pointwise[:, None], t, axis=0)
    expected = np.linalg.solve(BASIS.gram(), rhs)
    np.testing.assert_allclose(se.beta[1].coefs[:, 0], expected, rtol=1e-6)


def test_weights_match_replicated_observations() -> None:
    y, z1, z2 = scalar_data(n=12)
    counts = np.array([1, 2, 3, 1, 1, 2, 1, 1, 4, 1, 1, 2], dtype=float)
    weighted = fregress(y, [1.0, z1, z2], weights=counts)
    reps = np.repeat(np.arange(12), counts.astype(int))
    replicated = fregress(y[reps], [1.0, z1[reps], z2[reps]])
    np.testing.assert_allclose(weighted.coefficients, replicated.coefficients, rtol=1e-10)


@settings(max_examples=25, deadline=None)
@given(
    scale=st.floats(min_value=-5.0, max_value=5.0).filter(lambda v: abs(v) > 1e-3),
    seed=st.integers(min_value=0, max_value=10_000),
)
def test_fit_is_linear_in_the_response(scale: float, seed: int) -> None:
    rng = np.random.default_rng(seed)
    x = curves(20, seed=seed)
    y = rng.standard_normal(20)
    base = fregress(y, [1.0, x], beta=[None, (BASIS, 1e-3)])
    scaled = fregress(scale * y, [1.0, x], beta=[None, (BASIS, 1e-3)])
    np.testing.assert_allclose(
        scaled.coefficients, scale * base.coefficients, rtol=1e-8, atol=1e-10
    )
    assert scaled.df == pytest.approx(base.df)


# --------------------------------------------------------------------------- #
# functional response
# --------------------------------------------------------------------------- #


def test_function_on_scalar_recovers_exact_coefficients() -> None:
    beta0 = curves(1, seed=7)
    beta1 = curves(1, seed=8)
    z = RNG.standard_normal(15)
    y = FData(beta0.coefs + beta1.coefs * z[None, :], BASIS)
    model = fregress(y, [1.0, z])
    np.testing.assert_allclose(model.beta[0].coefs, beta0.coefs, atol=1e-10)
    np.testing.assert_allclose(model.beta[1].coefs, beta1.coefs, atol=1e-10)
    np.testing.assert_allclose(model.fitted.coefs, y.coefs, atol=1e-10)
    assert model.functional_response
    assert model.df is None
    assert model.gcv is None
    assert model.ocv is None


def test_function_on_scalar_normal_equations_match_kronecker_form() -> None:
    y = curves(10, seed=9)
    z = RNG.standard_normal(10)
    model = fregress(y, [1.0, z], lam=1e-3)
    zmat = np.column_stack([np.ones(10), z])
    gram = BASIS.gram()
    cmat = np.kron(zmat.T @ zmat, gram) + np.kron(np.eye(2), 1e-3 * BASIS.penalty(2))
    dmat = np.concatenate([gram @ y.coefs @ zmat[:, j] for j in range(2)])
    np.testing.assert_allclose(model.cmat, cmat, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(model.dmat, dmat, rtol=1e-10, atol=1e-12)


def test_concurrent_model_recovers_exact_coefficient() -> None:
    xbasis = BSpline(domain=(0.0, 1.0), n_basis=5, order=3)
    bbasis = BSpline(domain=(0.0, 1.0), n_basis=5, order=3)
    ybasis = xbasis * bbasis
    x = curves(12, seed=10, basis=xbasis)
    beta = curves(1, seed=11, basis=bbasis)
    nodes = np.linspace(0.0, 1.0, 200)
    values = x(nodes) * beta(nodes)
    fitted = smooth(values, nodes, basis=ybasis, lam=0.0).fd
    model = fregress(fitted, [x], beta=[bbasis])
    np.testing.assert_allclose(model.beta[0].coefs, beta.coefs, atol=1e-8)
    grid = np.linspace(0.0, 1.0, 17)
    np.testing.assert_allclose(model.fitted(grid), fitted(grid), atol=1e-8)


def test_mixed_functional_and_scalar_covariates_predict() -> None:
    y = curves(14, seed=12)
    x = curves(14, seed=13)
    z = RNG.standard_normal(14)
    model = fregress(y, {"const": 1.0, "z": z, "x": x}, lam=1e-4)
    again = model.predict({"const": 1.0, "z": z, "x": x})
    np.testing.assert_allclose(again.coefs, model.fitted.coefs, atol=1e-10)
    new = model.predict([1.0, z[:3], x[:3]])
    np.testing.assert_allclose(new.coefs, model.fitted.coefs[:, :3], atol=1e-10)
    assert model.predict() is model.fitted


def test_functional_cv_matches_brute_force_refit() -> None:
    y = curves(8, seed=14)
    z = RNG.standard_normal(8)
    model = fregress(y, [1.0, z], lam=1e-3)
    cv = model.cv()
    for i in (0, 5):
        keep = np.arange(8) != i
        sub = fregress(y[keep], [1.0, z[keep]], lam=1e-3)
        pred = sub.predict([1.0, z[i : i + 1]])
        np.testing.assert_allclose(
            cv.errors.coefs[:, i], y.coefs[:, i] - pred.coefs[:, 0], atol=1e-9
        )
    gram = BASIS.gram()
    assert cv.sse == pytest.approx(float(np.sum(cv.errors.coefs * (gram @ cv.errors.coefs))))


def test_functional_cv_with_functional_covariate() -> None:
    y = curves(6, seed=15)
    x = curves(6, seed=16)
    model = fregress(y, [1.0, x], lam=1e-2)
    cv = model.cv()
    keep = np.arange(6) != 2
    sub = fregress(y[keep], [1.0, x[keep]], lam=1e-2)
    pred = sub.predict([1.0, x[2:3]])
    np.testing.assert_allclose(cv.errors.coefs[:, 2], y.coefs[:, 2] - pred.coefs[:, 0], atol=1e-9)


def _observed_design() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.linspace(0.0, 1.0, 31)
    z = np.random.default_rng(17).standard_normal(9)
    obs = np.random.default_rng(18).standard_normal((31, 9))
    return t, z, obs


def test_functional_stderr_is_the_exact_sandwich() -> None:
    t, z, obs = _observed_design()
    sm = smooth(obs, t, basis=BASIS, lam=1e-4)
    model = fregress(sm, [1.0, z], lam=1e-3, penalty=LDO(2))
    sigma = np.diag(np.linspace(0.2, 1.0, 31))
    se = model.stderr(sigma_e=sigma)
    # b is linear in the observations: build its Jacobian column by column.
    jac = np.zeros((model.cmat.shape[0], obs.size))
    for k in range(obs.size):
        unit = np.zeros(obs.size)
        unit[k] = 1.0
        smk = smooth(unit.reshape(obs.shape, order="F"), t, basis=BASIS, lam=1e-4)
        jac[:, k] = fregress(smk, [1.0, z], lam=1e-3).coefficients
    full = np.kron(np.eye(9), sigma)
    np.testing.assert_allclose(se.cov, jac @ full @ jac.T, rtol=1e-8, atol=1e-14)
    assert len(se.beta) == 2
    assert se.fitted.n_curves == 9


def test_functional_stderr_explicit_map_and_scalar_sigma() -> None:
    t, z, obs = _observed_design()
    sm = smooth(obs, t, basis=BASIS, lam=1e-4)
    x = curves(9, seed=19)
    model = fregress(sm.fd, [1.0, z, x], lam=1e-3)
    se = model.stderr(sigma_e=0.3, y2c_map=sm.y2c_map)
    assert se.cov.shape == (21, 21)
    assert np.all(np.linalg.eigvalsh(se.cov) > -1e-12)
    assert np.all(se.fitted.coefs > -1e-8)


def test_functional_stderr_needs_map_and_covariance() -> None:
    t, z, obs = _observed_design()
    sm = smooth(obs, t, basis=BASIS, lam=1e-4)
    with pytest.raises(ValueError, match="y2c_map"):
        fregress(sm.fd, [1.0, z]).stderr(sigma_e=1.0)
    with pytest.raises(ValueError, match="sigma_e"):
        fregress(sm, [1.0, z]).stderr()
    with pytest.raises(ValueError, match="rows"):
        fregress(sm, [1.0, z]).stderr(sigma_e=1.0, y2c_map=np.eye(3))


def test_fourier_response_with_harmonic_penalty() -> None:
    basis = Fourier(domain=(0.0, 365.0), n_basis=9)
    y = FData(RNG.standard_normal((9, 10)), basis)
    region = ["a", "b"] * 5
    model = fregress(
        "y ~ region", {"y": y, "region": region}, lam=10.0, penalty=LDO.harmonic(365.0)
    )
    assert model.names == ("const", "region.b")
    assert model.beta[0].basis == basis


# --------------------------------------------------------------------------- #
# formula interface
# --------------------------------------------------------------------------- #


def test_formula_expands_categoricals_like_r() -> None:
    y, z1, _ = scalar_data(n=6)
    data = {"y": y, "g": np.array(["b", "a", "c", "a", "b", "c"], dtype=object), "z": z1}
    with_const = fregress("y ~ g + z", data)
    assert with_const.names == ("const", "g.b", "g.c", "z")
    no_const = fregress("y ~ g - 1", data)
    assert no_const.names == ("g.a", "g.b", "g.c")
    zero = fregress("y ~ 0 + z", data)
    assert zero.names == ("z",)
    dup = fregress("y ~ z + z", data)
    assert dup.names == ("const", "z")


def test_formula_beta_by_source_variable() -> None:
    y = curves(6, seed=20)
    other = BSpline(domain=(0.0, 1.0), n_basis=5)
    data = {"y": y, "g": ["p", "q", "r", "p", "q", "r"]}
    model = fregress("y ~ g", data, beta={"g": other, "const": (BASIS, 1e-3)})
    assert [b.basis for b in model.beta] == [BASIS, other, other]
    assert model.terms[0].lam == 1e-3


@pytest.mark.parametrize(
    ("formula", "error"),
    [
        ("no tilde", ValueError),
        ("y ~ z - w", ValueError),
        ("y ~ z*w", ValueError),
        ("y ~ 0", ValueError),
        ("y ~ missing", KeyError),
    ],
)
def test_formula_errors(formula: str, error: type[Exception]) -> None:
    y, z1, _ = scalar_data(n=5)
    with pytest.raises(error):
        fregress(formula, {"y": y, "z": z1, "w": z1})


def test_formula_needs_mapping_and_rejects_duplicate_names() -> None:
    y, z1, _ = scalar_data(n=4)
    with pytest.raises(TypeError, match="mapping"):
        fregress("y ~ z", [z1])
    data = {"y": y, "g": ["a", "b", "a", "b"], "g.b": z1}
    with pytest.raises(ValueError, match="unique"):
        fregress("y ~ g + g.b", data)


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #


def test_beta_specification_forms() -> None:
    x = curves(12, seed=21)
    y = RNG.standard_normal(12)
    model = fregress(y, [1.0, x], beta=[None, (BASIS,)], lam=0.5, penalty=1)
    assert model.terms[1].lam == 0.5
    assert model.terms[1].penalty == LDO(1)
    model = fregress(y, [1.0, x], beta=[None, (BASIS, 0.1, LDO(3))])
    assert model.terms[1].penalty == LDO(3)
    model = fregress(y, [1.0, x], beta={"x1": BASIS})
    assert model.terms[1].basis == BASIS
    with pytest.raises(TypeError, match="beta entry"):
        fregress(y, [1.0, x], beta=[None, "spline"])
    with pytest.raises(ValueError, match="entries"):
        fregress(y, [1.0, x], beta=[None])
    with pytest.raises(ValueError, match="Constant"):
        fregress(y, [x, 1.0], beta=[None, BASIS])
    with pytest.raises(ValueError, match="non-negative"):
        fregress(y, [1.0, x], beta=[None, (BASIS, -1.0)])
    with pytest.raises(ValueError, match="lives on"):
        fregress(y, [1.0, x], beta=[None, BSpline(domain=(0.0, 2.0), n_basis=5)])


def test_single_basis_applies_to_every_functional_response_term() -> None:
    y = curves(8, seed=22)
    other = BSpline(domain=(0.0, 1.0), n_basis=4)
    model = fregress(y, [1.0, RNG.standard_normal(8)], beta=other)
    assert all(b.basis == other for b in model.beta)


@pytest.mark.parametrize(
    ("x", "match"),
    [
        ([np.ones(3)], "hold 10 values"),
        ([curves(4)], "curves, expected 10"),
        ([FData(np.ones((7, 10, 2)), BASIS)], "variable axis"),
        ([curves(10, basis=BSpline(domain=(0.0, 2.0), n_basis=7))], "lives on"),
    ],
)
def test_covariate_validation(x: list[Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        fregress(curves(10, seed=23), x)


def test_response_validation() -> None:
    x = RNG.standard_normal(6)
    with pytest.raises(ValueError, match="vector"):
        fregress(np.ones((6, 2)), [x])
    with pytest.raises(ValueError, match="variable axis"):
        fregress(FData(np.ones((7, 6, 2)), BASIS), [x])
    with pytest.raises(ValueError, match="covariates x"):
        fregress(np.ones(6))
    with pytest.raises(ValueError, match="at least one"):
        fregress(np.ones(6), [])
    with pytest.raises(ValueError, match="positive"):
        fregress(np.ones(6), [x], weights=-np.ones(6))
    column = fregress(np.arange(6.0)[:, None], [1.0, x])
    assert column.fitted.shape == (6,)
    single = fregress(np.arange(6.0), x)
    assert single.names == ("x0",)


def test_constrained_smooth_is_rejected() -> None:
    t = np.linspace(0.0, 1.0, 20)
    positive = smooth(np.exp(np.sin(3 * t))[:, None].repeat(3, axis=1), t, constraint="positive")
    with pytest.raises(ValueError, match="constrained"):
        fregress(positive, [1.0])


def test_predict_validation() -> None:
    y = curves(10, seed=24)
    z = RNG.standard_normal(10)
    x = curves(10, seed=25)
    model = fregress(y, {"c": 1.0, "z": z, "x": x})
    with pytest.raises(ValueError, match="missing"):
        model.predict({"c": 1.0})
    with pytest.raises(ValueError, match="expected 3"):
        model.predict([1.0, z])
    with pytest.raises(ValueError, match="must be a curve"):
        model.predict([1.0, z, z])
    with pytest.raises(ValueError, match="number per observation"):
        model.predict([1.0, x, x])
    with pytest.raises(ValueError, match="disagree"):
        model.predict([1.0, z[:4], x[:3]])


def test_scalar_predict_with_numbers_only() -> None:
    y, z1, _ = scalar_data(n=8)
    model = fregress(y, [1.0, z1])
    np.testing.assert_allclose(
        model.predict([1.0, 2.0]), [model.coefficients[0] + 2.0 * model.coefficients[1]]
    )


# --------------------------------------------------------------------------- #
# scikit-learn estimator
# --------------------------------------------------------------------------- #


def test_estimator_on_matrix_matches_function() -> None:
    y, z1, z2 = scalar_data()
    matrix = np.column_stack([z1, z2])
    est = FRegress().fit(matrix, y)
    ref = fregress(y, [1.0, z1, z2])
    np.testing.assert_allclose(est.predict(matrix), ref.fitted, rtol=1e-12)
    assert isinstance(est.result_, FRegressResult)
    assert len(est.coef_) == 3
    assert est.n_features_in_ == 2
    assert est.score(matrix, y) > 0.9
    no_int = FRegress(fit_intercept=False).fit(matrix, y)
    assert len(no_int.coef_) == 2


def test_estimator_with_functional_covariates() -> None:
    x = curves(30, seed=26)
    truth = FData(np.linspace(1.0, -1.0, 7), BASIS)
    y = np.asarray(inprod(x, truth))[:, 0]
    est = FRegress().fit(x, y)
    np.testing.assert_allclose(est.predict(x), y, atol=1e-8)
    z = RNG.standard_normal(30)
    mixed = FRegress(beta=[None, None, (BASIS, 1e-6)]).fit([z, x], y)
    assert mixed.predict([z, x]).shape == (30,)


def test_estimator_functional_response_and_errors() -> None:
    y = curves(12, seed=27)
    matrix = RNG.standard_normal((12, 1))
    est = FRegress(lam=1e-3).fit(matrix, y)
    assert est.predict(matrix).n_curves == 12
    with pytest.raises(ValueError, match="samples"):
        FRegress().fit(np.ones((5, 1)), np.ones(6))
