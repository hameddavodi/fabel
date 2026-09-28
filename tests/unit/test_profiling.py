"""Unit tests for ``fdatools.profiling`` (generalized profiling of ODE parameters)."""

from __future__ import annotations

import importlib.util
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import fdatools as fdt
from fdatools.profiling import (
    CSTR_CONDITIONS,
    CSTR_PARAMETERS,
    InnerFit,
    ODEModel,
    ProfiledODE,
    ProfileResult,
    cstr_inputs,
    cstr_model,
    fitzhugh_nagumo_model,
    profile_ode,
    simpson_rule,
)

HAS_TORCH = importlib.util.find_spec("torch") is not None
FHN_TRUE = np.array([0.2, 0.2, 3.0])


def decay_model(**kw: Any) -> ODEModel:
    return ODEModel(lambda x, t, th: -th[0] * x, n_states=1, n_params=1, **kw)


@pytest.fixture(scope="module")
def fhn_data() -> tuple[np.ndarray, np.ndarray]:
    fhn = fitzhugh_nagumo_model()
    t = np.linspace(0.0, 20.0, 201)
    return t, fhn.simulate(t, [-1.0, 1.0], FHN_TRUE)


@pytest.fixture(scope="module")
def fhn_basis() -> fdt.BSpline:
    return fdt.BSpline(domain=(0.0, 20.0), breaks=np.linspace(0.0, 20.0, 201).tolist())


# --------------------------------------------------------------------- model


def test_model_validation() -> None:
    with pytest.raises(ValueError, match="positive"):
        ODEModel(lambda x, t, th: x, n_states=0, n_params=1)
    with pytest.raises(ValueError, match="state names"):
        ODEModel(lambda x, t, th: x, n_states=2, n_params=1, state_names=("a",))
    m = decay_model()
    assert m.state_names == ("x0",)
    assert m.param_names == ("theta0",)


def test_model_argument_checks() -> None:
    m = decay_model()
    with pytest.raises(ValueError, match="x must have shape"):
        m(np.ones(3), np.zeros(3), [1.0])
    with pytest.raises(ValueError, match="t must have shape"):
        m(np.ones((3, 1)), np.zeros(2), [1.0])
    with pytest.raises(ValueError, match="theta must have shape"):
        m(np.ones((3, 1)), np.zeros(3), [1.0, 2.0])
    bad = ODEModel(lambda x, t, th: x[:, 0], n_states=1, n_params=1)
    with pytest.raises(ValueError, match="rhs returned shape"):
        bad(np.ones((3, 1)), np.zeros(3), [1.0])


def _points(seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    return rng.uniform(-2, 2, (7, 2)), rng.uniform(0, 5, 7)


def test_finite_difference_derivatives_match_analytic() -> None:
    exact = fitzhugh_nagumo_model()
    numeric = ODEModel(exact.rhs, n_states=2, n_params=3)
    x, t = _points()
    for got, want in zip(
        numeric.jacobians(x, t, FHN_TRUE), exact.jacobians(x, t, FHN_TRUE), strict=True
    ):
        np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-8)
    for got, want in zip(
        numeric.hessians(x, t, FHN_TRUE), exact.hessians(x, t, FHN_TRUE), strict=True
    ):
        np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-5)


def test_partial_derivatives_are_mixed() -> None:
    exact = fitzhugh_nagumo_model()
    x, t = _points(1)
    only_x = ODEModel(exact.rhs, 2, 3, jac_x=exact.jac_x)
    only_theta = ODEModel(exact.rhs, 2, 3, jac_theta=exact.jac_theta)
    for m in (only_x, only_theta):
        for got, want in zip(
            m.jacobians(x, t, FHN_TRUE), exact.jacobians(x, t, FHN_TRUE), strict=True
        ):
            np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-8)


def test_bad_jacobian_shape_is_reported() -> None:
    m = ODEModel(lambda x, t, th: x, 1, 1, jac_x=lambda x, t, th: np.ones((len(t), 2)))
    with pytest.raises(ValueError, match="jac_x returned shape"):
        m.jacobians(np.ones((2, 1)), np.zeros(2), [1.0])


def test_simulate_decay_and_breaks() -> None:
    m = decay_model()
    t = np.linspace(0.0, 2.0, 9)
    x = m.simulate(t, [2.0], [1.5])
    np.testing.assert_allclose(x[:, 0], 2.0 * np.exp(-1.5 * t), rtol=1e-8)
    xb = m.simulate(t, [2.0], [1.5], breaks=[0.5, 1.0, 1.1, 5.0])
    np.testing.assert_allclose(xb, x, rtol=1e-8)
    single = m.simulate([0.0], [2.0], [1.5])
    assert single.tolist() == [[2.0]]


def test_simulate_errors() -> None:
    m = decay_model()
    with pytest.raises(ValueError, match="strictly increasing"):
        m.simulate([1.0, 0.0], [1.0], [1.0])
    with pytest.raises(ValueError, match="x0 must have shape"):
        m.simulate([0.0, 1.0], [1.0, 2.0], [1.0])
    blowup = ODEModel(lambda x, t, th: x**2, 1, 1)
    with pytest.raises(ValueError, match="solver failed"):
        blowup.simulate([0.0, 2.0], [1.0], [0.0], method="RK45")


@pytest.mark.skipif(not HAS_TORCH, reason="torch not installed")
def test_from_torch_matches_analytic() -> None:
    import torch

    def rhs(x: Any, t: Any, th: Any) -> Any:
        v, r = x[:, 0], x[:, 1]
        a, b, c = th[0], th[1], th[2]
        return torch.stack([c * (v - v**3 / 3 + r), -(v - a + b * r) / c], dim=1)

    auto = ODEModel.from_torch(rhs, 2, 3, state_names=("V", "R"), param_names=("a", "b", "c"))
    exact = fitzhugh_nagumo_model()
    x, t = _points(2)
    np.testing.assert_allclose(auto(x, t, FHN_TRUE), exact(x, t, FHN_TRUE), rtol=1e-12)
    for got, want in zip(
        auto.jacobians(x, t, FHN_TRUE), exact.jacobians(x, t, FHN_TRUE), strict=True
    ):
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
    for got, want in zip(
        auto.hessians(x, t, FHN_TRUE), exact.hessians(x, t, FHN_TRUE), strict=True
    ):
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
    assert auto.state_names == ("V", "R")


def test_import_fdatools_profiling_does_not_import_torch() -> None:
    import subprocess
    import sys

    code = "import sys, fdatools.profiling; print('torch' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


# ---------------------------------------------------------------------- CSTR


@pytest.mark.parametrize("condition", CSTR_CONDITIONS)
def test_cstr_inputs_shape_and_nominal(condition: str) -> None:
    u = cstr_inputs([70.0, -3.0], condition)
    assert u.shape == (2, 5)
    np.testing.assert_array_equal(u[0], u[1])


def test_cstr_inputs_unknown() -> None:
    with pytest.raises(ValueError, match="unknown CSTR condition"):
        cstr_inputs([0.0], "all.cool.ramp")


def test_cstr_model_validation() -> None:
    with pytest.raises(ValueError, match="unknown CSTR value"):
        cstr_model(Vol=2.0)
    with pytest.raises(ValueError, match="estimate must be"):
        cstr_model(estimate=("kref", "kref"))
    with pytest.raises(ValueError, match="estimate must be"):
        cstr_model(estimate=())
    with pytest.raises(ValueError, match="unknown CSTR condition"):
        cstr_model("nope")


@pytest.mark.parametrize("estimate", [tuple(CSTR_PARAMETERS), ("EoverR",), ("b", "kref", "a")])
def test_cstr_derivatives_match_finite_differences(estimate: tuple[str, ...]) -> None:
    exact = cstr_model("all.hot.step", estimate=estimate, V=1.3, delH=-120.0)
    numeric = ODEModel(exact.rhs, 2, len(estimate))
    rng = np.random.default_rng(3)
    x = np.column_stack([rng.uniform(0.5, 2, 9), rng.uniform(330, 370, 9)])
    t = rng.uniform(0, 64, 9)
    theta = np.array([CSTR_PARAMETERS[n] for n in estimate]) * 1.1
    for got, want in zip(numeric.jacobians(x, t, theta), exact.jacobians(x, t, theta), strict=True):
        np.testing.assert_allclose(got, want, rtol=1e-7, atol=1e-7 * np.abs(want).max())
    # second derivatives by differencing Jacobians whose entries span ~1e4 in
    # scale (C ~ 1, T ~ 350): the difference noise is ~1e-3 of the largest entry
    for got, want in zip(numeric.hessians(x, t, theta), exact.hessians(x, t, theta), strict=True):
        np.testing.assert_allclose(got, want, rtol=1e-4, atol=1e-3 * np.abs(want).max())
    assert exact.param_names == estimate


def test_cstr_custom_inputs() -> None:
    def constant(t: np.ndarray) -> np.ndarray:
        return np.tile([1.0, 2.0, 323.0, 335.0, 15.0], (len(t), 1))

    custom = cstr_model(constant)
    builtin = cstr_model("all.cool.step")
    x = np.array([[1.5, 341.0]])
    theta = np.array(list(CSTR_PARAMETERS.values()))
    np.testing.assert_allclose(custom(x, [1.0], theta), builtin(x, [1.0], theta), rtol=1e-14)
    wrong = cstr_model(lambda t: np.ones((len(t), 4)))
    with pytest.raises(ValueError, match="CSTR inputs returned shape"):
        wrong(x, [1.0], theta)


# ---------------------------------------------------------------- quadrature


def test_simpson_errors() -> None:
    with pytest.raises(ValueError, match="odd integer"):
        simpson_rule([0.0, 1.0], 4)
    with pytest.raises(ValueError, match="strictly increasing"):
        simpson_rule([1.0, 0.0])


@settings(max_examples=40, deadline=None)
@given(
    edges=st.lists(st.floats(-5, 5), min_size=2, max_size=6, unique=True),
    half=st.integers(1, 4),
    coefs=st.lists(st.floats(-3, 3), min_size=4, max_size=4),
)
def test_simpson_is_exact_for_cubics(edges: list[float], half: int, coefs: list[float]) -> None:
    breaks = sorted(edges)
    if min(np.diff(breaks)) < 1e-3:
        return
    nodes, weights = simpson_rule(breaks, 2 * half + 1)
    poly = np.polynomial.Polynomial(coefs)
    exact = poly.integ()(breaks[-1]) - poly.integ()(breaks[0])
    assert float(weights @ poly(nodes)) == pytest.approx(exact, rel=1e-10, abs=1e-10)
    assert float(weights.sum()) == pytest.approx(breaks[-1] - breaks[0], rel=1e-12)


# ------------------------------------------------------------------- problem


def test_problem_input_forms_agree(fhn_data: Any) -> None:
    t, x = fhn_data
    basis = fdt.BSpline(domain=(0.0, 20.0), n_basis=30)
    y = x.copy()
    y[::3, 1] = np.nan
    a = ProfiledODE(fitzhugh_nagumo_model(), t, y, basis, lam=10.0)
    keep = ~np.isnan(y[:, 1])
    b = ProfiledODE(
        fitzhugh_nagumo_model(), [t, t[keep]], [y[:, 0], y[keep, 1]], [basis, basis], 10.0
    )
    assert a.n_obs == b.n_obs == 201 + int(keep.sum())
    c = np.linspace(-1, 1, a.n_coefs)
    assert a.criterion(c, FHN_TRUE) == pytest.approx(b.criterion(c, FHN_TRUE), rel=1e-14)
    unobserved = ProfiledODE(fitzhugh_nagumo_model(), t, [x[:, 0], None], basis, 10.0)
    assert unobserved.n_obs == 201


def test_problem_validation(fhn_data: Any) -> None:
    t, x = fhn_data
    fhn = fitzhugh_nagumo_model()
    b = fdt.BSpline(domain=(0.0, 20.0), n_basis=10)
    with pytest.raises(ValueError, match="one entry per state"):
        ProfiledODE(fhn, t, x, [b])
    with pytest.raises(ValueError, match="Basis objects"):
        ProfiledODE(fhn, t, x, [b, "b"])  # type: ignore[list-item]
    with pytest.raises(ValueError, match="share one domain"):
        ProfiledODE(fhn, t, x, [b, fdt.BSpline(domain=(0.0, 21.0), n_basis=10)])
    with pytest.raises(ValueError, match="lam must be"):
        ProfiledODE(fhn, t, x, b, lam=-1.0)
    with pytest.raises(ValueError, match="state_weights must be"):
        ProfiledODE(fhn, t, x, b, state_weights=[1.0, 0.0])
    with pytest.raises(ValueError, match="basis domain"):
        ProfiledODE(fhn, t + 1.0, x, b)
    with pytest.raises(ValueError, match="no state has any observation"):
        ProfiledODE(fhn, t, np.full_like(x, np.nan), b)
    with pytest.raises(ValueError, match="y must have shape"):
        ProfiledODE(fhn, t, x[:, :1], b)
    with pytest.raises(ValueError, match="y must have one entry"):
        ProfiledODE(fhn, t, [x[:, 0]], b)
    with pytest.raises(ValueError, match="differ"):
        ProfiledODE(fhn, t, [x[:-1, 0], None], b)
    with pytest.raises(ValueError, match="equal-length"):
        ProfiledODE(fhn, t, x, b, quadrature=([0.0, 1.0], [1.0]))
    with pytest.raises(ValueError, match="non-negative"):
        ProfiledODE(fhn, t, x, b, quadrature=([0.0, 1.0], [1.0, -1.0]))
    p = ProfiledODE(fhn, t, x, b)
    with pytest.raises(ValueError, match="coefficient sizes"):
        p.criterion([np.ones(3), np.ones(10)], FHN_TRUE)
    with pytest.raises(ValueError, match="expected 20 coefficients"):
        p.criterion(np.ones(5), FHN_TRUE)
    with pytest.raises(ValueError, match="theta must have 3"):
        p.criterion(np.ones(20), [1.0])


def test_one_dimensional_y_for_one_state() -> None:
    t = np.linspace(0, 1, 11)
    p = ProfiledODE(decay_model(), t, np.exp(-t), fdt.BSpline(n_basis=6))
    assert p.n_obs == 11


def test_non_spline_basis_uses_even_panels() -> None:
    t = np.linspace(0, 1, 11)
    p = ProfiledODE(decay_model(), t, np.exp(-t)[:, None], fdt.Monomial(n_basis=4))
    assert p.nodes.shape == (4 * 4 * 5,)
    assert float(p.weights.sum()) == pytest.approx(1.0)


def test_custom_quadrature() -> None:
    t = np.linspace(0, 1, 11)
    nodes, weights = np.polynomial.legendre.leggauss(20)
    p = ProfiledODE(
        decay_model(),
        t,
        np.exp(-t)[:, None],
        fdt.BSpline(n_basis=8),
        lam=1e2,
        quadrature=((nodes + 1) / 2, weights / 2),
    )
    inner = p.fit_states([1.0])
    assert inner.converged
    assert inner.sse < 1e-9


def test_linear_ode_inner_fit_is_exact_smoother() -> None:
    """For f = 0 the inner problem is a penalised least-squares smoother."""
    t = np.linspace(0, 1, 15)
    y = np.sin(3 * t)
    basis = fdt.BSpline(n_basis=9)
    zero = ODEModel(lambda x, t, th: 0.0 * x, 1, 1)
    p = ProfiledODE(zero, t, y[:, None], basis, lam=0.01)
    inner = p.fit_states([0.0])
    phi = np.asarray(basis(t))
    nodes, w = p.nodes, p.weights
    dphi = np.asarray(basis(nodes, 1))
    c = np.linalg.solve(phi.T @ phi + 0.01 * dphi.T @ (w[:, None] * dphi), phi.T @ y)
    np.testing.assert_allclose(inner.coefs[0], c, rtol=1e-10)
    s = phi @ np.linalg.solve(phi.T @ phi + 0.01 * dphi.T @ (w[:, None] * dphi), phi.T)
    assert inner.df == pytest.approx(np.trace(s), rel=1e-10)
    assert inner.gcv == pytest.approx(15 * inner.sse / (15 - inner.df) ** 2, rel=1e-12)


def test_dcoefs_dtheta_matches_finite_differences(fhn_data: Any, fhn_basis: Any) -> None:
    t, x = fhn_data
    rng = np.random.default_rng(4)
    y = x + 0.05 * rng.standard_normal(x.shape)
    p = ProfiledODE(fitzhugh_nagumo_model(), t, y, fhn_basis, lam=100.0)
    theta = np.array([0.25, 0.3, 2.8])
    inner = p.fit_states(theta)
    for j in range(3):
        h = 1e-6
        e = np.zeros(3)
        e[j] = h
        plus = np.concatenate(p.fit_states(theta + e, inner.coefs).coefs)
        minus = np.concatenate(p.fit_states(theta - e, inner.coefs).coefs)
        fd = (plus - minus) / (2 * h)
        np.testing.assert_allclose(inner.dcoefs_dtheta[:, j], fd, atol=1e-6 * np.abs(fd).max())


def test_inner_iteration_cap_reports_non_convergence(fhn_data: Any, fhn_basis: Any) -> None:
    t, x = fhn_data
    p = ProfiledODE(fitzhugh_nagumo_model(), t, x, fhn_basis, lam=100.0)
    inner = p.fit_states([1.0, 1.0, 1.0], max_iter=1)
    assert isinstance(inner, InnerFit)
    assert not inner.converged
    assert inner.n_iter == 1


def test_fhn_recovery_with_noise(fhn_data: Any, fhn_basis: Any) -> None:
    t, x = fhn_data
    rng = np.random.default_rng(20260906)
    y = x + 0.05 * rng.standard_normal(x.shape)
    fit = profile_ode(fitzhugh_nagumo_model(), t, y, fhn_basis, lam=1e3, theta0=[0.4, 0.4, 2.0])
    assert isinstance(fit, ProfileResult)
    assert fit.converged
    assert np.all(np.abs(fit.theta - FHN_TRUE) < 3 * fit.stderr)
    assert np.all(fit.stderr < 0.05)
    assert fit.param_names == ("a", "b", "c")
    assert fit.state_names == ("V", "R")
    assert fit.theta_path.shape == (fit.n_iter + 1, 3)
    assert np.all(np.diff(fit.sse_path) <= 1e-9 * fit.sse_path[0])
    assert fit.sigma2 == pytest.approx(fit.inner.sse / (2 * 201 - 3))
    states = fit(t)
    assert states.shape == (201, 2)
    assert np.sqrt(np.mean((states - x) ** 2)) < 0.03
    np.testing.assert_allclose(fit(t, 1)[:, 0], fit.states[0](t, 1)[:, 0])
    np.testing.assert_allclose(fit.cov, fit.cov.T)


def test_fhn_recovery_with_unobserved_state(fhn_data: Any, fhn_basis: Any) -> None:
    t, x = fhn_data
    rng = np.random.default_rng(7)
    y = np.column_stack([x[:, 0] + 0.05 * rng.standard_normal(201), np.full(201, np.nan)])
    fit = profile_ode(fitzhugh_nagumo_model(), t, y, fhn_basis, lam=1e3, theta0=[0.3, 0.3, 2.5])
    assert fit.converged
    assert np.all(np.abs(fit.theta - FHN_TRUE) < 3 * fit.stderr)
    # the unobserved R is reconstructed from the equations alone
    assert np.sqrt(np.mean((fit(t)[:, 1] - x[:, 1]) ** 2)) < 0.05


def test_cstr_recovery_from_temperature_only() -> None:
    full = cstr_model("all.cool.step")
    t = np.arange(0.0, 24.001, 0.25)
    x = full.simulate(t, [1.5965, 341.3754], list(CSTR_PARAMETERS.values()), breaks=range(4, 24, 4))
    rng = np.random.default_rng(11)
    y = x + rng.standard_normal(x.shape) * [0.02, 0.5]
    y[:, 0] = np.nan
    basis = fdt.BSpline(domain=(0.0, 24.0), breaks=np.arange(0.0, 24.001, 0.5).tolist())
    wt = [0.02, float(np.nanvar(y[:, 1]))]
    fit = profile_ode(
        cstr_model("all.cool.step", estimate=("kref", "EoverR")),
        t,
        y,
        basis,
        lam=[100.0, 100.0],
        theta0=[0.4, 0.8],
        state_weights=[1 / w for w in wt],
    )
    assert fit.converged
    truth = np.array([CSTR_PARAMETERS["kref"], CSTR_PARAMETERS["EoverR"]])
    assert np.all(np.abs(fit.theta - truth) < 3 * fit.stderr)


def test_fit_needs_more_observations_than_parameters() -> None:
    p = ProfiledODE(decay_model(), [0.0], [[1.0]], fdt.BSpline(n_basis=4))
    with pytest.raises(ValueError, match="more observations"):
        p.fit([1.0])


def test_exact_data_converges_by_offset() -> None:
    t = np.linspace(0.0, 1.0, 21)
    basis = fdt.BSpline(domain=(0.0, 1.0), breaks=np.linspace(0.0, 1.0, 11).tolist())
    p = ProfiledODE(decay_model(), t, (2 * np.exp(-0.7 * t))[:, None], basis, lam=1e4)
    fit = p.fit([0.2], coef0=np.zeros(p.n_coefs))
    assert fit.converged
    assert fit.theta[0] == pytest.approx(0.7, rel=1e-5)
    capped = p.fit([0.2], max_iter=1)
    assert not capped.converged
    assert capped.n_iter == 1
