"""Unit tests for :mod:`fdatools.registration`."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.exceptions import NotFittedError

from fdatools import Basis, BSpline, FData, Fourier, Monomial
from fdatools.registration import (
    AmpPhaseDecomposition,
    RegistrationResult,
    Registrator,
    _CurveProblem,
    _fine_grid,
    _minimise,
    _newton_step,
    _penalty_matrix,
    _warp_values,
    _WarpQuadrature,
    landmark_register,
    register,
)

DOMAIN = (0.0, 1.0)
GRID = np.linspace(*DOMAIN, 401)


def fit_curves(values: np.ndarray, basis: Any, grid: np.ndarray = GRID) -> FData:
    """Least-squares fit of sampled curves (columns of ``values``)."""
    return FData(np.linalg.lstsq(basis(grid), values, rcond=None)[0], basis)


def bumps(centres: list[float], width: float = 0.1, n_basis: int = 20) -> FData:
    """Gaussian bumps centred at ``centres``."""
    values = np.stack([np.exp(-(((GRID - c) / width) ** 2)) for c in centres], axis=1)
    return fit_curves(values, BSpline(domain=DOMAIN, n_basis=n_basis))


def peak_times(fd: FData) -> np.ndarray:
    fine = np.linspace(*fd.domain, 4001)
    return fine[np.argmax(np.asarray(fd(fine)), axis=0)]


# --------------------------------------------------------------------------- #
# derivatives of the criterion
# --------------------------------------------------------------------------- #


def _problem(criterion: str, periodic: bool, lam: float = 0.3) -> tuple[_CurveProblem, np.ndarray]:
    curve_basis: Basis
    if periodic:
        curve_basis = Fourier(domain=DOMAIN, n_basis=7)
        values = np.stack(
            [np.sin(2 * np.pi * (GRID - 0.07)) + 0.3 * np.cos(4 * np.pi * GRID)], axis=1
        )
        target = np.sin(2 * np.pi * GRID)
    else:
        curve_basis = BSpline(domain=DOMAIN, n_basis=12, order=5)
        values = np.stack([np.exp(-(((GRID - 0.45) / 0.15) ** 2))], axis=1)
        target = np.exp(-(((GRID - 0.55) / 0.15) ** 2))
    curve = fit_curves(values, curve_basis)
    grid = _fine_grid(DOMAIN, curve_basis.n_basis)
    goal = np.interp(grid, GRID, target)
    warp_basis = BSpline(domain=DOMAIN, n_basis=6)
    problem = _CurveProblem(
        curve=curve,
        target=goal,
        quadrature=_WarpQuadrature(warp_basis, grid),
        penalty=_penalty_matrix(warp_basis, 2, lam),
        lam=lam,
        criterion=criterion,
        periodic=periodic,
        has_curvature=True,
    )
    size = warp_basis.n_basis - 1 + (1 if periodic else 0)
    params = np.random.default_rng(7).normal(scale=0.3, size=size)
    return problem, params


@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
@pytest.mark.parametrize("periodic", [False, True])
def test_gradient_and_hessian_match_finite_differences(criterion: str, periodic: bool) -> None:
    problem, params = _problem(criterion, periodic)
    value, grad, hess = problem.evaluate(params)
    step = 1e-6
    num_grad = np.zeros_like(params)
    num_hess = np.zeros((params.size, params.size))
    for k in range(params.size):
        shift = np.zeros_like(params)
        shift[k] = step
        up = problem.evaluate(params + shift)
        down = problem.evaluate(params - shift)
        num_grad[k] = (up[0] - down[0]) / (2 * step)
        num_hess[:, k] = (up[1] - down[1]) / (2 * step)
    np.testing.assert_allclose(grad, num_grad, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(hess, num_hess, rtol=1e-4, atol=1e-6)
    assert np.isfinite(value)


def test_warp_quadrature_matches_exact_warp() -> None:
    basis = BSpline(domain=DOMAIN, n_basis=6)
    coefs = np.array([0.0, 0.4, -0.3, 0.8, 0.1, -0.5])
    grid = _fine_grid(DOMAIN, 10)
    h, h_c, _ = _WarpQuadrature(basis, grid).warp(coefs)
    np.testing.assert_allclose(h, _warp_values(basis, coefs[:, None], grid)[:, 0], atol=1e-13)
    step = 1e-6
    for k in range(coefs.size):
        bump = np.zeros_like(coefs)
        bump[k] = step
        up = _warp_values(basis, (coefs + bump)[:, None], grid)[:, 0]
        down = _warp_values(basis, (coefs - bump)[:, None], grid)[:, 0]
        np.testing.assert_allclose(h_c[:, k], (up - down) / (2 * step), atol=1e-8)


# --------------------------------------------------------------------------- #
# the optimiser
# --------------------------------------------------------------------------- #


def test_newton_step_reflects_negative_curvature() -> None:
    hess = np.diag([2.0, -4.0])
    grad = np.array([2.0, 4.0])
    step = _newton_step(hess, grad)
    np.testing.assert_allclose(step, [-1.0, -1.0])
    assert float(grad @ step) < 0.0


def test_newton_step_zero_hessian_is_finite() -> None:
    step = _newton_step(np.zeros((2, 2)), np.array([1.0, -1.0]))
    assert np.all(np.isfinite(step))


def test_minimise_quadratic_in_one_step() -> None:
    target = np.array([1.0, -2.0])

    def quadratic(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        d = p - target
        return float(d @ d), 2.0 * d, 2.0 * np.eye(2)

    params, value, used, converged = _minimise(quadratic, np.zeros(2), 10, 1e-12)
    np.testing.assert_allclose(params, target)
    assert converged
    assert used == 1
    assert value == pytest.approx(0.0)


def test_minimise_reports_non_convergence() -> None:
    def slow(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        # A wrong (too large) Hessian makes every Newton step short.
        return float(p @ p), 2.0 * p, 2e6 * np.eye(1)

    params, _, used, converged = _minimise(slow, np.ones(1), 3, 1e-12)
    assert not converged
    assert used == 3
    assert params[0] > 0.9


def test_minimise_stops_when_no_step_decreases() -> None:
    def ascent(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        # The gradient points the wrong way, so every step increases the value.
        return float(p @ p), -2.0 * p - 1.0, 2.0 * np.eye(1)

    params, _, used, converged = _minimise(ascent, np.ones(1), 5, 1e-12)
    assert not converged
    assert used == 0
    np.testing.assert_allclose(params, [1.0])


def test_minimise_empty_parameter_vector() -> None:
    def constant(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        return 1.0, np.zeros(0), np.zeros((0, 0))

    _, value, used, converged = _minimise(constant, np.zeros(0), 5, 1e-12)
    assert converged
    assert used == 0
    assert value == 1.0


def test_newton_step_non_finite_hessian_is_steepest_descent() -> None:
    grad = np.array([1.0, -2.0])
    np.testing.assert_array_equal(_newton_step(np.full((2, 2), np.nan), grad), -grad)
    np.testing.assert_array_equal(_newton_step(np.diag([np.inf, 1.0]), grad), -grad)


def test_newton_step_falls_back_when_eigh_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    from fdatools._backend import default_namespace

    def broken(_: Any) -> Any:
        raise np.linalg.LinAlgError("Eigenvalues did not converge")

    monkeypatch.setattr(default_namespace().linalg, "eigh", broken)
    grad = np.array([3.0, -1.0])
    np.testing.assert_array_equal(_newton_step(np.eye(2), grad), -grad)


def test_minimise_rejects_non_finite_steps() -> None:
    def cliff(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        # Finite below 1.5, overflowing beyond: the full Newton step (to 2.0)
        # lands on the non-finite side and must be shortened, not accepted.
        if p[0] > 1.5:
            return float("nan"), np.full(1, np.nan), np.full((1, 1), np.nan)
        d = p - 2.0
        return float(d @ d), 2.0 * d, 2.0 * np.eye(1)

    params, value, _, _ = _minimise(cliff, np.zeros(1), 50, 1e-12)
    assert np.isfinite(value)
    assert 1.0 <= params[0] <= 1.5


def test_minimise_rejects_a_non_finite_start() -> None:
    def broken(p: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        return float("inf"), np.zeros(1), np.eye(1)

    with pytest.raises(ValueError, match="not finite"):
        _minimise(broken, np.zeros(1), 5, 1e-12)


def test_warp_quadrature_survives_huge_latent_coefficients() -> None:
    # exp(900) overflows float64; the warp only depends on W up to a constant,
    # so the moments are taken relative to the largest latent coefficient.
    basis = BSpline(domain=DOMAIN, n_basis=5)
    grid = _fine_grid(DOMAIN, 10)
    coefs = np.array([0.0, 300.0, 900.0, 600.0, 0.0])
    h, grad, hess = _WarpQuadrature(basis, grid).warp(coefs)
    assert np.all(np.isfinite(h))
    assert np.all(np.isfinite(grad))
    assert np.all(np.isfinite(hess))
    shifted, _, _ = _WarpQuadrature(basis, grid).warp(coefs - 900.0)
    np.testing.assert_allclose(h, shifted, rtol=0, atol=1e-12)
    assert np.all(np.diff(h) >= 0.0)
    values = _warp_values(basis, coefs[:, None], grid)
    np.testing.assert_allclose(values[:, 0], h, rtol=0, atol=1e-10)


# --------------------------------------------------------------------------- #
# continuous registration
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
def test_register_aligns_shifted_bumps(criterion: str) -> None:
    centres = [0.42, 0.5, 0.58]
    fd = bumps(centres)
    res = register(fd, criterion=criterion, warp_basis=BSpline(domain=DOMAIN, n_basis=5), lam=1e-4)
    before = np.ptp(peak_times(fd))
    after = np.ptp(peak_times(res.registered))
    assert after < 0.2 * before
    assert res.criterion is not None
    assert res.n_iter is not None
    assert np.all(res.n_iter < 50)
    assert res.target is not None
    assert res.target.n_curves == 1
    np.testing.assert_array_equal(res.shift, np.zeros(3))


def test_register_lowers_the_criterion_from_the_identity() -> None:
    fd = bumps([0.45, 0.55])
    options: dict[str, Any] = {"warp_basis": BSpline(domain=DOMAIN, n_basis=4), "lam": 1e-3}
    start = register(fd, max_iter=0, **options)
    end = register(fd, **options)
    assert start.criterion is not None
    assert end.criterion is not None
    assert np.all(end.criterion < start.criterion)


def test_register_default_warp_basis_is_linear() -> None:
    fd = bumps([0.45, 0.55])
    res = register(fd)
    assert isinstance(res.latent.basis, BSpline)
    assert res.latent.basis.n_basis == 2
    assert res.latent.basis.order == 2
    np.testing.assert_array_equal(np.asarray(res.latent.coefs)[0], [0.0, 0.0])


def test_register_zero_iterations_is_the_identity() -> None:
    fd = bumps([0.4, 0.6])
    res = register(fd, max_iter=0)
    np.testing.assert_allclose(res.registered(GRID), fd(GRID), atol=1e-10)
    np.testing.assert_allclose(res.warp(GRID), np.tile(GRID[:, None], (1, 2)), atol=1e-12)
    np.testing.assert_array_equal(res.n_iter, [0, 0])


def test_register_per_curve_targets() -> None:
    fd = bumps([0.45, 0.55])
    targets = bumps([0.5, 0.5])
    res = register(fd, targets, warp_basis=BSpline(domain=DOMAIN, n_basis=4), lam=1e-4)
    np.testing.assert_allclose(peak_times(res.registered), [0.5, 0.5], atol=0.02)


def test_register_init_is_normalised_and_broadcast() -> None:
    fd = bumps([0.45, 0.55])
    basis = BSpline(domain=DOMAIN, n_basis=4)
    init = np.array([1.0, 1.2, 0.8, 1.1])
    res = register(fd, warp_basis=basis, init=init, max_iter=0)
    coefs = np.asarray(res.latent.coefs)
    np.testing.assert_allclose(coefs, np.tile((init - init[0])[:, None], (1, 2)))
    via_fd = register(fd, init=FData(init, basis), max_iter=0)
    np.testing.assert_allclose(np.asarray(via_fd.latent.coefs), coefs)


def test_register_warns_when_not_converged() -> None:
    fd = bumps([0.3, 0.7])
    with pytest.warns(RuntimeWarning, match="did not converge"):
        register(fd, warp_basis=BSpline(domain=DOMAIN, n_basis=6), lam=1e-6, max_iter=1)


def test_register_low_order_curves_use_gauss_newton() -> None:
    basis = BSpline(domain=DOMAIN, n_basis=30, order=2)
    values = np.stack([np.sin(np.pi * GRID**p) for p in (0.9, 1.1)], axis=1)
    fd = fit_curves(values, basis)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = register(fd, criterion="least_squares", max_iter=200)
    assert res.criterion is not None
    assert np.all(np.isfinite(res.criterion))


def test_register_accepts_non_spline_curve_bases() -> None:
    fd = fit_curves(np.stack([GRID**2, GRID**2.2], axis=1), Monomial(domain=DOMAIN, n_basis=5))
    res = register(fd, criterion="least_squares")
    assert res.registered.basis == fd.basis


def test_register_zero_curves_are_stationary() -> None:
    fd = FData(np.zeros((6, 2)), BSpline(domain=DOMAIN, n_basis=6))
    res = register(fd)
    assert res.n_iter is not None
    np.testing.assert_array_equal(res.n_iter, [0, 0])
    assert res.criterion is not None
    np.testing.assert_allclose(res.criterion, [0.0, 0.0])


def test_periodic_registration_recovers_shifts() -> None:
    basis = Fourier(domain=DOMAIN, n_basis=7)
    offsets = [-0.06, 0.0, 0.06]
    values = np.stack([np.sin(2 * np.pi * (GRID - s)) for s in offsets], axis=1)
    fd = fit_curves(values, basis)
    target = fit_curves(np.sin(2 * np.pi * GRID)[:, None], basis)
    res = register(fd, target, periodic=True, criterion="least_squares")
    np.testing.assert_allclose(res.shift, offsets, atol=1e-6)
    np.testing.assert_allclose(res.registered(GRID), target(GRID) * np.ones((1, 3)), atol=1e-6)
    warps = res.warp_values(GRID)
    np.testing.assert_allclose(warps - GRID[:, None], np.tile(offsets, (GRID.size, 1)), atol=1e-6)


def test_periodic_registration_with_spline_curves_wraps() -> None:
    basis = BSpline(domain=DOMAIN, n_basis=15)
    values = np.stack([np.sin(2 * np.pi * (GRID - s)) for s in (-0.05, 0.05)], axis=1)
    fd = fit_curves(values, basis)
    res = register(fd, periodic=True, init_shift=0.01, criterion="least_squares")
    assert np.all(np.isfinite(np.asarray(res.registered.coefs)))
    assert res.shift[0] < res.shift[1]


# --------------------------------------------------------------------------- #
# landmark registration
# --------------------------------------------------------------------------- #


def test_landmarks_move_to_their_targets() -> None:
    fd = bumps([0.4, 0.5, 0.62])
    marks = peak_times(fd)
    res = register(fd, landmarks=marks)
    target = float(np.mean(marks))
    np.testing.assert_allclose(res.warp_values(np.array([target]))[0], marks, atol=1e-3)
    np.testing.assert_allclose(peak_times(res.registered), target, atol=2e-3)
    assert res.warp_inverse is not None
    np.testing.assert_allclose(np.asarray(res.latent.coefs).mean(axis=0), 0.0, atol=1e-14)


def test_landmark_inverse_composes_to_identity() -> None:
    fd = bumps([0.35, 0.5, 0.6])
    res = register(fd, landmarks=[[0.3, 0.7], [0.45, 0.8], [0.35, 0.75]])
    assert res.warp_inverse is not None
    values = res.warp_values(GRID)
    for i in range(3):
        back = np.asarray(res.warp_inverse[i](np.clip(values[:, i], *DOMAIN)))[:, 0]
        np.testing.assert_allclose(back, GRID, atol=1e-6)


def test_landmark_register_wrapper_and_explicit_targets() -> None:
    fd = bumps([0.4, 0.6])
    basis = BSpline(domain=DOMAIN, n_basis=6)
    res = landmark_register(fd, [0.4, 0.6], [0.5], warp_basis=basis, lam=1e-3)
    assert res.latent.basis == basis
    np.testing.assert_allclose(res.warp_values(np.array([0.5]))[0], [0.4, 0.6], atol=5e-3)
    assert res.target is None
    assert res.criterion is None


# --------------------------------------------------------------------------- #
# warps are monotone (property)
# --------------------------------------------------------------------------- #


@settings(max_examples=40, deadline=None)
@given(
    st.lists(st.floats(-3.0, 3.0, allow_nan=False), min_size=5, max_size=5),
    st.floats(0.5, 20.0),
)
def test_warps_are_strictly_increasing_bijections(coefs: list[float], width: float) -> None:
    basis = BSpline(domain=(1.0, 1.0 + width), n_basis=5)
    grid = np.linspace(1.0, 1.0 + width, 257)
    values = _warp_values(basis, np.asarray(coefs)[:, None], grid)[:, 0]
    assert values[0] == pytest.approx(1.0)
    assert values[-1] == pytest.approx(1.0 + width)
    assert np.all(np.diff(values) > 0.0)


@settings(max_examples=15, deadline=None)
@given(st.lists(st.floats(0.3, 0.7), min_size=2, max_size=4))
def test_registration_warps_are_monotone(centres: list[float]) -> None:
    fd = bumps(centres)
    res = register(fd, warp_basis=BSpline(domain=DOMAIN, n_basis=5), lam=1e-3)
    assert np.all(np.diff(res.warp_values(GRID), axis=0) > 0.0)


# --------------------------------------------------------------------------- #
# amplitude / phase decomposition
# --------------------------------------------------------------------------- #


def _result(x: FData, y: FData, h: FData) -> RegistrationResult:
    latent = FData(np.zeros((2, x.n_curves)), BSpline(domain=x.domain, n_basis=2, order=2))
    return RegistrationResult(
        registered=y, warp=h, unregistered=x, latent=latent, shift=np.zeros(x.n_curves)
    )


def test_decompose_known_values() -> None:
    basis = Monomial(domain=DOMAIN, n_basis=3)
    lines = FData(np.array([[0.0, 0.0], [1.0, -1.0], [0.0, 0.0]]), basis)
    identity = FData(np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]]), basis)
    amp, phase, rsq, c = _result(lines, lines, identity).decompose()
    # Trapezoidal rule on 201 points: the integral of t^2 is 1/3 + 1/240000.
    assert amp == pytest.approx(1.0 / 3.0 + 1.0 / 240000.0, rel=1e-12)
    assert phase == pytest.approx(0.0, abs=1e-15)
    assert rsq == pytest.approx(0.0, abs=1e-12)
    assert c == 1.0


def test_decompose_constant_uses_sample_covariance() -> None:
    basis = Monomial(domain=DOMAIN, n_basis=3)
    levels = FData(np.array([[1.0, 2.0], [0.0, 0.0], [0.0, 0.0]]), basis)
    warps = FData(np.array([[0.0, 0.0], [1.0, 1.0], [0.1, -0.1]]), basis)
    decomposition = _result(levels, levels, warps).decompose()
    assert isinstance(decomposition, AmpPhaseDecomposition)
    assert decomposition.c == pytest.approx(0.88, rel=1e-12)
    assert decomposition.amp_mse == pytest.approx(0.22, rel=1e-12)
    assert decomposition.phase_mse == pytest.approx(-0.27, rel=1e-12)
    assert decomposition.rsq == pytest.approx(5.4, rel=1e-12)


def test_decompose_sub_interval_and_errors() -> None:
    fd = bumps([0.45, 0.55])
    res = register(fd, landmarks=[0.45, 0.55])
    whole = res.decompose()
    part = res.decompose((0.2, 0.8))
    assert part.amp_mse != whole.amp_mse
    with pytest.raises(ValueError, match="inside"):
        res.decompose((0.5, 1.5))
    single = _result(fd[0], fd[0], res.warp[0])
    with pytest.raises(ValueError, match="two curves"):
        single.decompose()


def _growth_accelerations() -> FData:
    from fdatools.datasets import load_growth
    from fdatools.smoothing import smooth

    growth = load_growth()
    age = np.asarray(growth.age, dtype=float)
    basis = BSpline(domain=(float(age[0]), float(age[-1])), order=6, breaks=age.tolist())
    heights = np.asarray(growth.hgtf, dtype=float)[:, :10]
    return smooth(heights, age, basis=basis, lam=0.01, penalty=4).fd.derivative(2)


@pytest.mark.parametrize("n_warp", [6, 8])
def test_register_unpenalised_eigen_on_growth_accelerations(n_warp: int) -> None:
    # Regression: lam=0 with the eigen criterion used to crash with a raw
    # numpy LinAlgError ("Eigenvalues did not converge") after exp(W) overflowed.
    acc = _growth_accelerations()
    warp_basis = BSpline(domain=acc.domain, n_basis=n_warp)
    with warnings.catch_warnings():
        warnings.filterwarnings("error", category=RuntimeWarning, message="(?!registration)")
        warnings.filterwarnings("ignore", "registration did not converge", RuntimeWarning)
        res = register(acc, warp_basis=warp_basis, lam=0.0)
    assert res.criterion is not None
    assert np.all(np.isfinite(res.criterion))
    assert np.all(np.isfinite(np.asarray(res.latent.coefs)))
    assert np.all(np.isfinite(np.asarray(res.registered.coefs)))
    t = np.linspace(*acc.domain, 501)
    warps = res.warp_values(t)
    assert np.all(np.isfinite(warps))
    assert np.all(np.diff(warps, axis=0) >= 0.0)
    start = register(acc, warp_basis=warp_basis, lam=0.0, max_iter=0)
    assert start.criterion is not None
    assert np.all(res.criterion <= start.criterion + 1e-12)


def test_register_rejects_non_finite_input() -> None:
    fd = bumps([0.45, 0.55])
    coefs = np.array(fd.coefs, dtype=float)
    coefs[3, 1] = np.nan
    bad = FData(coefs, fd.basis)
    with pytest.raises(ValueError, match="finite"):
        register(bad)
    with pytest.raises(ValueError, match="finite"):
        register(fd, bad[1])
    with pytest.raises(ValueError, match="finite"):
        register(fd, init=np.array([0.0, np.inf]))
    with pytest.raises(ValueError, match="finite"):
        register(fd, periodic=True, init_shift=np.nan)
    with pytest.raises(ValueError, match="finite"):
        register(bad, landmarks=[0.4, 0.5])
    with pytest.raises(ValueError, match="finite"):
        register(fd, landmarks=[0.4, np.nan])


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #


def test_register_rejects_bad_arguments() -> None:
    fd = bumps([0.45, 0.55])
    wide = FData(np.zeros((20, 2, 2)), BSpline(domain=DOMAIN, n_basis=20))
    with pytest.raises(ValueError, match="variables per curve"):
        register(fd, wide)
    with pytest.raises(ValueError, match="variables per curve"):
        register(wide, fd)
    with pytest.raises(ValueError, match="target must hold"):
        register(fd, bumps([0.4, 0.5, 0.6]))
    with pytest.raises(ValueError, match="criterion"):
        register(fd, criterion="median")
    with pytest.raises(ValueError, match="non-negative"):
        register(fd, lam=-1.0)
    with pytest.raises(ValueError, match="max_iter"):
        register(fd, max_iter=-1)
    with pytest.raises(ValueError, match="B-spline"):
        register(fd, warp_basis=Fourier(domain=DOMAIN, n_basis=3))
    with pytest.raises(ValueError, match="domain"):
        register(fd, warp_basis=BSpline(domain=(0.0, 2.0), n_basis=4))
    with pytest.raises(ValueError, match="rows"):
        register(fd, warp_basis=BSpline(domain=DOMAIN, n_basis=4), init=np.zeros(3))
    with pytest.raises(ValueError, match="columns"):
        register(fd, warp_basis=BSpline(domain=DOMAIN, n_basis=4), init=np.zeros((4, 3)))


def test_landmarks_reject_bad_arguments() -> None:
    fd = bumps([0.45, 0.55])
    with pytest.raises(ValueError, match="one row per curve"):
        register(fd, landmarks=[0.4, 0.5, 0.6])
    with pytest.raises(ValueError, match="strictly inside"):
        register(fd, landmarks=[0.0, 0.5])
    with pytest.raises(ValueError, match="increasing"):
        register(fd, landmarks=[[0.6, 0.4], [0.5, 0.7]])
    with pytest.raises(ValueError, match="entries"):
        register(fd, landmarks=[0.4, 0.5], target_landmarks=[0.4, 0.5])
    with pytest.raises(ValueError, match="strictly inside"):
        register(fd, landmarks=[0.4, 0.5], target_landmarks=[1.5])
    with pytest.raises(ValueError, match="non-negative"):
        register(fd, landmarks=[0.4, 0.5], lam=-1.0)
    with pytest.raises(ValueError, match="B-spline"):
        register(fd, landmarks=[0.4, 0.5], warp_basis=Fourier(domain=DOMAIN, n_basis=3))


# --------------------------------------------------------------------------- #
# scikit-learn estimator
# --------------------------------------------------------------------------- #


def test_registrator_fit_transform_fdata() -> None:
    fd = bumps([0.45, 0.5, 0.55])
    est = Registrator(warp_basis=BSpline(domain=DOMAIN, n_basis=4), lam=1e-4)
    coefs = est.fit_transform(fd)
    assert coefs.shape == (3, fd.basis.n_basis)
    assert est.n_features_in_ == fd.basis.n_basis
    np.testing.assert_allclose(coefs, np.asarray(est.result_.registered.coefs).T, atol=1e-10)
    again = est.transform(fd[:2])
    np.testing.assert_allclose(again, coefs[:2], atol=1e-8)


def test_registrator_accepts_coefficient_rows() -> None:
    rng = np.random.default_rng(3)
    rows = rng.normal(size=(4, 8))
    est = Registrator(criterion="least_squares").fit(rows)
    assert est.transform(rows).shape == (4, 8)
    basis = BSpline(domain=DOMAIN, n_basis=8)
    assert Registrator(basis=basis).fit(rows).target_.basis == basis


def test_registrator_errors() -> None:
    with pytest.raises(NotFittedError):
        Registrator().transform(np.zeros((2, 8)))
    est = Registrator().fit(bumps([0.45, 0.55]))
    with pytest.raises(ValueError, match="expected"):
        est.transform(FData(np.zeros((5, 1)), BSpline(domain=DOMAIN, n_basis=5)))
    with pytest.raises(ValueError, match="columns"):
        Registrator(basis=BSpline(domain=DOMAIN, n_basis=5)).fit(np.zeros((2, 8)))


# --------------------------------------------------------------------------- #
# multivariate curves
# --------------------------------------------------------------------------- #


def two_variable_bumps(centres: list[float], n_basis: int = 20) -> FData:
    """Curves with two variables that share their phase: a bump and a wave."""
    first = np.stack([np.exp(-(((GRID - c) / 0.1) ** 2)) for c in centres], axis=1)
    second = np.stack(
        [np.sin(2 * np.pi * (GRID - c)) * first[:, i] for i, c in enumerate(centres)], axis=1
    )
    basis = BSpline(domain=DOMAIN, n_basis=n_basis)
    stacked = np.stack([first, second], axis=2).reshape(GRID.size, -1)
    coefs = np.linalg.lstsq(basis(GRID), stacked, rcond=None)[0]
    return FData(coefs.reshape(n_basis, len(centres), 2), basis)


def _multi_problem(
    criterion: str, periodic: bool, weights: tuple[float, ...]
) -> tuple[_CurveProblem, np.ndarray]:
    problem, params = _problem(criterion, periodic)
    curve = problem.curve
    coefs = np.stack([curve.coefs, 0.5 * curve.coefs + 0.2], axis=2)
    target = np.stack([problem.target, 0.8 * problem.target - 0.1], axis=1)
    multi = _CurveProblem(
        curve=FData(coefs, curve.basis),
        target=target,
        quadrature=problem.quadrature,
        penalty=problem.penalty,
        lam=problem.lam,
        criterion=criterion,
        periodic=periodic,
        has_curvature=True,
        weights=weights,
    )
    return multi, params


@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
@pytest.mark.parametrize("periodic", [False, True])
def test_multivariate_gradient_and_hessian_match_finite_differences(
    criterion: str, periodic: bool
) -> None:
    problem, params = _multi_problem(criterion, periodic, (1.0, 0.7))
    value, grad, hess = problem.evaluate(params)
    step = 1e-6
    num_grad = np.zeros_like(params)
    num_hess = np.zeros((params.size, params.size))
    for k in range(params.size):
        shift = np.zeros_like(params)
        shift[k] = step
        up = problem.evaluate(params + shift)
        down = problem.evaluate(params - shift)
        num_grad[k] = (up[0] - down[0]) / (2 * step)
        num_hess[:, k] = (up[1] - down[1]) / (2 * step)
    np.testing.assert_allclose(grad, num_grad, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(hess, num_hess, rtol=1e-4, atol=1e-6)
    assert np.isfinite(value)


@settings(max_examples=25, deadline=None)
@given(
    st.floats(min_value=0.0, max_value=3.0),
    st.floats(min_value=0.01, max_value=3.0),
    st.sampled_from(["eigen", "least_squares"]),
)
def test_multivariate_criterion_is_the_weighted_sum(w1: float, w2: float, criterion: str) -> None:
    problem, params = _multi_problem(criterion, False, (w1, w2))
    total = problem.evaluate(params)
    parts = []
    for var in range(2):
        single = _CurveProblem(
            curve=FData(problem.curve.coefs[:, :, var], problem.curve.basis),
            target=problem.target[:, var],
            quadrature=problem.quadrature,
            penalty=problem.penalty,
            lam=0.0,
            criterion=criterion,
            periodic=False,
            has_curvature=True,
        )
        parts.append(single.evaluate(params))
    pen = problem.penalty[1:, 1:]
    rough = problem.lam * float(params @ pen @ params)
    expected = w1 * parts[0][0] + w2 * parts[1][0] + rough
    assert total[0] == pytest.approx(expected, rel=1e-12, abs=1e-14)
    np.testing.assert_allclose(
        total[1] - 2 * problem.lam * pen @ params,
        w1 * parts[0][1] + w2 * parts[1][1],
        rtol=1e-10,
        atol=1e-12,
    )


@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
def test_multivariate_registration_aligns_every_variable(criterion: str) -> None:
    centres = [0.44, 0.5, 0.56]
    fd = two_variable_bumps(centres)
    res = register(fd, criterion=criterion, warp_basis=BSpline(domain=DOMAIN, n_basis=5), lam=1e-4)
    assert res.registered.coefs.shape == fd.coefs.shape
    assert res.latent.coefs.shape == (5, 3)
    assert res.warp.n_curves == 3
    assert res.target is not None
    assert res.target.n_vars == 2
    first = FData(res.registered.coefs[:, :, 0], fd.basis)
    assert np.ptp(peak_times(first)) < 0.3 * np.ptp(centres)
    before = FData(fd.coefs[:, :, 1], fd.basis)
    after = FData(res.registered.coefs[:, :, 1], fd.basis)
    spread_before = float(np.mean(np.var(np.asarray(before(GRID)), axis=1)))
    spread_after = float(np.mean(np.var(np.asarray(after(GRID)), axis=1)))
    assert spread_after < 0.3 * spread_before


def test_var_weights_choose_the_variables_that_drive_the_warp() -> None:
    fd = two_variable_bumps([0.45, 0.55])
    zeros = np.zeros_like(fd.coefs[:, :, 0])
    flat = FData(np.stack([zeros, fd.coefs[:, :, 1]], axis=2), fd.basis)
    wbasis = BSpline(domain=DOMAIN, n_basis=4)
    ignored = register(flat, criterion="least_squares", warp_basis=wbasis, var_weights=[1, 0])
    np.testing.assert_allclose(ignored.latent.coefs, 0.0, atol=1e-14)
    np.testing.assert_allclose(ignored.registered.coefs, flat.coefs, atol=1e-10)
    driven = register(flat, criterion="least_squares", warp_basis=wbasis, lam=1e-4)
    assert float(np.max(np.abs(driven.latent.coefs))) > 0.1
    alone = register(
        FData(flat.coefs[:, :, 1], fd.basis),
        criterion="least_squares",
        warp_basis=wbasis,
        lam=1e-4,
    )
    np.testing.assert_allclose(driven.latent.coefs, alone.latent.coefs, atol=1e-8)


def test_var_weights_scale_the_criterion() -> None:
    fd = two_variable_bumps([0.46, 0.54])
    options: dict[str, Any] = {
        "criterion": "least_squares",
        "warp_basis": BSpline(domain=DOMAIN, n_basis=4),
    }
    once = register(fd, **options)
    twice = register(fd, var_weights=[2.0, 2.0], **options)
    np.testing.assert_allclose(twice.latent.coefs, once.latent.coefs, atol=1e-7)
    assert once.criterion is not None
    assert twice.criterion is not None
    np.testing.assert_allclose(twice.criterion, 2.0 * np.asarray(once.criterion), rtol=1e-8)


def test_single_variable_axis_matches_univariate() -> None:
    fd = bumps([0.45, 0.5, 0.55])
    column = FData(fd.coefs[:, :, None], fd.basis)
    options: dict[str, Any] = {"warp_basis": BSpline(domain=DOMAIN, n_basis=4), "lam": 1e-3}
    plain = register(fd, **options)
    stacked = register(column, **options)
    np.testing.assert_array_equal(stacked.latent.coefs, plain.latent.coefs)
    np.testing.assert_array_equal(stacked.registered.coefs[:, :, 0], plain.registered.coefs)


def test_multivariate_per_curve_targets_and_periodic_shifts() -> None:
    basis = Fourier(domain=DOMAIN, n_basis=7)
    shifts = [-0.06, 0.0, 0.06]
    first = np.stack([np.sin(2 * np.pi * (GRID - s)) for s in shifts], axis=1)
    second = np.stack([np.cos(4 * np.pi * (GRID - s)) for s in shifts], axis=1)
    both = np.concatenate([first, second], axis=1)
    coefs = np.linalg.lstsq(basis(GRID), both, rcond=None)[0]
    fd = FData(np.stack([coefs[:, :3], coefs[:, 3:]], axis=2), basis)
    target = FData(np.repeat(fd.coefs[:, 1:2, :], 3, axis=1), basis)
    res = register(
        fd,
        target,
        criterion="least_squares",
        periodic=True,
        warp_basis=BSpline(domain=DOMAIN, n_basis=3, order=3),
        lam=1e-2,
    )
    np.testing.assert_allclose(res.shift, shifts, atol=1e-3)
    np.testing.assert_allclose(res.apply(fd).coefs, res.registered.coefs, atol=1e-12)


def test_var_weights_are_validated() -> None:
    fd = two_variable_bumps([0.45, 0.55])
    with pytest.raises(ValueError, match="2 entries"):
        register(fd, var_weights=[1.0])
    with pytest.raises(ValueError, match="non-negative"):
        register(fd, var_weights=[1.0, -1.0])
    with pytest.raises(ValueError, match="positive"):
        register(fd, var_weights=[0.0, 0.0])
    with pytest.raises(ValueError, match="finite"):
        register(fd, var_weights=[1.0, np.nan])
    # Landmark registration ignores the weights.
    register(fd, landmarks=[0.45, 0.55], var_weights=[1.0])


def test_multivariate_landmarks_share_one_warp() -> None:
    fd = two_variable_bumps([0.45, 0.5, 0.55])
    marks = [[0.45], [0.5], [0.55]]
    res = register(fd, landmarks=marks)
    assert res.registered.coefs.shape == fd.coefs.shape
    assert res.warp_inverse is not None
    assert res.warp_inverse.n_curves == 3
    for var in range(2):
        single = register(FData(fd.coefs[:, :, var], fd.basis), landmarks=marks)
        np.testing.assert_array_equal(res.latent.coefs, single.latent.coefs)
        np.testing.assert_allclose(
            res.registered.coefs[:, :, var], single.registered.coefs, atol=1e-13
        )


def test_decompose_multivariate_sums_over_variables() -> None:
    fd = bumps([0.45, 0.5, 0.55])
    res = register(fd, landmarks=[0.45, 0.5, 0.55])
    uni = res.decompose()
    doubled = RegistrationResult(
        registered=FData(np.stack([res.registered.coefs] * 2, axis=2), fd.basis),
        warp=res.warp,
        unregistered=FData(np.stack([fd.coefs] * 2, axis=2), fd.basis),
        latent=res.latent,
        shift=res.shift,
    )
    multi = doubled.decompose()
    assert multi.amp_mse == pytest.approx(2.0 * uni.amp_mse, rel=1e-12)
    assert multi.phase_mse == pytest.approx(2.0 * uni.phase_mse, rel=1e-10)
    assert multi.rsq == pytest.approx(uni.rsq, rel=1e-10)
    assert multi.c == pytest.approx(uni.c, rel=1e-12)


def test_registrator_rejects_multivariate_curves() -> None:
    with pytest.raises(ValueError, match="univariate"):
        Registrator().fit(two_variable_bumps([0.45, 0.55]))


# --------------------------------------------------------------------------- #
# applying warps to new curves (register.newfd)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("criterion", ["eigen", "least_squares"])
def test_apply_reproduces_the_registered_curves(criterion: str) -> None:
    fd = bumps([0.45, 0.5, 0.55])
    res = register(fd, criterion=criterion, warp_basis=BSpline(domain=DOMAIN, n_basis=4), lam=1e-3)
    np.testing.assert_array_equal(res.apply(fd).coefs, res.registered.coefs)
    multi = two_variable_bumps([0.45, 0.55])
    mres = register(multi, landmarks=[0.45, 0.55])
    np.testing.assert_array_equal(mres.apply(multi).coefs, mres.registered.coefs)


def test_apply_warps_derivatives_and_other_bases() -> None:
    fd = bumps([0.45, 0.5, 0.55])
    res = register(fd, landmarks=[0.45, 0.5, 0.55])
    slope = fd.derivative()
    velocity = res.apply(slope)
    assert velocity.basis == slope.basis
    # By definition: x_i'(h_i(t)) on the registration grid, fitted in the basis.
    grid = _fine_grid(DOMAIN, slope.basis.n_basis)
    h = np.asarray(res.warp_values(grid))
    values = np.stack([np.asarray(slope[i](h[:, i]))[:, 0] for i in range(3)], axis=1)
    expected = np.linalg.lstsq(slope.basis(grid), values, rcond=None)[0]
    np.testing.assert_allclose(velocity.coefs, expected, rtol=1e-10, atol=1e-10)
    other = fit_curves(np.asarray(fd(GRID)), BSpline(domain=DOMAIN, n_basis=30))
    assert res.apply(other).basis.n_basis == 30


def test_apply_periodic_shift_wraps_around() -> None:
    basis = Fourier(domain=DOMAIN, n_basis=5)
    shifts = [-0.05, 0.05]
    values = np.stack([np.sin(2 * np.pi * (GRID - s)) for s in shifts], axis=1)
    fd = fit_curves(values, basis)
    res = register(
        fd,
        criterion="least_squares",
        periodic=True,
        warp_basis=BSpline(domain=DOMAIN, n_basis=2, order=2),
    )
    np.testing.assert_allclose(res.apply(fd).coefs, res.registered.coefs, atol=1e-12)
    waves = np.stack([np.cos(2 * np.pi * (GRID - s)) for s in shifts], axis=1)
    aligned = np.asarray(res.apply(fit_curves(waves, basis))(GRID))
    np.testing.assert_allclose(aligned[:, 0], aligned[:, 1], atol=1e-6)


def test_apply_rejects_bad_curves() -> None:
    fd = bumps([0.45, 0.55])
    res = register(fd, landmarks=[0.45, 0.55])
    with pytest.raises(ValueError, match="one per warp"):
        res.apply(bumps([0.45, 0.5, 0.55]))
    with pytest.raises(ValueError, match="domain"):
        res.apply(FData(np.ones((4, 2)), BSpline(domain=(0.0, 2.0), n_basis=4)))
    with pytest.raises(ValueError, match="finite"):
        res.apply(FData(np.full((20, 2), np.nan), fd.basis))


# --------------------------------------------------------------------------- #
# multivariate PyTorch input
# --------------------------------------------------------------------------- #


def test_multivariate_torch_matches_numpy() -> None:
    torch = pytest.importorskip("torch")
    fd = two_variable_bumps([0.46, 0.54])
    coefs = torch.tensor(np.asarray(fd.coefs), dtype=torch.float64, requires_grad=True)
    options: dict[str, Any] = {"warp_basis": BSpline(domain=DOMAIN, n_basis=4), "lam": 1e-3}
    plain = register(fd, **options)
    tensor = register(FData(coefs, fd.basis), **options)
    assert isinstance(tensor.registered.coefs, torch.Tensor)
    assert tuple(tensor.registered.coefs.shape) == fd.coefs.shape
    np.testing.assert_allclose(
        tensor.latent.coefs.detach().numpy(), plain.latent.coefs, rtol=1e-7, atol=1e-9
    )
    np.testing.assert_allclose(
        tensor.registered.coefs.detach().numpy(), plain.registered.coefs, rtol=1e-6, atol=1e-8
    )
    tensor.registered(GRID).sum().backward()
    assert coefs.grad is not None
    assert tuple(coefs.grad.shape) == fd.coefs.shape
    marks = register(FData(coefs.detach(), fd.basis), landmarks=[0.46, 0.54])
    np.testing.assert_allclose(
        marks.registered.coefs.numpy(),
        register(fd, landmarks=[0.46, 0.54]).registered.coefs,
        atol=1e-12,
    )


def test_multivariate_autograd_objective_matches_analytic() -> None:
    pytest.importorskip("torch")
    from fdatools._internal.registration_torch import AutogradObjective

    for criterion in ("eigen", "least_squares"):
        problem, params = _multi_problem(criterion, True, (0.6, 1.3))
        auto = AutogradObjective(problem)(params)
        exact = problem.evaluate(params)
        assert auto[0] == pytest.approx(exact[0], rel=1e-12)
        np.testing.assert_allclose(auto[1], exact[1], rtol=1e-9, atol=1e-11)
        np.testing.assert_allclose(auto[2], exact[2], rtol=1e-8, atol=1e-9)


def test_apply_to_tensor_curves_is_differentiable() -> None:
    torch = pytest.importorskip("torch")
    fd = two_variable_bumps([0.45, 0.55])
    res = register(fd, landmarks=[0.45, 0.55])
    coefs = torch.tensor(np.asarray(fd.coefs), dtype=torch.float64, requires_grad=True)
    out = res.apply(FData(coefs, fd.basis))
    assert isinstance(out.coefs, torch.Tensor)
    np.testing.assert_allclose(out.coefs.detach().numpy(), res.registered.coefs, atol=1e-12)
    out(GRID).sum().backward()
    assert coefs.grad is not None
    assert bool(torch.all(torch.isfinite(coefs.grad)))
    uni = bumps([0.45, 0.55])
    ures = register(uni, landmarks=[0.45, 0.55])
    tensor = ures.apply(FData(torch.tensor(np.asarray(uni.coefs)), uni.basis))
    assert tuple(tensor.coefs.shape) == uni.coefs.shape
