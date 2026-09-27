"""Unit tests for :mod:`fabel.registration`."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.exceptions import NotFittedError

from fabel import Basis, BSpline, FData, Fourier, Monomial
from fabel.registration import (
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
    from fabel._backend import default_namespace

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
    from fabel.datasets import load_growth
    from fabel.smoothing import smooth

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
    with pytest.raises(ValueError, match="univariate"):
        register(wide)
    with pytest.raises(ValueError, match="univariate"):
        register(fd, wide)
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
