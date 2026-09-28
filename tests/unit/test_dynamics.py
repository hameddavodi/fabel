"""Unit tests for :mod:`fabel.dynamics`: principal differential analysis."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any, cast

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.base import clone
from sklearn.exceptions import NotFittedError

from fabel import LDO, BSpline, Constant, FData, Fourier, Monomial
from fabel.dynamics import PDA, PDAStability, phase_plane
from fabel.smoothing import smooth

TWO_PI = 2.0 * np.pi


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


def harmonic_curves(n_curves: int = 4, seed: int = 0) -> FData:
    """Curves ``a sin t + b cos t`` on one period: all solve ``D²x + x = 0``."""
    rng = np.random.default_rng(seed)
    basis = Fourier(domain=(0.0, TWO_PI), n_basis=3)
    coefs = np.zeros((3, n_curves))
    coefs[1:, :] = rng.normal(size=(2, n_curves))
    return FData(coefs, basis)


def oscillator_system(n_curves: int = 3, seed: int = 1) -> FData:
    """Two variables ``x1 = r sin(t + p)``, ``x2 = r cos(t + p)``: ``Dx1 = x2``, ``Dx2 = -x1``."""
    rng = np.random.default_rng(seed)
    basis = Fourier(domain=(0.0, TWO_PI), n_basis=3)
    coefs = np.zeros((3, n_curves, 2))
    scale = np.sqrt(np.pi)  # Fourier basis functions are normalised: sin(t) / sqrt(pi)
    for n in range(n_curves):
        r, p = rng.uniform(0.5, 2.0), rng.uniform(0.0, TWO_PI)
        # sin(t + p) = cos p sin t + sin p cos t ; basis order is (1, sin, cos)
        coefs[1, n, 0], coefs[2, n, 0] = r * np.cos(p) * scale, r * np.sin(p) * scale
        # cos(t + p) = cos p cos t - sin p sin t
        coefs[1, n, 1], coefs[2, n, 1] = -r * np.sin(p) * scale, r * np.cos(p) * scale
    return FData(coefs, basis)


def forced_curves(n_curves: int = 4, seed: int = 3) -> tuple[FData, FData]:
    """Curves ``x = a sin t + b cos t`` and inputs ``u`` with ``Dx = -x + u`` exactly.

    ``Dx + x = (a - b) sin t + (a + b) cos t``, which is ``u`` for each curve.
    """
    rng = np.random.default_rng(seed)
    basis = Fourier(domain=(0.0, TWO_PI), n_basis=3)
    a, b = rng.normal(size=(2, n_curves))
    x = np.zeros((3, n_curves))
    u = np.zeros((3, n_curves))
    x[1], x[2] = a, b
    u[1], u[2] = a - b, a + b
    return FData(x, basis), FData(u, basis)


def quadratic_curve() -> FData:
    """The single curve ``x(t) = 1 + t²`` on ``[0, 1]`` in a monomial basis."""
    return FData(np.array([[1.0], [0.0], [1.0]]), Monomial(domain=(0.0, 1.0), n_basis=3))


# --------------------------------------------------------------------------- #
# estimation
# --------------------------------------------------------------------------- #


def test_harmonic_curves_recover_the_harmonic_equation() -> None:
    pda = PDA(order=2, n_grid=None).fit(harmonic_curves())
    beta0, beta1 = (float(w.coefs[0, 0]) for w in pda.weights_)
    assert beta0 == pytest.approx(1.0, abs=1e-12)
    assert beta1 == pytest.approx(0.0, abs=1e-12)
    grid = np.linspace(0.0, TWO_PI, 9)
    np.testing.assert_allclose(pda.residuals_(grid), 0.0, atol=1e-10)


def test_the_default_weight_basis_is_constant_on_the_domain() -> None:
    pda = PDA(order=2).fit(harmonic_curves())
    assert len(pda.weights_) == 2
    for weight in pda.weights_:
        assert isinstance(weight.basis, Constant)
        assert weight.domain == (0.0, TWO_PI)
        assert weight.n_curves == 1


def test_exact_quadrature_gives_the_closed_form() -> None:
    # For x = 1 + t² with a constant weight, β = -∫x Dx / ∫x² = -(3/2) / (28/15).
    pda = PDA(order=1, n_grid=None).fit(quadratic_curve())
    assert float(pda.weights_[0].coefs[0, 0]) == pytest.approx(-45.0 / 56.0, rel=1e-13)


def test_trapezoid_grid_converges_to_the_exact_answer() -> None:
    exact = -45.0 / 56.0
    errors = [
        abs(float(PDA(order=1, n_grid=n).fit(quadratic_curve()).weights_[0].coefs[0, 0]) - exact)
        for n in (11, 101, 1001)
    ]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 1e-6


def test_the_operator_applies_the_fitted_equation() -> None:
    curves = harmonic_curves()
    pda = PDA(order=2, n_grid=None).fit(curves)
    assert isinstance(pda.operator_, LDO)
    assert pda.operator_.order == 2
    grid = np.linspace(0.0, TWO_PI, 7)
    np.testing.assert_allclose(curves(grid, pda.operator_), 0.0, atol=1e-10)


def test_functional_weights_with_a_roughness_penalty() -> None:
    curves = harmonic_curves(6)
    weight_basis = BSpline(domain=(0.0, TWO_PI), n_basis=7, order=4)
    rough = [
        PDA(order=2, weight_basis=weight_basis, lam=lam, n_grid=None).fit(curves)
        for lam in (0.0, 1e3)
    ]
    grid = np.linspace(0.0, TWO_PI, 11)
    # the harmonic equation lies in the weight space, so either way β0 = 1, β1 = 0
    for pda in rough:
        np.testing.assert_allclose(pda.weights_[0](grid)[:, 0], 1.0, atol=1e-8)
        np.testing.assert_allclose(pda.weights_[1](grid)[:, 0], 0.0, atol=1e-8)


def test_the_penalty_smooths_the_weights() -> None:
    rng = np.random.default_rng(5)
    basis = BSpline(domain=(0.0, 1.0), n_basis=12, order=6)
    curves = FData(rng.normal(size=(12, 8)), basis)
    weight_basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=4)
    penalty = np.asarray(weight_basis.penalty(2))

    def roughness(lam: float) -> float:
        pda = PDA(order=2, weight_basis=weight_basis, lam=lam, n_grid=None).fit(curves)
        return sum(float(w.coefs[:, 0] @ penalty @ w.coefs[:, 0]) for w in pda.weights_)

    assert roughness(0.0) > roughness(1e-2) > roughness(1e2)


def test_per_derivative_weight_bases_and_lambdas() -> None:
    curves = harmonic_curves(5)
    bases = [Constant(domain=(0.0, TWO_PI)), BSpline(domain=(0.0, TWO_PI), n_basis=5)]
    pda = PDA(order=2, weight_basis=bases, lam=[0.0, 1.0], penalty=LDO(2), n_grid=None).fit(curves)
    assert isinstance(pda.weights_[0].basis, Constant)
    assert isinstance(pda.weights_[1].basis, BSpline)
    assert float(pda.weights_[0].coefs[0, 0]) == pytest.approx(1.0, abs=1e-10)


def test_system_recovers_the_oscillator() -> None:
    curves = oscillator_system()
    pda = PDA(order=1, n_grid=None).fit(curves)
    # D x1 + β11 x1 + β12 x2 = 0 with D x1 = x2  -> β11 = 0, β12 = -1
    # D x2 + β21 x1 + β22 x2 = 0 with D x2 = -x1 -> β21 = 1, β22 = 0
    expected = [[0.0, -1.0], [1.0, 0.0]]
    for i in range(2):
        for k in range(2):
            assert float(pda.weights_[i][k][0].coefs[0, 0]) == pytest.approx(
                expected[i][k], abs=1e-10
            )
    assert pda.operator_ is None
    assert pda.residuals_.n_vars == 2
    np.testing.assert_allclose(pda.residuals_(np.linspace(0, TWO_PI, 5)), 0.0, atol=1e-10)


def test_system_trapezoid_matches_the_exact_fit_on_an_exact_problem() -> None:
    curves = oscillator_system()
    pda = PDA(order=1).fit(curves)
    assert float(pda.weights_[0][1][0].coefs[0, 0]) == pytest.approx(-1.0, abs=1e-10)


@settings(max_examples=25, deadline=None)
@given(
    scale=st.floats(min_value=1e-3, max_value=1e3),
    seed=st.integers(min_value=0, max_value=10_000),
)
def test_weights_are_invariant_to_scaling_the_curves(scale: float, seed: int) -> None:
    rng = np.random.default_rng(seed)
    basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=5)
    curves = FData(rng.normal(size=(9, 5)), basis)
    scaled = FData(scale * np.asarray(curves.coefs), basis)
    a = PDA(order=2).fit(curves)
    b = PDA(order=2).fit(scaled)
    for wa, wb in zip(a.weights_, b.weights_, strict=True):
        np.testing.assert_allclose(wb.coefs, wa.coefs, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(
        b.residuals_.coefs, scale * np.asarray(a.residuals_.coefs), rtol=1e-7, atol=1e-9 * scale
    )


@settings(max_examples=15, deadline=None)
@given(seed=st.integers(min_value=0, max_value=10_000))
def test_weights_do_not_depend_on_curve_order(seed: int) -> None:
    rng = np.random.default_rng(seed)
    basis = BSpline(domain=(0.0, 1.0), n_basis=8, order=5)
    coefs = rng.normal(size=(8, 6))
    order = rng.permutation(6)
    a = PDA(order=1, n_grid=None).fit(FData(coefs, basis))
    b = PDA(order=1, n_grid=None).fit(FData(coefs[:, order], basis))
    np.testing.assert_allclose(b.weights_[0].coefs, a.weights_[0].coefs, rtol=1e-10)


# --------------------------------------------------------------------------- #
# transform and solve
# --------------------------------------------------------------------------- #


def test_transform_applies_the_fitted_equation_to_new_curves() -> None:
    pda = PDA(order=2, n_grid=None).fit(harmonic_curves(4, seed=0))
    fresh = harmonic_curves(3, seed=9)
    residuals = pda.transform(fresh)
    assert residuals.n_curves == 3
    np.testing.assert_allclose(residuals(np.linspace(0.0, TWO_PI, 7)), 0.0, atol=1e-10)


def test_fit_transform_returns_the_residuals() -> None:
    curves = harmonic_curves()
    residuals = PDA(order=2).fit_transform(curves)
    np.testing.assert_allclose(residuals.coefs, PDA(order=2).fit(curves).residuals_.coefs)


def test_solve_integrates_the_harmonic_equation() -> None:
    pda = PDA(order=2, n_grid=None).fit(harmonic_curves())
    grid = np.linspace(0.0, TWO_PI, 50)
    np.testing.assert_allclose(pda.solve(grid, [0.0, 1.0]), np.sin(grid), atol=1e-8)


def test_solve_runs_backwards_from_the_first_point() -> None:
    pda = PDA(order=2, n_grid=None).fit(harmonic_curves())
    grid = np.linspace(TWO_PI, 0.0, 20)
    np.testing.assert_allclose(pda.solve(grid, [0.0, 1.0]), np.sin(grid - TWO_PI), atol=1e-8)


def test_solve_a_system() -> None:
    pda = PDA(order=1, n_grid=None).fit(oscillator_system())
    grid = np.linspace(0.0, TWO_PI, 25)
    solution = pda.solve(grid, [[0.0], [1.0]])
    assert solution.shape == (25, 2)
    np.testing.assert_allclose(solution[:, 0], np.sin(grid), atol=1e-8)
    np.testing.assert_allclose(solution[:, 1], np.cos(grid), atol=1e-8)


def test_solve_a_single_point_returns_the_initial_value() -> None:
    pda = PDA(order=2, n_grid=None).fit(harmonic_curves())
    np.testing.assert_allclose(pda.solve([1.0], [0.25, 3.0]), [0.25])


# --------------------------------------------------------------------------- #
# forcing functions
# --------------------------------------------------------------------------- #


def test_forcing_weights_recover_the_forced_equation() -> None:
    x, u = forced_curves()
    pda = PDA(order=1, n_grid=None).fit(x, forcing=u)
    assert float(pda.weights_[0].coefs[0, 0]) == pytest.approx(1.0, abs=1e-12)
    assert len(pda.forcing_weights_) == 1
    assert float(pda.forcing_weights_[0].coefs[0, 0]) == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(pda.residuals_(np.linspace(0.0, TWO_PI, 9)), 0.0, atol=1e-10)


def test_unforced_fit_has_no_forcing_weights() -> None:
    assert PDA(order=2).fit(harmonic_curves()).forcing_weights_ == ()
    assert PDA(order=1).fit(oscillator_system()).forcing_weights_ == ()


def test_forcing_weight_scales_inversely_with_the_input() -> None:
    x, u = forced_curves()
    half = FData(0.5 * np.asarray(u.coefs), u.basis)
    pda = PDA(order=1, n_grid=None).fit(x, forcing=[half])
    assert float(pda.forcing_weights_[0].coefs[0, 0]) == pytest.approx(2.0, abs=1e-11)


@settings(max_examples=20, deadline=None)
@given(
    scale=st.floats(min_value=1e-2, max_value=1e2),
    seed=st.integers(min_value=0, max_value=10_000),
)
def test_forced_weights_are_invariant_to_scaling_curves_and_inputs(scale: float, seed: int) -> None:
    rng = np.random.default_rng(seed)
    basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=5)
    curves = FData(rng.normal(size=(9, 4)), basis)
    inputs = FData(rng.normal(size=(9, 4)), basis)
    a = PDA(order=1).fit(curves, forcing=inputs)
    b = PDA(order=1).fit(
        FData(scale * np.asarray(curves.coefs), basis),
        forcing=FData(scale * np.asarray(inputs.coefs), basis),
    )
    np.testing.assert_allclose(b.weights_[0].coefs, a.weights_[0].coefs, rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(
        b.forcing_weights_[0].coefs, a.forcing_weights_[0].coefs, rtol=1e-7, atol=1e-10
    )


def test_a_single_forcing_curve_is_shared_by_all_curves() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=5)
    rng = np.random.default_rng(4)
    curves = FData(rng.normal(size=(9, 3)), basis)
    one = FData(rng.normal(size=(9, 1)), basis)
    three = FData(np.repeat(np.asarray(one.coefs), 3, axis=1), basis)
    shared = PDA(order=1).fit(curves, forcing=one)
    repeated = PDA(order=1).fit(curves, forcing=three)
    np.testing.assert_allclose(
        shared.forcing_weights_[0].coefs, repeated.forcing_weights_[0].coefs, rtol=1e-12
    )
    np.testing.assert_allclose(shared.residuals_.coefs, repeated.residuals_.coefs, atol=1e-12)


def test_several_forcing_functions_with_their_own_bases_and_penalties() -> None:
    x, u = forced_curves(6)
    rng = np.random.default_rng(7)
    extra = FData(rng.normal(size=(3, 6)), u.basis)
    bspline = BSpline(domain=(0.0, TWO_PI), n_basis=6)
    pda = PDA(
        order=1,
        n_grid=None,
        forcing_basis=[Constant(domain=(0.0, TWO_PI)), bspline],
        forcing_lam=[0.0, 1e-2],
    ).fit(x, forcing=[u, extra])
    first, second = pda.forcing_weights_
    assert isinstance(first.basis, Constant)
    assert second.basis is bspline
    assert float(first.coefs[0, 0]) == pytest.approx(1.0, abs=1e-8)
    np.testing.assert_allclose(second(np.linspace(0.0, TWO_PI, 7)), 0.0, atol=1e-8)


def test_the_forcing_penalty_smooths_the_forcing_weight() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=9, order=5)
    rng = np.random.default_rng(5)
    curves = FData(rng.normal(size=(9, 4)), basis)
    inputs = FData(rng.normal(size=(9, 4)), basis)
    wbasis = BSpline(domain=(0.0, 1.0), n_basis=7)
    rough = PDA(order=1, forcing_basis=wbasis).fit(curves, forcing=inputs)
    smoothed = PDA(order=1, forcing_basis=wbasis, forcing_lam=10.0).fit(curves, forcing=inputs)
    grid = np.linspace(0.0, 1.0, 201)
    rough_curvature = np.sum(rough.forcing_weights_[0](grid, 2) ** 2)
    smooth_curvature = np.sum(smoothed.forcing_weights_[0](grid, 2) ** 2)
    assert smooth_curvature < rough_curvature


def test_a_forced_system_with_an_unforced_equation() -> None:
    # x1 + phi0 (phi0 the constant Fourier function) and x2 solve
    # D x1 = x2 and D x2 = -(x1 + phi0) + phi0: equation 2 is forced by phi0.
    curves = oscillator_system()
    coefs = np.asarray(curves.coefs).copy()
    coefs[0, :, 0] = 1.0
    phi0 = FData(np.array([[1.0], [0.0], [0.0]]), curves.basis)
    pda = PDA(order=1, n_grid=None).fit(FData(coefs, curves.basis), forcing=[None, phi0])
    assert pda.forcing_weights_[0] == ()
    (alpha,) = pda.forcing_weights_[1]
    assert float(alpha.coefs[0, 0]) == pytest.approx(1.0, abs=1e-9)
    assert float(pda.weights_[1][0][0].coefs[0, 0]) == pytest.approx(1.0, abs=1e-9)
    assert float(pda.weights_[0][1][0].coefs[0, 0]) == pytest.approx(-1.0, abs=1e-9)


def test_transform_and_fit_transform_with_forcing() -> None:
    x, u = forced_curves()
    pda = PDA(order=1, n_grid=None).fit(x, forcing=u)
    fresh_x, fresh_u = forced_curves(2, seed=11)
    residuals = pda.transform(fresh_x, forcing=fresh_u)
    np.testing.assert_allclose(residuals(np.linspace(0.0, TWO_PI, 7)), 0.0, atol=1e-10)
    both = PDA(order=1).fit_transform(x, forcing=u)
    np.testing.assert_allclose(both.coefs, PDA(order=1).fit(x, forcing=u).residuals_.coefs)


def test_transform_needs_matching_forcing() -> None:
    x, u = forced_curves()
    forced = PDA(order=1).fit(x, forcing=u)
    with pytest.raises(ValueError, match="forcing functions per equation"):
        forced.transform(x)
    with pytest.raises(ValueError, match="forcing functions per equation"):
        PDA(order=1).fit(x).transform(x, forcing=u)


def decaying_to_a_level() -> tuple[PDA, FData]:
    """PDA of solutions of ``x'' + 1.5 x' + 2 x = 3`` (level 1.5) forced by ``u = 1``."""
    domain = (0.0, 30.0)
    t = np.linspace(0.0, 30.0, 601)
    omega = np.sqrt(2.0 - 0.75**2)
    decay = np.exp(-0.75 * t)
    y = np.column_stack(
        [
            1.5 + decay * np.cos(omega * t),
            1.5 + decay * np.sin(omega * t),
            1.5 - decay * (np.cos(omega * t) + 0.5 * np.sin(omega * t)),
        ]
    )
    basis = BSpline(domain=domain, n_basis=120, order=6)
    curves = smooth(y, t, basis=basis, lam=0.0).fd
    one = FData(np.array([1.0]), Constant(domain=domain))
    return PDA(order=2, n_grid=1001).fit(curves, forcing=one), one


def test_solve_with_forcing_follows_the_inhomogeneous_equation() -> None:
    pda, one = decaying_to_a_level()
    beta0, beta1 = (float(w.coefs[0, 0]) for w in pda.weights_)
    alpha = float(pda.forcing_weights_[0].coefs[0, 0])
    assert (beta0, beta1, alpha) == pytest.approx((2.0, 1.5, 3.0), rel=1e-4)
    grid = np.linspace(0.0, 30.0, 7)
    level = alpha / beta0
    forced = pda.solve(grid, [level, 0.0], forcing=one)
    np.testing.assert_allclose(forced, level, rtol=1e-8)
    free = pda.solve(grid, [level, 0.0])
    assert abs(free[-1]) < 1e-6
    with pytest.raises(ValueError, match="1 curves"):
        pda.solve(grid, [0.0, 0.0], forcing=FData(np.ones((1, 2)), one.basis))


@pytest.mark.parametrize(
    ("forcing", "error", "match"),
    [
        ("not data", TypeError, "FData"),
        ([np.ones(3)], TypeError, "FData"),
        ("three curves", ValueError, "curves"),
        ("two variables", ValueError, "one variable"),
        ("other domain", ValueError, "domain"),
    ],
)
def test_invalid_forcing(forcing: Any, error: type[Exception], match: str) -> None:
    x, u = forced_curves()
    value = (
        {
            "three curves": FData(np.ones((3, 3)), u.basis),
            "two variables": FData(np.ones((3, 4, 2)), u.basis),
            "other domain": FData(np.ones((3, 1)), Fourier(domain=(0.0, 1.0), n_basis=3)),
        }.get(forcing, forcing)
        if isinstance(forcing, str)
        else forcing
    )
    with pytest.raises(error, match=match):
        PDA(order=1).fit(x, forcing=value)


def test_invalid_system_forcing() -> None:
    curves = oscillator_system()
    u = FData(np.ones((3, 1)), curves.basis)
    with pytest.raises(ValueError, match="one forcing entry per equation"):
        PDA(order=1).fit(curves, forcing=u)
    with pytest.raises(ValueError, match="2 entries"):
        PDA(order=1).fit(curves, forcing=[u])


@pytest.mark.parametrize(
    ("params", "error", "match"),
    [
        ({"forcing_lam": -1.0}, ValueError, "forcing_lam"),
        ({"forcing_lam": [0.0, 0.0]}, ValueError, "forcing_lam"),
        ({"forcing_lam": ["x"]}, TypeError, "forcing_lam"),
        ({"forcing_basis": Constant(domain=(0.0, 1.0))}, ValueError, "domain"),
        ({"forcing_basis": [None, None]}, ValueError, "forcing_basis"),
    ],
)
def test_invalid_forcing_parameters(
    params: dict[str, Any], error: type[Exception], match: str
) -> None:
    x, u = forced_curves()
    with pytest.raises(error, match=match):
        PDA(order=1, **params).fit(x, forcing=u)


def test_system_forcing_parameters_per_equation() -> None:
    curves = oscillator_system()
    u = FData(np.ones((3, 1)), curves.basis)
    bspline = BSpline(domain=(0.0, TWO_PI), n_basis=5)
    pda = PDA(order=1, n_grid=None, forcing_basis=[None, [bspline]], forcing_lam=[0.0, [1.0]]).fit(
        curves, forcing=[u, [u]]
    )
    assert isinstance(pda.forcing_weights_[0][0].basis, Constant)
    assert pda.forcing_weights_[1][0].basis is bspline
    with pytest.raises(ValueError, match="one per equation"):
        PDA(order=1, forcing_basis=[None]).fit(curves, forcing=[u, u])


# --------------------------------------------------------------------------- #
# stability
# --------------------------------------------------------------------------- #


def test_stability_of_the_harmonic_equation() -> None:
    result = PDA(order=2, n_grid=None).fit(harmonic_curves()).stability(n_points=5)
    assert isinstance(result, PDAStability)
    np.testing.assert_allclose(result.t, np.linspace(0.0, TWO_PI, 5))
    np.testing.assert_allclose(result.eigenvalues, np.tile([1j, -1j], (5, 1)), atol=1e-10)
    np.testing.assert_array_equal(result.limits, 0.0)


def test_stability_of_a_system_at_given_times() -> None:
    result = PDA(order=1, n_grid=None).fit(oscillator_system()).stability([0.0, 1.0, 2.0])
    assert result.eigenvalues.shape == (3, 2)
    np.testing.assert_allclose(np.abs(result.eigenvalues), 1.0, atol=1e-10)


def test_eigenvalues_are_ordered_by_decreasing_modulus() -> None:
    t = np.linspace(0.0, 1.0, 201)
    basis = BSpline(domain=(0.0, 1.0), n_basis=30, order=6)
    # x'' + 3x' + 2x = 0 has the roots -1 and -2
    y = np.column_stack([np.exp(-t), np.exp(-2 * t), np.exp(-t) - 0.5 * np.exp(-2 * t)])
    pda = PDA(order=2).fit(smooth(y, t, basis=basis, lam=0.0).fd)
    values = pda.stability(n_points=3).eigenvalues
    np.testing.assert_allclose(values, np.tile([-2.0, -1.0], (3, 1)), atol=1e-4)


def test_the_equilibrium_is_where_the_forced_solution_settles() -> None:
    pda, one = decaying_to_a_level()
    beta0 = float(pda.weights_[0].coefs[0, 0])
    alpha = float(pda.forcing_weights_[0].coefs[0, 0])
    limits = pda.stability(n_points=4).limits
    np.testing.assert_allclose(limits, np.tile([alpha / beta0, 0.0], (4, 1)), atol=1e-12)
    settled = pda.solve(np.linspace(0.0, 30.0, 3), [0.0, 0.0], forcing=one)[-1]
    assert settled == pytest.approx(limits[0, 0], rel=1e-8)


def test_limits_are_nan_where_the_system_matrix_is_singular() -> None:
    x, u = forced_curves()
    pda = PDA(order=1).fit(x, forcing=FData(np.ones((3, 1)), u.basis))
    pda._nested = [[[FData(np.array([0.0]), Constant(domain=pda.domain_))]]]
    assert np.all(np.isnan(pda.stability(n_points=3).limits))


def test_stability_uses_the_single_curve_fitting_forcing_by_default() -> None:
    x, u = forced_curves()
    one = FData(np.ones((3, 1)), u.basis)
    pda = PDA(order=1).fit(x, forcing=one)
    np.testing.assert_allclose(
        pda.stability(n_points=4).limits, pda.stability(n_points=4, forcing=one).limits
    )
    several = PDA(order=1).fit(x, forcing=u)
    with pytest.raises(ValueError, match="one curve per forcing function"):
        several.stability()
    with pytest.raises(ValueError, match="forcing functions per equation"):
        PDA(order=1).fit(x).stability(forcing=one)
    with pytest.raises(ValueError, match="n_points"):
        pda.stability(n_points=0)
    with pytest.raises(ValueError, match="domain"):
        pda.stability([100.0])


def test_stability_plot_draws_real_and_imaginary_parts() -> None:
    result = PDA(order=2).fit(harmonic_curves()).stability(n_points=11)
    axes = result.plot(color="k")
    assert len(axes.lines) == 5
    _, given_ax = plt.subplots()
    assert result.plot(ax=given_ax) is given_ax


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #


def test_fit_rejects_non_functional_data() -> None:
    with pytest.raises(TypeError, match="FData"):
        PDA().fit(cast(Any, np.ones((5, 3))))


@pytest.mark.parametrize(
    ("params", "match"),
    [
        ({"order": 0}, "order"),
        ({"n_grid": 1}, "n_grid"),
        ({"lam": -1.0}, "lam"),
        ({"lam": [0.0]}, "lam"),
        ({"weight_basis": [Constant(domain=(0.0, TWO_PI))]}, "weight_basis"),
        ({"weight_basis": Constant(domain=(0.0, 1.0))}, "domain"),
    ],
)
def test_invalid_parameters(params: dict[str, Any], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        PDA(**params).fit(harmonic_curves())


def test_grid_must_resolve_the_curve_basis() -> None:
    curves = FData(np.eye(12)[:, :3], BSpline(domain=(0.0, 1.0), n_basis=12))
    with pytest.raises(ValueError, match="n_grid"):
        PDA(order=1, n_grid=5).fit(curves)


def test_singular_problem_is_reported() -> None:
    zero = FData(np.zeros((3, 2)), Fourier(domain=(0.0, TWO_PI), n_basis=3))
    with pytest.raises(ValueError, match="singular"):
        PDA(order=1).fit(zero)


def test_linearly_dependent_derivatives_are_reported() -> None:
    # In the oscillator D x1 = x2, so (x1, x2, D x1, D x2) cannot all carry a weight.
    with pytest.raises(ValueError, match="linearly dependent"):
        PDA(order=2).fit(oscillator_system())


def test_unfitted_estimator_raises() -> None:
    with pytest.raises(NotFittedError):
        PDA().transform(harmonic_curves())
    with pytest.raises(NotFittedError):
        PDA().solve([0.0, 1.0], [1.0, 0.0])


def test_transform_checks_variables_and_domain() -> None:
    pda = PDA(order=1).fit(oscillator_system())
    with pytest.raises(ValueError, match="variables"):
        pda.transform(harmonic_curves())
    other = FData(np.ones((3, 1)), Fourier(domain=(0.0, 1.0), n_basis=3))
    with pytest.raises(ValueError, match="domain"):
        PDA(order=1).fit(harmonic_curves()).transform(other)
    with pytest.raises(TypeError, match="FData"):
        pda.transform(cast(Any, np.ones(3)))


def test_solve_validates_its_arguments() -> None:
    pda = PDA(order=2).fit(harmonic_curves())
    with pytest.raises(ValueError, match="initial"):
        pda.solve([0.0, 1.0], [1.0])
    with pytest.raises(ValueError, match="domain"):
        pda.solve([0.0, 100.0], [1.0, 0.0])
    with pytest.raises(ValueError, match="monotone"):
        pda.solve([0.0, 2.0, 1.0], [1.0, 0.0])
    with pytest.raises(ValueError, match="at least one"):
        pda.solve([], [1.0, 0.0])


def test_estimator_is_clonable() -> None:
    pda = PDA(order=3, lam=0.5, n_grid=None)
    twin = clone(pda)
    assert twin.get_params() == pda.get_params()


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #


def test_plot_overlay_draws_the_weight_trajectory() -> None:
    curves = harmonic_curves(5)
    weight_basis = BSpline(domain=(0.0, TWO_PI), n_basis=5)
    pda = PDA(order=2, weight_basis=weight_basis, lam=1.0).fit(curves)
    ax = pda.plot_overlay(n_points=31)
    trajectory = np.asarray(ax.lines[0].get_xydata())
    grid = np.linspace(0.0, TWO_PI, 31)
    np.testing.assert_allclose(trajectory[:, 0], pda.weights_[1](grid)[:, 0])
    np.testing.assert_allclose(trajectory[:, 1], pda.weights_[0](grid)[:, 0])
    boundary = np.asarray(ax.lines[1].get_xydata())
    np.testing.assert_allclose(boundary[:, 1], boundary[:, 0] ** 2 / 4.0)
    assert ax.get_xlabel()
    assert ax.get_ylabel()


def test_plot_overlay_marks_time_labels() -> None:
    pda = PDA(order=2).fit(harmonic_curves())
    _, ax = plt.subplots()
    out = pda.plot_overlay(ax=ax, labels={1.0: "a", 2.0: "b"})
    assert out is ax
    assert [text.get_text() for text in ax.texts] == ["a", "b"]


def test_plot_overlay_needs_a_single_second_order_equation() -> None:
    with pytest.raises(ValueError, match="second-order"):
        PDA(order=1).fit(harmonic_curves()).plot_overlay()
    rng = np.random.default_rng(3)
    system = FData(rng.normal(size=(8, 6, 2)), BSpline(domain=(0.0, 1.0), n_basis=8, order=5))
    with pytest.raises(ValueError, match="second-order"):
        PDA(order=2).fit(system).plot_overlay()


def test_phase_plane_plots_velocity_against_acceleration() -> None:
    curves = harmonic_curves(2)
    ax = phase_plane(curves)
    assert len(ax.lines) == 2
    xy = np.asarray(ax.lines[0].get_xydata())
    grid = np.linspace(0.0, TWO_PI, xy.shape[0])
    np.testing.assert_allclose(xy[:, 0], curves(grid, 1)[:, 0])
    np.testing.assert_allclose(xy[:, 1], curves(grid, 2)[:, 0])


def test_phase_plane_with_custom_grid_derivatives_and_labels() -> None:
    curves = harmonic_curves(1)
    t = np.linspace(0.0, np.pi, 11)
    _, ax = plt.subplots()
    out = phase_plane(curves, t, deriv=(0, 1), labels={0.5: "m"}, ax=ax, color="k")
    assert out is ax
    np.testing.assert_allclose(np.asarray(ax.lines[0].get_xydata())[:, 0], curves(t)[:, 0])
    assert [text.get_text() for text in ax.texts] == ["m"]
    assert ax.get_xlabel() == "D0 x"
    assert ax.get_ylabel() == "D1 x"


def test_phase_plane_rejects_multivariate_curves() -> None:
    with pytest.raises(ValueError, match="variable"):
        phase_plane(oscillator_system())
