"""Unit tests for :mod:`fdatools.density` (penalised density and intensity fits)."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fdatools import LDO, BSpline, FData, Fourier, Monomial
from fdatools._linalg import composite_gauss_legendre
from fdatools.density import (
    DensityResult,
    IntensityResult,
    _null_direction,
    _Problem,
    fit_density,
    fit_intensity,
)

DOMAIN = (0.0, 1.0)


def _integral(
    values_of: Callable[[np.ndarray], Any],
    domain: tuple[float, float] = DOMAIN,
    panels: int = 120,
) -> float:
    nodes, weights = composite_gauss_legendre(np.linspace(domain[0], domain[1], panels + 1), 12)
    return float(weights @ np.asarray(values_of(nodes)))


@pytest.fixture
def sample() -> np.ndarray:
    rng = np.random.default_rng(7)
    return rng.beta(2.0, 5.0, size=150)


# --------------------------------------------------------------------------- #
# derivatives of the criteria
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", ["density", "intensity"])
def test_gradient_and_hessian_match_finite_differences(kind: str, sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=7)
    problem = _Problem(kind, basis, sample, basis.penalty(2), 0.3)
    coefs = np.linspace(-0.5, 0.7, 7)
    gradient, hessian = problem.derivatives(coefs)
    step = 1e-6
    for k in range(7):
        e = np.zeros(7)
        e[k] = step
        numeric = (problem.value(coefs + e) - problem.value(coefs - e)) / (2 * step)
        assert numeric == pytest.approx(gradient[k], rel=1e-6, abs=1e-6)
        g_plus, _ = problem.derivatives(coefs + e)
        g_minus, _ = problem.derivatives(coefs - e)
        np.testing.assert_allclose((g_plus - g_minus) / (2 * step), hessian[:, k], atol=1e-5)


def test_intensity_value_is_infinite_on_overflow() -> None:
    basis = BSpline(domain=DOMAIN, n_basis=5)
    problem = _Problem("intensity", basis, np.array([0.5]), basis.penalty(2), 0.0)
    assert problem.value(np.full(5, 800.0)) == float("inf")


# --------------------------------------------------------------------------- #
# density
# --------------------------------------------------------------------------- #


def test_density_integrates_to_one_and_is_centred(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=9)
    res = fit_density(sample, basis=basis, lam=1e-3)
    assert isinstance(res, DensityResult)
    assert res.converged
    assert res.gradient_norm < 1e-9
    assert _integral(res) == pytest.approx(1.0, abs=1e-12)
    # W is only fixed up to a constant here: the returned member has int W = 0.
    assert _integral(lambda t: res.fd(t)[:, 0]) == pytest.approx(0.0, abs=1e-12)
    assert res.normaliser == pytest.approx(1.0 / _integral(lambda t: np.exp(res.fd(t)[:, 0])))
    assert res.log_likelihood == pytest.approx(float(np.sum(res.log_density(sample))))
    np.testing.assert_allclose(np.exp(res.log_density(sample)), res(sample))
    assert res.penalty_matrix.shape == (9, 9)
    assert res.lam == 1e-3


def test_density_is_normalised_intensity(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=8)
    dens = fit_density(sample, basis=basis, lam=0.05, penalty=1)
    inten = fit_intensity(sample, basis=basis, lam=0.05, penalty=1)
    assert inten.expected_count == pytest.approx(sample.size, rel=1e-12)
    t = np.linspace(0.0, 1.0, 17)
    np.testing.assert_allclose(dens(t), inten(t) / sample.size, rtol=1e-9)


def test_density_criterion_is_minimised(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=8)
    res = fit_density(sample, basis=basis, lam=0.1)
    problem = _Problem("density", basis, sample, basis.penalty(2), 0.1)
    coefs = res.fd.coefs[:, 0]
    assert problem.value(coefs) == pytest.approx(res.criterion)
    rng = np.random.default_rng(0)
    for _ in range(5):
        assert problem.value(coefs + 1e-3 * rng.normal(size=8)) > res.criterion


def test_density_without_constants_is_identified() -> None:
    # No constant function in the span: W is unique and is not shifted.
    rng = np.random.default_rng(3)
    x = np.clip(rng.normal(0.4, 0.3, size=200), -1.0, 1.0)
    basis = Monomial(domain=(-1.0, 1.0), exponents=[1, 2])
    res = fit_density(x, basis=basis)
    problem = _Problem("density", basis, x, basis.penalty(2), 0.0)
    assert _null_direction(problem) is None
    assert res.converged
    assert _integral(res, (-1.0, 1.0)) == pytest.approx(1.0, abs=1e-12)
    # MLE of a truncated normal: the mean of the sample matches the model mean.
    model_mean = _integral(lambda t: t * res(t), (-1.0, 1.0))
    assert model_mean == pytest.approx(float(np.mean(x)), abs=1e-10)


def test_density_penalty_charging_constants_is_identified(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=6)
    operator = LDO(weights=[1.0, 0.0])
    problem = _Problem("density", basis, sample, basis.penalty(operator), 1.0)
    assert _null_direction(problem) is None
    res = fit_density(sample, basis=basis, lam=1.0, penalty=operator)
    assert res.converged
    assert _integral(res) == pytest.approx(1.0, abs=1e-12)


def test_density_zero_lambda_keeps_null_direction(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=6)
    problem = _Problem("density", basis, sample, basis.penalty(LDO(weights=[1.0, 0.0])), 0.0)
    assert _null_direction(problem) is not None


def test_heavy_third_derivative_penalty_gives_normal_density() -> None:
    rng = np.random.default_rng(11)
    x = rng.normal(size=400)
    x = x[np.abs(x) < 4]
    basis = BSpline(domain=(-4.0, 4.0), n_basis=12)
    res = fit_density(x, basis=basis, lam=1e8, penalty=LDO(3))
    # log p is (numerically) quadratic: a normal with the sample's moments.
    t = np.linspace(-3.0, 3.0, 13)
    second = np.diff(res.log_density(t), 2)
    np.testing.assert_allclose(second, second[0], atol=1e-4)


def test_heavy_first_derivative_penalty_gives_uniform(sample: np.ndarray) -> None:
    res = fit_density(sample, basis=BSpline(domain=DOMAIN, n_basis=6), lam=1e10, penalty=1)
    np.testing.assert_allclose(res(np.linspace(0, 1, 5)), 1.0, atol=1e-6)


def test_density_default_basis_and_domain() -> None:
    rng = np.random.default_rng(5)
    x = rng.gamma(3.0, size=250)
    res = fit_density(x, lam=1e-2)
    assert res.fd.basis.domain == (float(x.min()), float(x.max()))
    assert res.fd.basis.n_basis == 20
    wide = fit_density(x, domain=(0.0, 20.0), lam=1e-2)
    assert wide.fd.basis.domain == (0.0, 20.0)
    small = fit_density(x[:12], lam=1e-2)
    assert small.fd.basis.n_basis == 4


def test_density_with_fourier_basis() -> None:
    rng = np.random.default_rng(9)
    angles = np.mod(rng.vonmises(1.0, 2.0, size=300), 2 * np.pi)
    basis = Fourier(domain=(0.0, 2 * np.pi), n_basis=5)
    res = fit_density(angles, basis=basis, lam=1e-3)
    assert res.converged
    assert _integral(res, (0.0, 2 * np.pi)) == pytest.approx(1.0, abs=1e-12)
    assert float(res(np.array([1.0]))[0]) > float(res(np.array([1.0 + np.pi]))[0])


def test_density_start_values(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=7)
    base = fit_density(sample, basis=basis, lam=1e-2)
    from_fd = fit_density(sample, basis=basis, lam=1e-2, start=base.fd)
    from_array = fit_density(sample, basis=basis, lam=1e-2, start=np.full(7, 0.3))
    t = np.linspace(0, 1, 9)
    np.testing.assert_allclose(from_fd(t), base(t), rtol=1e-10)
    np.testing.assert_allclose(from_array(t), base(t), rtol=1e-10)
    assert from_fd.n_iter <= 2


def test_nonconvergence_warns(sample: np.ndarray) -> None:
    basis = BSpline(domain=DOMAIN, n_basis=9)
    with pytest.warns(RuntimeWarning, match="did not converge"):
        res = fit_density(sample, basis=basis, lam=1e-4, max_iter=1)
    assert not res.converged
    assert res.n_iter == 1


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"lam": -1.0}, "lam"),
        ({"lam": float("nan")}, "lam"),
        ({"tol": 0.0}, "tol"),
        ({"max_iter": 0}, "max_iter"),
        ({"start": np.zeros(3)}, "start"),
        ({"start": np.full(5, 800.0)}, "overflow"),
        ({"domain": (1.0, 1.0)}, "domain"),
    ],
)
def test_invalid_settings(kwargs: dict[str, Any], message: str, sample: np.ndarray) -> None:
    fit = fit_intensity if message == "overflow" else fit_density
    if "domain" not in kwargs:
        kwargs = {"basis": BSpline(domain=DOMAIN, n_basis=5), **kwargs}
    with pytest.raises(ValueError, match=message):
        fit(sample, **kwargs)


@pytest.mark.parametrize(
    ("x", "message"),
    [
        (np.array([]), "at least one"),
        (np.ones((3, 2)), "one-dimensional"),
        (np.array([0.1, np.nan]), "finite"),
        (np.array([0.1, 1.5]), "basis domain"),
        (np.array([-0.1, 0.5]), "basis domain"),
    ],
)
@pytest.mark.parametrize("fit", [fit_density, fit_intensity])
def test_invalid_samples(x: np.ndarray, message: str, fit: object) -> None:
    with pytest.raises(ValueError, match=message):
        fit(x, basis=BSpline(domain=DOMAIN, n_basis=5))  # type: ignore[operator]


# --------------------------------------------------------------------------- #
# intensity
# --------------------------------------------------------------------------- #


def test_homogeneous_limit_is_the_event_rate() -> None:
    rng = np.random.default_rng(1)
    events = np.cumsum(rng.exponential(scale=0.25, size=120))
    window = (0.0, float(events[-1]))
    basis = BSpline(domain=window, n_basis=8)
    res = fit_intensity(events, basis=basis, lam=1e6, penalty=1)
    assert isinstance(res, IntensityResult)
    rate = events.size / window[1]
    np.testing.assert_allclose(res(np.linspace(*window, 7)), rate, rtol=1e-4)
    assert res.expected_count == pytest.approx(events.size, rel=1e-10)
    assert res.log_likelihood == pytest.approx(
        float(np.sum(res.fd(events)[:, 0])) - res.expected_count
    )


def test_intensity_tracks_a_varying_rate() -> None:
    rng = np.random.default_rng(4)
    candidates = np.cumsum(rng.exponential(scale=1 / 30.0, size=900))
    candidates = candidates[candidates < 20.0]
    keep = rng.uniform(size=candidates.size) < (15.0 + 13.0 * np.sin(candidates)) / 30.0
    events = candidates[keep]
    basis = BSpline(domain=(0.0, 20.0), n_basis=23)
    res = fit_intensity(events, basis=basis, lam=1.0)
    assert res.converged
    peak, trough = res(np.array([np.pi / 2 + 2 * np.pi]))[0], res(np.array([3 * np.pi / 2]))[0]
    assert peak > 3 * trough
    assert res.gradient_norm < 1e-9


def test_intensity_default_window_starts_at_zero() -> None:
    events = np.linspace(1.0, 9.0, 50)
    res = fit_intensity(events, lam=1.0)
    assert res.fd.basis.domain == (0.0, 9.0)
    assert res.fd.basis.n_basis == 5
    shifted = fit_intensity(events - 5.0, lam=1.0)
    assert shifted.fd.basis.domain == (-4.0, 4.0)
    custom = fit_intensity(events, domain=(0.0, 10.0), lam=1.0)
    assert custom.fd.basis.domain == (0.0, 10.0)


def test_intensity_order_of_events_is_irrelevant() -> None:
    rng = np.random.default_rng(2)
    events = np.sort(rng.uniform(0.0, 5.0, size=60))
    basis = BSpline(domain=(0.0, 5.0), n_basis=6)
    a = fit_intensity(events, basis=basis, lam=0.1)
    b = fit_intensity(events[::-1], basis=basis, lam=0.1)
    np.testing.assert_allclose(a.fd.coefs, b.fd.coefs, rtol=1e-12)


def test_intensity_start_as_fdata() -> None:
    events = np.linspace(0.5, 4.5, 30)
    basis = BSpline(domain=(0.0, 5.0), n_basis=5)
    start = FData(np.full(5, np.log(6.0)), basis)
    res = fit_intensity(events, basis=basis, lam=1.0, start=start)
    assert res.converged


# --------------------------------------------------------------------------- #
# properties and backends
# --------------------------------------------------------------------------- #


@settings(max_examples=25, deadline=None)
@given(
    seed=st.integers(min_value=0, max_value=2**31 - 1),
    size=st.integers(min_value=20, max_value=200),
    log_lam=st.floats(min_value=-4.0, max_value=4.0),
)
def test_density_properties(seed: int, size: int, log_lam: float) -> None:
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.0, 1.0, size=size) ** 2
    basis = BSpline(domain=DOMAIN, n_basis=7)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = fit_density(x, basis=basis, lam=10.0**log_lam)
    assert res.converged
    assert _integral(res) == pytest.approx(1.0, abs=1e-10)
    inten = fit_intensity(x, basis=basis, lam=10.0**log_lam)
    t = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose(res(t), inten(t) / size, rtol=1e-7)


def test_torch_input_and_evaluation() -> None:
    torch = pytest.importorskip("torch", reason="torch extra not installed")
    x = torch.tensor([0.1, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.8], dtype=torch.float64)
    basis = BSpline(domain=DOMAIN, n_basis=6)
    dens = fit_density(x, basis=basis, lam=1e-2)
    ref = fit_density(x.numpy(), basis=basis, lam=1e-2)
    t = torch.linspace(0.0, 1.0, 5, dtype=torch.float64)
    values = dens(t)
    assert isinstance(values, torch.Tensor)
    np.testing.assert_allclose(values.numpy(), ref(t.numpy()), rtol=1e-12)
    inten = fit_intensity(x, basis=basis, lam=1e-2)
    assert isinstance(inten(t), torch.Tensor)
