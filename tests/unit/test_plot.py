"""Smoke tests for the plotting mixin (Agg backend, no display needed)."""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from fabel import BSpline, FData

RNG = np.random.default_rng(7)


@pytest.fixture
def fd() -> FData:
    """Return three random cubic spline curves."""
    basis = BSpline(domain=(0.0, 1.0), n_basis=7)
    return FData(RNG.normal(size=(7, 3)), basis)


@pytest.fixture(autouse=True)
def _close_figures() -> None:
    """Close every figure the test opened."""
    yield
    plt.close("all")


def test_plot_returns_axes_with_one_line_per_curve(fd: FData) -> None:
    ax = fd.plot()
    assert isinstance(ax, plt.Axes)
    assert len(ax.lines) == fd.n_curves
    assert ax.get_xlabel() == "t"


def test_plot_draws_on_the_given_axes(fd: FData) -> None:
    _, ax = plt.subplots()
    assert fd.plot(ax=ax) is ax


def test_plot_accepts_a_derivative_and_style_kwargs(fd: FData) -> None:
    ax = fd.plot(deriv=1, linestyle="--")
    assert ax.lines[0].get_linestyle() == "--"


def test_plot_handles_multivariate_coefficients() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    ax = FData(RNG.normal(size=(5, 2, 3)), basis).plot()
    assert len(ax.lines) == 6


def test_plot_fit_overlays_the_observations(fd: FData) -> None:
    t = np.linspace(0.0, 1.0, 15)
    ax = fd.plot_fit(fd(t), t)
    assert len(ax.lines) == 2 * fd.n_curves
    assert ax.lines[-1].get_linestyle() == "None"


def test_plot_fit_accepts_one_dimensional_data() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=5)
    curve = FData(RNG.normal(size=(5, 1)), basis)
    t = np.linspace(0.0, 1.0, 9)
    ax = curve.plot_fit(curve(t)[:, 0], t)
    assert len(ax.lines) == 2
