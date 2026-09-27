"""Smoke tests for the plotting mixin (Agg backend, no display needed)."""

from __future__ import annotations

from collections.abc import Iterator

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from fabel import BSpline, FData

RNG = np.random.default_rng(7)


@pytest.fixture
def fd() -> FData:
    """Return three random cubic spline curves."""
    basis = BSpline(domain=(0.0, 1.0), n_basis=7)
    return FData(RNG.normal(size=(7, 3)), basis)


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    """Close every figure the test opened."""
    yield
    plt.close("all")


def test_plot_returns_axes_with_one_line_per_curve(fd: FData) -> None:
    ax = fd.plot()
    assert isinstance(ax, Axes)
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


def test_plot_line_data_are_the_curve_values(fd: FData) -> None:
    """The drawn vertices reproduce the curve exactly.

    A stable substitute for an image hash: it pins the data the renderer is
    given without depending on the matplotlib version, fonts or the platform.
    """
    ax = fd.plot()
    for index, line in enumerate(ax.lines):
        t, values = np.asarray(line.get_xydata()).T
        np.testing.assert_allclose(values, fd(t)[:, index], rtol=1e-12, atol=1e-12)


def test_plot_deriv_line_data_are_the_derivative_values(fd: FData) -> None:
    ax = fd.plot(deriv=2)
    t = np.asarray(ax.lines[0].get_xydata())[:, 0]
    expected = fd(t, 2)
    for index, line in enumerate(ax.lines):
        np.testing.assert_allclose(
            np.asarray(line.get_ydata()), expected[:, index], rtol=1e-12, atol=1e-12
        )
