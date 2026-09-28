"""Smoke and vertex-data tests for the plotting helpers (Agg backend, no display needed)."""

from __future__ import annotations

from collections.abc import Iterator

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from fabel import BSpline, FData, Fourier
from fabel._plot import zerofind
from fabel.decomposition import FPCA
from fabel.regression import FRegressResult, fregress
from fabel.smoothing import smooth
from fabel.stats import ConfidenceBand, confidence_band, cycleplot, plot_beta, plot_scores

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


# --------------------------------------------------------------------------- #
# stats plot helpers: bands, plot_beta, cycleplot, plot_scores, zerofind
# --------------------------------------------------------------------------- #


def _model() -> FRegressResult:
    rng = np.random.default_rng(2)
    z = np.linspace(-1.0, 1.0, 30)
    return fregress(2.0 * z + 0.2 * rng.normal(size=30), {"const": 1.0, "z": z})


def test_band_plot_draws_estimate_limits_and_zero_line() -> None:
    t = np.linspace(0.0, 1.0, 6)
    band = ConfidenceBand(t, t - 0.5, 0.1 + 0.0 * t, t - 0.7, t - 0.3, 0.95, name="b")
    ax = band.plot(color="black")
    assert len(ax.lines) == 4
    np.testing.assert_allclose(
        np.asarray(ax.lines[0].get_xydata(), dtype=float), np.column_stack([t, t - 0.5])
    )
    np.testing.assert_allclose(np.asarray(ax.lines[1].get_ydata()), t - 0.7)
    np.testing.assert_allclose(np.asarray(ax.lines[2].get_ydata()), t - 0.3)
    assert ax.lines[1].get_linestyle() == "--"
    assert ax.lines[2].get_color() == ax.lines[0].get_color()
    np.testing.assert_allclose(np.asarray(ax.lines[3].get_ydata()), [0.0, 0.0])
    assert ax.get_title() == "b"


def test_band_plot_omits_zero_line_when_band_is_positive() -> None:
    t = np.linspace(0.0, 1.0, 4)
    band = ConfidenceBand(t, t + 2.0, 0.1 + 0.0 * t, t + 1.0, t + 3.0, 0.9)
    _, ax = plt.subplots()
    assert band.plot(ax=ax) is ax
    assert len(ax.lines) == 3
    assert ax.get_title() == ""


def test_band_plot_of_a_smooth_draws_every_curve(fd: FData) -> None:
    t = np.linspace(0.0, 1.0, 20)
    fit = smooth(fd(t), t, basis=fd.basis, lam=1e-6)
    band = confidence_band(fit, sigma_e=0.01)
    ax = band.plot()
    crosses = zerofind([band.lower.min(), band.upper.max()])
    assert len(ax.lines) == 3 * fd.n_curves + int(crosses)
    for index in range(fd.n_curves):
        np.testing.assert_allclose(np.asarray(ax.lines[index].get_ydata()), band.estimate[:, index])
        np.testing.assert_allclose(
            np.asarray(ax.lines[fd.n_curves + index].get_ydata()), band.lower[:, index]
        )


def test_plot_beta_one_panel_per_term() -> None:
    model = _model()
    panels = plot_beta(model, np.linspace(0.0, 1.0, 5), level=0.9)
    assert [ax.get_title() for ax in panels] == ["const", "z"]
    bands = confidence_band(model, np.linspace(0.0, 1.0, 5), level=0.9)
    np.testing.assert_allclose(np.asarray(panels[1].lines[1].get_ydata()), bands[1].lower)


def test_plot_beta_accepts_bands_and_given_axes() -> None:
    bands = confidence_band(_model())
    _, (first, second) = plt.subplots(2, 1)
    assert plot_beta(bands, axes=[first, second]) == [first, second]
    assert len(plot_beta(bands[1])) == 1
    with pytest.raises(ValueError, match="one axes per band"):
        plot_beta(bands, axes=[first])


def _loops() -> FData:
    basis = Fourier(domain=(0.0, 1.0), n_basis=3)
    coefs = np.zeros((3, 2, 2))
    coefs[1, :, 0] = [1.0, 2.0]
    coefs[2, :, 1] = [1.0, 0.5]
    return FData(coefs, basis)


def test_cycleplot_traces_each_bivariate_curve() -> None:
    loops = _loops()
    ax = cycleplot(loops, n_points=50)
    assert len(ax.lines) == 2
    grid = np.linspace(0.0, 1.0, 50)
    values = loops(grid)
    for index, line in enumerate(ax.lines):
        np.testing.assert_allclose(
            np.asarray(line.get_xydata(), dtype=float),
            np.column_stack([values[:, index, 0], values[:, index, 1]]),
        )
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("x(t)", "y(t)")
    xy = np.asarray(ax.lines[0].get_xydata())
    np.testing.assert_allclose(xy[0], xy[-1], atol=1e-12)


def test_cycleplot_accepts_two_univariate_coordinates(fd: FData) -> None:
    _, ax = plt.subplots()
    out = cycleplot(fd, fd * 2.0, ax=ax, n_points=11, linestyle=":")
    assert out is ax
    assert len(ax.lines) == fd.n_curves
    line = np.asarray(ax.lines[0].get_xydata())
    np.testing.assert_allclose(line[:, 1], 2.0 * line[:, 0], rtol=1e-12, atol=1e-12)
    assert ax.lines[0].get_linestyle() == ":"


def test_cycleplot_rejects_bad_input(fd: FData) -> None:
    with pytest.raises(ValueError, match="two variables"):
        cycleplot(fd)
    with pytest.raises(ValueError, match="one variable"):
        cycleplot(_loops(), fd)
    with pytest.raises(ValueError, match="same number of curves"):
        cycleplot(fd, fd[0])
    with pytest.raises(ValueError, match="n_points"):
        cycleplot(_loops(), n_points=1)


def test_plot_scores_scatters_the_chosen_components() -> None:
    scores = RNG.normal(size=(6, 3))
    ax = plot_scores(scores, (2, 0), labels=list("abcdef"), marker="x")
    np.testing.assert_allclose(
        np.asarray(ax.collections[0].get_offsets(), dtype=float), scores[:, [2, 0]]
    )
    assert [text.get_text() for text in ax.texts] == list("abcdef")
    assert ax.get_xlabel() == "Harmonic 3 score"


def test_plot_scores_reads_a_fitted_fpca() -> None:
    basis = BSpline(domain=(0.0, 1.0), n_basis=7)
    curves = FData(RNG.normal(size=(7, 12)), basis)
    fpca = FPCA(n=2).fit(curves)
    _, ax = plt.subplots()
    assert plot_scores(fpca, ax=ax) is ax
    np.testing.assert_allclose(
        np.asarray(ax.collections[0].get_offsets(), dtype=float), fpca.scores
    )
    share = 100.0 * float(fpca.varprop[0])
    assert ax.get_xlabel() == f"Harmonic 1 score ({share:.1f}%)"


def test_plot_scores_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="matrix"):
        plot_scores(np.ones(4))
    with pytest.raises(ValueError, match="out of range"):
        plot_scores(np.ones((4, 2)), (0, 2))
    with pytest.raises(ValueError, match="labels"):
        plot_scores(np.ones((4, 2)), labels=["a"])


def test_zerofind() -> None:
    assert zerofind([1.0, 5.0]) is False
    assert zerofind([0.0, 3.0]) is True
    assert zerofind(np.array([[-2.0], [-1.0]])) is False
    assert zerofind([-1.0, 1.0]) is True
    assert zerofind([]) is False
