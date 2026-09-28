"""Matplotlib helpers mixed into :class:`fabel.core.FData`.

Matplotlib is imported lazily inside the methods so that importing ``fabel``
never pulls in a plotting stack.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from fabel._backend import default_namespace, to_numpy

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

__all__ = ["PlotMixin", "zerofind"]

_PLOT_POINTS = 201


class PlotMixin:
    """Provide ``.plot()`` and ``.plot_fit()`` to a functional object."""

    def _plot_grid(self, n: int = _PLOT_POINTS) -> Any:
        """Return an evaluation grid spanning the domain."""
        xp = default_namespace()
        lower, upper = self.domain  # type: ignore[attr-defined]
        return xp.linspace(lower, upper, n, dtype=xp.float64)

    def plot(self, ax: Axes | None = None, *, deriv: int = 0, **kwargs: Any) -> Axes:
        """Plot every curve against its argument.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when omitted.
        deriv : int, optional
            Derivative to plot.  Defaults to ``0``.
        **kwargs
            Passed through to :meth:`matplotlib.axes.Axes.plot`.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((5, 2)), fb.BSpline(domain=(0.0, 1.0), n_basis=5))
        >>> type(fd.plot()).__name__
        'Axes'
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        grid = self._plot_grid()
        values = to_numpy(self(grid, deriv=deriv))  # type: ignore[operator]
        ax.plot(to_numpy(grid), values.reshape(values.shape[0], -1), **kwargs)
        ax.set_xlabel("t")
        ax.set_ylabel("value")
        return ax

    def plot_fit(self, y: Any, t: Any, ax: Axes | None = None, **kwargs: Any) -> Axes:
        """Plot the curves together with the data they were fitted to.

        Parameters
        ----------
        y : array_like
            Observed values, shape ``(len(t),)`` or ``(len(t), n_curves)``.
        t : array_like
            Observation points, shape ``(n_points,)``.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when omitted.
        **kwargs
            Passed through to :meth:`matplotlib.axes.Axes.plot` for the curves.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((5, 1)), fb.BSpline(domain=(0.0, 1.0), n_basis=5))
        >>> t = np.linspace(0.0, 1.0, 7)
        >>> len(fd.plot_fit(np.ones(7), t).lines) > 1
        True
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        points = to_numpy(t)
        observed = to_numpy(y).reshape(points.shape[0], -1)
        self.plot(ax=ax, **kwargs)
        ax.plot(points, observed, "o", markersize=3.0, fillstyle="none")
        return ax


# --------------------------------------------------------------------------- #
# standalone plot helpers behind fabel.stats (cycleplot, plot_scores, bands)
# --------------------------------------------------------------------------- #


def _axes(ax: Axes | None) -> Axes:
    """Return ``ax``, or the axes of a new figure when it is ``None``."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    return ax


def zerofind(values: Any) -> bool:
    """Return whether the range of ``values`` contains zero (R's ``zerofind``).

    Parameters
    ----------
    values : array_like
        Any collection of numbers.

    Returns
    -------
    bool
        ``True`` when ``min(values) <= 0 <= max(values)``.

    Examples
    --------
    >>> zerofind([1.0, 2.0, 5.0])
    False
    >>> zerofind([-1.0, 3.0])
    True
    """
    xp = default_namespace()
    flat = xp.reshape(xp.asarray(to_numpy(values), dtype=xp.float64), (-1,))
    if flat.shape[0] == 0:
        return False
    return bool(xp.min(flat) <= 0.0 <= xp.max(flat))


def cycle_lines(x: Any, y: Any, ax: Axes | None, **kwargs: Any) -> Axes:
    """Draw ``y`` against ``x`` column by column: one closed cycle per curve.

    ``x`` and ``y`` are ``(n_points, n_curves)`` arrays of values sampled on a
    common grid.
    """
    ax = _axes(ax)
    ax.plot(to_numpy(x), to_numpy(y), **kwargs)
    return ax


def score_scatter(
    scores: Any,
    components: tuple[int, int],
    ax: Axes | None,
    labels: Any,
    xlabel: str,
    ylabel: str,
    **kwargs: Any,
) -> Axes:
    """Scatter two columns of a score matrix, optionally labelling every point."""
    ax = _axes(ax)
    values = to_numpy(scores)
    first, second = components
    x, y = values[:, first], values[:, second]
    ax.scatter(x, y, **kwargs)
    if labels is not None:
        for label, px, py in zip(labels, x.tolist(), y.tolist(), strict=True):
            ax.annotate(str(label), (px, py), fontsize="small")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax


def band_lines(
    t: Any,
    estimate: Any,
    lower: Any,
    upper: Any,
    ax: Axes | None,
    title: str | None,
    **kwargs: Any,
) -> Axes:
    """Draw an estimate with its pointwise lower and upper limits.

    The estimate lines come first, then the lower limits, then the upper
    limits (dashed, in the estimate's colour); a dotted zero line is added last
    when the band crosses zero, as R's ``plotbeta`` does.
    """
    ax = _axes(ax)
    points = to_numpy(t)
    width = points.shape[0]
    centre = to_numpy(estimate).reshape(width, -1)
    low = to_numpy(lower).reshape(width, -1)
    high = to_numpy(upper).reshape(width, -1)
    drawn = ax.plot(points, centre, **kwargs)
    colours = [line.get_color() for line in drawn]
    for limits in (low, high):
        for column, colour in enumerate(colours):
            ax.plot(points, limits[:, column], linestyle="--", color=colour)
    if zerofind([low.min(), high.max()]):
        ax.axhline(0.0, color="grey", linestyle=":", linewidth=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel("value")
    if title is not None:
        ax.set_title(title)
    return ax


def panel_axes(count: int) -> list[Axes]:
    """Return ``count`` axes stacked in one new figure."""
    import matplotlib.pyplot as plt

    _, grid = plt.subplots(count, 1, squeeze=False, figsize=(6.4, 3.0 * count))
    return [grid[row, 0] for row in range(count)]
