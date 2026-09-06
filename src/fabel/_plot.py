"""Matplotlib helpers mixed into :class:`fabel.core.FData`.

Matplotlib is imported lazily inside the methods so that importing ``fabel``
never pulls in a plotting stack.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from fabel._backend import default_namespace, to_numpy

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

__all__ = ["PlotMixin"]

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
