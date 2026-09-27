r"""Functional statistics: covariance, correlation, depth, boxplots and tests.

Replaces R's ``var.fd``, ``cor.fd``, ``fbplot``, ``fdepth``, ``tperm.fd`` and
``Fperm.fd``.

* :func:`cov` and :func:`cor` give the (cross-)covariance surface and the
  correlation matrix of one or two sets of curves.
* :func:`depth` orders curves from the centre outwards (modified band depth,
  band depth or Fraiman--Muniz depth) and returns the deepest curve and a
  trimmed mean; :func:`boxplot` builds the functional boxplot of Sun and
  Genton (2011) on top of it.
* :func:`t_test` and :func:`f_test` are permutation tests whose statistic is
  the largest pointwise two-sample *t* or regression *F* statistic.

Depths and tests work on curve values sampled on a grid, so they return NumPy
arrays whatever the input namespace; :func:`cov` and :func:`cor` stay in the
namespace of their inputs.

Examples
--------
>>> import numpy as np
>>> from fabel import BSpline, FData
>>> from fabel.stats import depth
>>> fd = FData(np.random.default_rng(0).normal(size=(6, 9)), BSpline(n_basis=6))
>>> depth(fd).depth.shape
(9,)
"""

from __future__ import annotations

import math
import operator
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from fabel import _linalg
from fabel._backend import (
    array_namespace,
    asarray,
    default_namespace,
    result_namespace,
    to_numpy,
)
from fabel._operator import LDO
from fabel.basis import Basis, _same_domain
from fabel.core import BiFData, FData, _quadrature, inprod

if TYPE_CHECKING:  # pragma: no cover - typing only
    import numpy as np
    from matplotlib.axes import Axes

    Array = Any
    NDArray = np.ndarray[Any, np.dtype[Any]]
    RandomState = int | np.random.Generator | None
else:
    Array = Any
    NDArray = Any
    RandomState = Any

__all__ = [
    "BoxplotResult",
    "DepthResult",
    "PermutationTestResult",
    "boxplot",
    "cor",
    "cov",
    "depth",
    "f_test",
    "t_test",
]

#: Evaluation points used when curves are given as :class:`FData` and no grid
#: is passed -- the default of R's ``tperm.fd`` and ``Fperm.fd``.
_N_POINTS = 101

#: Depth notions understood by :func:`depth` and :func:`boxplot`.
_METHODS = ("MBD", "BD2", "BOTH", "FM")

#: Weight of band depth against modified band depth in the combined "both"
#: ordering: band depth rounded to four decimals, ties broken by MBD.
_BOTH_SCALE = 1e4


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, eq=False)
class DepthResult:
    """Depth of every curve, with the deepest curve and a trimmed mean.

    Attributes
    ----------
    depth : numpy.ndarray
        Depth of each curve, shape ``(n_curves,)``; larger is more central.
    median_index : int
        Position of the deepest curve (the first one if several tie).
    median : numpy.ndarray
        Values of the deepest curve at the evaluation points.
    trimmed_index : numpy.ndarray
        Positions, in increasing order, of the curves whose depth is at least
        the ``trim`` quantile of all depths.
    trimmed_mean : numpy.ndarray
        Pointwise mean of the curves in ``trimmed_index``.
    t : numpy.ndarray or None
        The evaluation points, or ``None`` when values were given without them.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.stats import depth
    >>> values = np.array([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])
    >>> depth(values).median_index
    1
    """

    depth: NDArray
    median_index: int
    median: NDArray
    trimmed_index: NDArray
    trimmed_mean: NDArray
    t: NDArray | None


@dataclass(frozen=True, eq=False)
class BoxplotResult:
    """A functional boxplot: depth ordering, central region, fences and outliers.

    Attributes
    ----------
    values : numpy.ndarray
        The curves as sampled, shape ``(n_points, n_curves)``.
    t : numpy.ndarray or None
        The evaluation points, or ``None`` when values were given without them.
    depth : numpy.ndarray
        Depth of each curve, shape ``(n_curves,)``.
    median_index : int
        Position of the deepest curve.
    median : numpy.ndarray
        Values of the deepest curve.
    outliers : numpy.ndarray
        Positions of the curves leaving the fences at some point, increasing.
    central_lower, central_upper : numpy.ndarray
        Pointwise envelope of the central region, the ``ceil(prob * n_curves)``
        deepest curves -- the box.
    fence_lower, fence_upper : numpy.ndarray
        The box inflated by ``factor`` times its height on either side.
    whisker_lower, whisker_upper : numpy.ndarray
        Pointwise envelope of every curve that is not an outlier.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.stats import boxplot
    >>> values = np.vstack([np.arange(7.0), np.arange(7.0)])
    >>> values[:, 6] = 40.0
    >>> boxplot(values).outliers.tolist()
    [6]
    """

    values: NDArray
    t: NDArray | None
    depth: NDArray
    median_index: int
    median: NDArray
    outliers: NDArray
    central_lower: NDArray
    central_upper: NDArray
    fence_lower: NDArray
    fence_upper: NDArray
    whisker_lower: NDArray
    whisker_upper: NDArray

    def plot(self, ax: Axes | None = None, **kwargs: Any) -> Axes:
        """Draw the boxplot: central region, median, whiskers and outliers.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when omitted.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.plot` for the median curve.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on.  The median is the first line, then the upper
            and lower whiskers, then one line per outlier; the central region
            is the single filled collection.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> from fabel.stats import boxplot
        >>> values = np.vstack([np.arange(7.0), np.arange(7.0) + 1.0])
        >>> len(boxplot(values).plot().lines)
        3
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        nxp = default_namespace()
        axis = self.t if self.t is not None else nxp.arange(1, self.values.shape[0] + 1)
        ax.fill_between(axis, self.central_lower, self.central_upper, color="violet", alpha=0.6)
        kwargs.setdefault("color", "black")
        ax.plot(axis, self.median, **kwargs)
        ax.plot(axis, self.whisker_upper, color="tab:blue")
        ax.plot(axis, self.whisker_lower, color="tab:blue")
        for index in self.outliers.tolist():
            ax.plot(axis, self.values[:, index], color="tab:red", linestyle="--")
        ax.set_xlabel("t")
        ax.set_ylabel("value")
        return ax


@dataclass(frozen=True, eq=False)
class PermutationTestResult:
    """Outcome of a permutation test built on a pointwise statistic.

    Attributes
    ----------
    statistic : float
        The observed statistic: the maximum of ``pointwise``.
    null : numpy.ndarray
        The statistic recomputed under each permutation, shape ``(n_perm,)``.
    pvalue : float
        Share of permutations whose statistic reaches the observed one.
    critical_value : float
        The ``1 - q`` quantile of ``null`` (R's ``qval``).
    pointwise : numpy.ndarray
        Observed statistic at each evaluation point, shape ``(n_points,)``; a
        single entry for a scalar response.
    pointwise_null : numpy.ndarray
        Pointwise statistic under each permutation, ``(n_perm, n_points)``.
    pointwise_pvalue : numpy.ndarray
        Pointwise permutation p-values.
    pointwise_critical_value : numpy.ndarray
        Pointwise ``1 - q`` quantiles of the null.
    t : numpy.ndarray or None
        The evaluation points; ``None`` for a scalar response.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.stats import t_test
    >>> rng = np.random.default_rng(0)
    >>> a = FData(rng.normal(size=(5, 6)), BSpline(n_basis=5))
    >>> b = FData(rng.normal(size=(5, 6)) + 3.0, BSpline(n_basis=5))
    >>> t_test(a, b, n_perm=50, random_state=1).pvalue
    0.0
    """

    statistic: float
    null: NDArray
    pvalue: float
    critical_value: float
    pointwise: NDArray
    pointwise_null: NDArray
    pointwise_pvalue: NDArray
    pointwise_critical_value: NDArray
    t: NDArray | None


# --------------------------------------------------------------------------- #
# covariance and correlation
# --------------------------------------------------------------------------- #


def _check_pair(x: FData, y: FData) -> None:
    """Validate two sets of curves for a cross-moment."""
    if x.n_vars != 1 or y.n_vars != 1:
        raise ValueError("cross moments need curves with a single variable each")
    if x.n_curves != y.n_curves:
        raise ValueError(
            f"cross moments need the same number of curves, got {x.n_curves} and {y.n_curves}"
        )
    if x.n_curves < 2:
        raise ValueError("a covariance needs at least two curves")
    if not _same_domain(x.domain, y.domain):
        raise ValueError(f"the curves must share a domain, got {x.domain} and {y.domain}")


def _default_grid(domain: tuple[float, float]) -> Array:
    """Return ``_N_POINTS`` equally spaced points spanning ``domain``."""
    nxp = default_namespace()
    return nxp.linspace(domain[0], domain[1], _N_POINTS, dtype=nxp.float64)


def cov(x: FData, y: FData | None = None) -> BiFData:
    """Return the sample (cross-)covariance surface.

    Replaces R's ``var.fd``.  With one argument this is :meth:`FData.cov`.  With
    two, it is the bivariate function

    ``c(s, t) = 1 / (n - 1) sum_i (x_i(s) - xbar(s)) (y_i(t) - ybar(t))``

    expanded in ``x.basis`` along ``s`` and ``y.basis`` along ``t``.  Both are
    exact: centring and the outer product act on the coefficients.

    Parameters
    ----------
    x : FData
        First set of curves.
    y : FData, optional
        Second set of curves, paired with ``x`` curve by curve.

    Returns
    -------
    BiFData
        The covariance surface.

    Raises
    ------
    ValueError
        If fewer than two curves are given, or ``x`` and ``y`` differ in their
        number of curves or domain, or carry several variables.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.stats import cov
    >>> rng = np.random.default_rng(0)
    >>> x = FData(rng.normal(size=(5, 8)), BSpline(n_basis=5))
    >>> y = FData(rng.normal(size=(4, 8)), BSpline(n_basis=4))
    >>> cov(x, y).coefs.shape
    (5, 4)
    """
    if y is None:
        return x.cov()
    _check_pair(x, y)
    xp = result_namespace(x.coefs, y.coefs)
    left = asarray(x.center().coefs, xp=xp)
    right = asarray(y.center().coefs, xp=xp)
    coefs = xp.matmul(left, xp.matrix_transpose(right)) / (x.n_curves - 1)
    return BiFData(coefs, x.basis, y.basis)


def _centred(values: Array) -> tuple[Array, Array]:
    """Return values centred over the curve axis, with variables moved first."""
    xp = array_namespace(values)
    stacked = values if len(values.shape) == 3 else values[:, :, None]
    stacked = xp.permute_dims(stacked, (2, 0, 1))  # (n_vars, n_points, n_curves)
    centred = stacked - xp.mean(stacked, axis=2, keepdims=True)
    n_curves = stacked.shape[2]
    variance = xp.sum(centred * centred, axis=2) / (n_curves - 1)
    return centred, variance


def cor(
    x: FData,
    y: FData | None = None,
    *,
    s: Any = None,
    t: Any = None,
) -> Array:
    """Return the sample correlation of ``x(s)`` with ``y(t)`` on a grid.

    Replaces R's ``cor.fd``.  Entry ``[a, b]`` is the correlation, across
    curves, of ``x(s[a])`` and ``y(t[b])``; with ``y`` omitted it is the
    autocorrelation of ``x``.  A correlation is a ratio of covariances, not a
    basis expansion, so it is returned as values rather than as a
    :class:`~fabel.core.BiFData`.

    Parameters
    ----------
    x : FData
        First set of curves.
    y : FData, optional
        Second set of curves, paired with ``x``.  Defaults to ``x``.
    s, t : array_like, optional
        Evaluation points for ``x`` and ``y``.  Each defaults to 101 equally
        spaced points over the domain.

    Returns
    -------
    array
        Shape ``(len(s), len(t))``.  When ``y`` is omitted and ``x`` holds
        several variables, each variable gets its own matrix, stacked on a
        trailing axis.

    Raises
    ------
    ValueError
        If fewer than two curves are given, or ``x`` and ``y`` cannot be paired.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.stats import cor
    >>> fd = FData(np.random.default_rng(0).normal(size=(5, 8)), BSpline(n_basis=5))
    >>> grid = np.linspace(0.0, 1.0, 3)
    >>> np.round(np.diag(cor(fd, s=grid, t=grid)), 12).tolist()
    [1.0, 1.0, 1.0]
    """
    if y is None:
        if x.n_curves < 2:
            raise ValueError("a correlation needs at least two curves")
        y_fd = x
    else:
        _check_pair(x, y)
        y_fd = y
    s_points = _default_grid(x.domain) if s is None else s
    t_points = _default_grid(y_fd.domain) if t is None else t
    left, left_var = _centred(x(s_points))
    right, right_var = _centred(y_fd(t_points))
    xp = result_namespace(left, right)
    left, left_var = asarray(left, xp=xp), asarray(left_var, xp=xp)
    right, right_var = asarray(right, xp=xp), asarray(right_var, xp=xp)
    covariance = xp.matmul(left, xp.matrix_transpose(right)) / (x.n_curves - 1)
    scale = xp.sqrt(left_var[:, :, None] * right_var[:, None, :])
    matrix = xp.permute_dims(covariance / scale, (1, 2, 0))
    return matrix[:, :, 0] if x.n_vars == 1 else matrix


# --------------------------------------------------------------------------- #
# depth
# --------------------------------------------------------------------------- #


def _sampled(data: FData | Any, t: Any) -> tuple[NDArray, NDArray | None]:
    """Return curve values ``(n_points, n_curves)`` in NumPy, and their points."""
    nxp = default_namespace()
    points: NDArray | None
    if isinstance(data, FData):
        if data.n_vars != 1:
            raise ValueError("depth is defined for curves with a single variable")
        points = _default_grid(data.domain) if t is None else to_numpy(t)
        values = to_numpy(data(points))
    else:
        values = nxp.asarray(to_numpy(data), dtype=nxp.float64)
        if len(values.shape) != 2:
            raise ValueError(
                f"values must be two-dimensional (n_points, n_curves), got shape {values.shape}"
            )
        points = None if t is None else to_numpy(t)
        if points is not None and points.shape != (values.shape[0],):
            raise ValueError(
                f"{values.shape[0]} rows of values need as many evaluation points, "
                f"got shape {points.shape}"
            )
    if not bool(nxp.all(nxp.isfinite(values))):
        raise ValueError("curve values must be finite")
    return nxp.asarray(values, dtype=nxp.float64), points


def _rank_counts(values: NDArray) -> tuple[NDArray, NDArray]:
    """Return, per entry, how many curves lie strictly below and at-or-below it.

    Counts are taken point by point across curves, as integers, so that every
    depth built from them is an exact ratio of integers and ties stay ties.
    """
    nxp = default_namespace()
    ordered = nxp.sort(values, axis=1)
    below = nxp.empty(values.shape, dtype=nxp.int64)
    at_or_below = nxp.empty(values.shape, dtype=nxp.int64)
    for row in range(values.shape[0]):
        below[row] = nxp.searchsorted(ordered[row], values[row], side="left")
        at_or_below[row] = nxp.searchsorted(ordered[row], values[row], side="right")
    return below, at_or_below


def _depth_values(values: NDArray, method: str) -> Array:
    """Return the depth of each column of ``values`` under ``method``."""
    nxp = default_namespace()
    n_points, n = values.shape
    below, at_or_below = _rank_counts(values)
    if method == "FM":
        # 1 - |1/2 - F_n(x)| with F_n(x) = at_or_below / n, times 2n: an integer.
        centrality = 2 * n - nxp.abs(n - 2 * at_or_below)
        return nxp.sum(centrality, axis=0) / (2.0 * n * n_points)
    if n < 2:
        raise ValueError("band depths need at least two curves")
    # Twice the average rank of each value among the n values at its point.
    twice_rank = below + at_or_below + 1
    if method == "BD2":
        return _band_depth(twice_rank, n)
    mbd = _modified_band_depth(twice_rank, n, n_points)
    if method == "MBD":
        return mbd
    return nxp.round(_band_depth(twice_rank, n) * _BOTH_SCALE) + mbd


def _modified_band_depth(twice_rank: NDArray, n: int, n_points: int) -> Array:
    """Mean share of two-curve bands containing each value (rank form).

    A value of average rank ``r`` is inside the band of every pair that has one
    member at or below it and one at or above it: ``(r - 1)(n - r)`` pairs of
    other curves plus the ``n - 1`` pairs it belongs to.  Four times that count
    is an integer, so the sum over points is exact.
    """
    nxp = default_namespace()
    count4 = (twice_rank - 2) * (2 * n - twice_rank) + 4 * (n - 1)
    return nxp.sum(count4, axis=0) / (2.0 * n_points * n * (n - 1))


def _band_depth(twice_rank: NDArray, n: int) -> Array:
    """Rank-based two-curve band depth of Sun, Genton and Nychka (2012)."""
    nxp = default_namespace()
    lowest = nxp.min(twice_rank, axis=0)
    highest = nxp.max(twice_rank, axis=0)
    count4 = (2 * n - highest) * (lowest - 2) + 4 * (n - 1)
    return count4 / (2.0 * n * (n - 1))


def _quantile(values: NDArray, prob: float) -> Array:
    """Return the type-7 quantile of ``values`` along axis 0.

    ``Q(p) = (1 - g) x_(j) + g x_(j+1)`` with ``j + g = (m - 1) p``, the sample
    quantile of Hyndman and Fan (1996) that R uses by default.  Where the two
    order statistics are equal the result is that value exactly, so a quantile
    falling on tied values compares equal to them.
    """
    nxp = default_namespace()
    ordered = nxp.sort(values, axis=0)
    size = ordered.shape[0]
    position = (size - 1) * prob
    lower = min(math.floor(position), size - 1)
    upper = min(lower + 1, size - 1)
    fraction = position - lower
    low, high = ordered[lower], ordered[upper]
    blended = (1.0 - fraction) * low + fraction * high
    return nxp.where(high == low, low, blended)


def _normalise_method(method: str) -> str:
    """Return the canonical spelling of a depth method name."""
    key = str(method).upper()
    if key not in _METHODS:
        raise ValueError(f"method must be one of 'MBD', 'BD2', 'both' or 'FM', got {method!r}")
    return key


def depth(
    data: FData | Any,
    *,
    method: str = "MBD",
    t: Any = None,
    trim: float = 0.25,
) -> DepthResult:
    r"""Rank curves from the centre outwards by a functional depth.

    Replaces R's ``fdepth`` (``type = "FM"``) and the depth half of ``fbplot``.
    With ``r_i(t)`` the rank of ``x_i(t)`` among the ``n`` values at ``t``
    (average ranks for ties):

    * ``"MBD"`` -- modified band depth (Lopez-Pintado and Romo, 2009): the
      share of bands spanned by pairs of curves that contain ``x_i(t)``,
      averaged over the points, ``mean_t [(r - 1)(n - r) + n - 1] / C(n, 2)``.
    * ``"BD2"`` -- the rank-based two-curve band depth of Sun, Genton and
      Nychka (2012), ``[(n - max_t r)(min_t r - 1) + n - 1] / C(n, 2)``.  It
      equals the share of bands containing the whole curve when the curves do
      not cross, and is what R's ``fbplot(method = "BD2")`` reports.
    * ``"both"`` -- band depth rounded to four decimals, ties broken by modified
      band depth (``round(1e4 BD2) + MBD``), as ``fbplot(method = "Both")``.
    * ``"FM"`` -- Fraiman and Muniz (2001) depth,
      ``mean_t [1 - |1/2 - F_n(x_i(t))|]`` with ``F_n`` the empirical
      distribution function at ``t``.

    Averages over ``t`` weight every evaluation point equally.

    Parameters
    ----------
    data : FData or array_like
        Curves, or their values as an ``(n_points, n_curves)`` array.
    method : {"MBD", "BD2", "both", "FM"}, optional
        Depth notion, case-insensitive.  Defaults to ``"MBD"``.
    t : array_like, optional
        Evaluation points.  For :class:`FData` they default to 101 equally
        spaced points over the domain; for an array they are only recorded.
    trim : float, optional
        Share of least deep curves left out of the trimmed mean: curves whose
        depth reaches the ``trim`` quantile are kept.  Default ``0.25``.

    Returns
    -------
    DepthResult
        Depths, the deepest curve and the trimmed mean.

    Raises
    ------
    ValueError
        If ``method`` or ``trim`` is invalid, the values are not a finite
        two-dimensional array matching ``t``, or a band depth is asked of fewer
        than two curves.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.stats import depth
    >>> values = np.array([[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]])
    >>> depth(values, method="FM").depth.tolist()
    [0.75, 1.0, 0.75, 0.5]
    """
    key = _normalise_method(method)
    if not 0.0 <= trim < 1.0:
        raise ValueError(f"trim must lie in [0, 1), got {trim}")
    values, points = _sampled(data, t)
    scores = _depth_values(values, key)
    nxp = default_namespace()
    median_index = int(nxp.argmax(scores))
    kept = nxp.nonzero(scores >= _quantile(scores, trim))[0]
    return DepthResult(
        depth=scores,
        median_index=median_index,
        median=values[:, median_index],
        trimmed_index=kept,
        trimmed_mean=nxp.mean(values[:, kept], axis=1),
        t=points,
    )


def boxplot(
    data: FData | Any,
    *,
    method: str = "MBD",
    t: Any = None,
    prob: float = 0.5,
    factor: float = 1.5,
) -> BoxplotResult:
    """Build the functional boxplot of Sun and Genton (2011).

    Replaces R's ``fbplot`` and ``boxplot.fd``.  Curves are ordered by
    :func:`depth`.  The central region is the pointwise envelope of the
    ``ceil(prob * n_curves)`` deepest curves; inflating it by ``factor`` times
    its height on both sides gives the fences, and a curve that leaves the
    fences anywhere is an outlier.  The whiskers are the envelope of the curves
    that are not outliers.

    Parameters
    ----------
    data : FData or array_like
        Curves, or their values as an ``(n_points, n_curves)`` array.
    method : {"MBD", "BD2", "both", "FM"}, optional
        Depth used to order the curves.  Defaults to ``"MBD"``.
    t : array_like, optional
        Evaluation points, as for :func:`depth`.
    prob : float, optional
        Share of curves in the central region.  Default ``0.5``.
    factor : float, optional
        Fence inflation factor.  Default ``1.5``, as in a classical boxplot.

    Returns
    -------
    BoxplotResult
        Depths, median, envelopes and outliers; call ``.plot()`` to draw it.

    Raises
    ------
    ValueError
        If ``prob`` is not in ``(0, 1]``, ``factor`` is negative, or the data
        are invalid (see :func:`depth`).

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.stats import boxplot
    >>> rng = np.random.default_rng(1)
    >>> values = rng.normal(size=(20, 15))
    >>> values[:, 4] += 25.0
    >>> boxplot(values).outliers.tolist()
    [4]
    """
    key = _normalise_method(method)
    if not 0.0 < prob <= 1.0:
        raise ValueError(f"prob must lie in (0, 1], got {prob}")
    if factor < 0.0:
        raise ValueError(f"factor must be non-negative, got {factor}")
    values, points = _sampled(data, t)
    scores = _depth_values(values, key)
    nxp = default_namespace()
    n = values.shape[1]
    order = nxp.argsort(-scores, stable=True)
    central = values[:, order[: max(1, math.ceil(prob * n))]]
    lower = nxp.min(central, axis=1)
    upper = nxp.max(central, axis=1)
    spread = upper - lower
    fence_lower = lower - factor * spread
    fence_upper = upper + factor * spread
    outside = (values < fence_lower[:, None]) | (values > fence_upper[:, None])
    outlying = nxp.any(outside, axis=0)
    regular = values[:, ~outlying]
    median_index = int(nxp.argmax(scores))
    return BoxplotResult(
        values=values,
        t=points,
        depth=scores,
        median_index=median_index,
        median=values[:, median_index],
        outliers=nxp.nonzero(outlying)[0],
        central_lower=lower,
        central_upper=upper,
        fence_lower=fence_lower,
        fence_upper=fence_upper,
        whisker_lower=nxp.min(regular, axis=1),
        whisker_upper=nxp.max(regular, axis=1),
    )


# --------------------------------------------------------------------------- #
# permutation tests
# --------------------------------------------------------------------------- #


def _generator(random_state: RandomState) -> Any:
    """Return an object drawing permutations through ``.permutation(n)``."""
    if random_state is None or (
        isinstance(random_state, int) and not isinstance(random_state, bool)
    ):
        return default_namespace().random.default_rng(random_state)
    if hasattr(random_state, "permutation"):
        return random_state
    try:
        seed = operator.index(random_state)
    except TypeError:
        raise TypeError(
            "random_state must be None, an integer seed or a numpy.random.Generator, "
            f"got {type(random_state).__name__}"
        ) from None
    return default_namespace().random.default_rng(seed)


def _check_test_args(n_perm: int, q: float) -> None:
    if n_perm < 1:
        raise ValueError(f"n_perm must be at least 1, got {n_perm}")
    if not 0.0 < q < 1.0:
        raise ValueError(f"q must lie strictly between 0 and 1, got {q}")


def _summarise(
    observed: NDArray, null: NDArray, q: float, points: NDArray | None
) -> PermutationTestResult:
    """Collapse pointwise statistics into a max-statistic permutation test."""
    nxp = default_namespace()
    statistic = float(nxp.max(observed))
    maxima = nxp.max(null, axis=1)
    return PermutationTestResult(
        statistic=statistic,
        null=maxima,
        pvalue=float(nxp.mean(maxima >= statistic)),
        critical_value=float(_quantile(maxima, 1.0 - q)),
        pointwise=observed,
        pointwise_null=null,
        pointwise_pvalue=nxp.mean(null >= observed[None, :], axis=0),
        pointwise_critical_value=_quantile(null, 1.0 - q),
        t=points,
    )


def _welch(first: NDArray, second: NDArray) -> Array:
    """Pointwise ``|mean_1 - mean_2| / sqrt(var_1 / n_1 + var_2 / n_2)``."""
    nxp = default_namespace()
    n1, n2 = first.shape[1], second.shape[1]
    m1 = nxp.mean(first, axis=1)
    m2 = nxp.mean(second, axis=1)
    v1 = nxp.sum((first - m1[:, None]) ** 2, axis=1) / (n1 - 1)
    v2 = nxp.sum((second - m2[:, None]) ** 2, axis=1) / (n2 - 1)
    return nxp.abs(m1 - m2) / nxp.sqrt(v1 / n1 + v2 / n2)


def t_test(
    x1: FData,
    x2: FData,
    *,
    n_perm: int = 200,
    q: float = 0.05,
    t: Any = None,
    random_state: RandomState = None,
) -> PermutationTestResult:
    """Permutation test for a difference between two groups of curves.

    Replaces R's ``tperm.fd``.  The pointwise statistic is Welch's

    ``T(t) = |xbar_1(t) - xbar_2(t)| / sqrt(s_1^2(t) / n_1 + s_2^2(t) / n_2)``

    and the test statistic is its maximum over the evaluation points.  Its null
    distribution is obtained by pooling the curves, drawing a random
    permutation, and assigning the first ``n_1`` to the first group.

    Parameters
    ----------
    x1, x2 : FData
        The two groups, each with at least two curves, on one domain.
    n_perm : int, optional
        Number of permutations.  Default ``200``.
    q : float, optional
        Upper-tail probability for the critical value: ``critical_value`` is
        the ``1 - q`` quantile of the null distribution.  Default ``0.05``.
    t : array_like, optional
        Evaluation points; 101 equally spaced points over the domain by default.
    random_state : int, numpy.random.Generator or None, optional
        Seed or generator for the permutations; each one is drawn with the
        generator's ``permutation(n)``.

    Returns
    -------
    PermutationTestResult
        Observed statistic, null distribution, p-value and pointwise detail.

    Raises
    ------
    ValueError
        If a group has fewer than two curves, the groups differ in domain or
        carry several variables, or ``n_perm`` or ``q`` is out of range.
    TypeError
        If ``random_state`` is not a seed or a generator.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.stats import t_test
    >>> rng = np.random.default_rng(0)
    >>> a = FData(rng.normal(size=(5, 8)), BSpline(n_basis=5))
    >>> b = FData(rng.normal(size=(5, 7)), BSpline(n_basis=5))
    >>> result = t_test(a, b, n_perm=100, random_state=0)
    >>> 0.0 <= result.pvalue <= 1.0
    True
    """
    _check_test_args(n_perm, q)
    if x1.n_vars != 1 or x2.n_vars != 1:
        raise ValueError("t_test needs curves with a single variable")
    if x1.n_curves < 2 or x2.n_curves < 2:
        raise ValueError("each group needs at least two curves")
    if not _same_domain(x1.domain, x2.domain):
        raise ValueError(f"the groups must share a domain, got {x1.domain} and {x2.domain}")
    nxp = default_namespace()
    points = _default_grid(x1.domain) if t is None else to_numpy(t)
    pooled = nxp.concat([to_numpy(x1(points)), to_numpy(x2(points))], axis=1)
    n1, n = x1.n_curves, pooled.shape[1]
    observed = _welch(pooled[:, :n1], pooled[:, n1:])
    rng = _generator(random_state)
    null = nxp.empty((n_perm, points.shape[0]), dtype=nxp.float64)
    for k in range(n_perm):
        order = nxp.asarray(rng.permutation(n))
        null[k] = _welch(pooled[:, order[:n1]], pooled[:, order[n1:]])
    return _summarise(observed, null, q, points)


# --------------------------------------------------------------------------- #
# the regression behind f_test
# --------------------------------------------------------------------------- #


def _per_covariate(value: Any, n_cov: int, name: str) -> list[Any]:
    """Broadcast one setting to every covariate, or check a per-covariate list.

    Only a list or a tuple is read as one entry per covariate; anything else
    (a number, a NumPy scalar, a basis, an operator, ``None``) is shared.
    """
    if not isinstance(value, (list, tuple)):
        return [value] * n_cov
    entries = list(value)
    if len(entries) != n_cov:
        raise ValueError(
            f"{name} must be a single value or one entry per covariate "
            f"({n_cov}), got {len(entries)}"
        )
    return entries


def _covariates(x: Any, n: int, domain: tuple[float, float] | None) -> list[FData | NDArray]:
    """Validate the covariates: scalar vectors of length ``n`` or ``n`` curves."""
    nxp = default_namespace()
    items = list(x) if isinstance(x, (list, tuple)) else [x]
    if not items:
        raise ValueError("f_test needs at least one covariate")
    out: list[FData | NDArray] = []
    for item in items:
        if isinstance(item, FData):
            if item.n_vars != 1:
                raise ValueError("functional covariates must have a single variable")
            if item.n_curves != n:
                raise ValueError(f"every covariate needs {n} observations, got {item.n_curves}")
            if domain is not None and not _same_domain(item.domain, domain):
                raise ValueError(
                    f"functional covariates must share the response domain {domain}, "
                    f"got {item.domain}"
                )
            out.append(item)
            continue
        vector = nxp.asarray(to_numpy(item), dtype=nxp.float64)
        if len(vector.shape) != 1:
            raise ValueError(
                f"a scalar covariate must be one-dimensional, got shape {vector.shape}"
            )
        if vector.shape[0] != n:
            raise ValueError(f"every covariate needs {n} observations, got {vector.shape[0]}")
        out.append(vector)
    return out


def _penalty_block(bases: Sequence[Basis | None], lams: Sequence[float], ops: Sequence[Any]) -> Any:
    """Block-diagonal roughness penalty, zero for plain scalar coefficients."""
    nxp = default_namespace()
    sizes = [1 if b is None else b.n_basis for b in bases]
    total = sum(sizes)
    out = nxp.zeros((total, total), dtype=nxp.float64)
    start = 0
    for basis, lam, op, size in zip(bases, lams, ops, sizes, strict=True):
        if basis is not None and lam > 0.0:
            matrix = to_numpy(basis.penalty(op if isinstance(op, LDO) else int(op)))
            out[start : start + size, start : start + size] = lam * matrix
        start += size
    return out


def _solve_normal(matrix: NDArray, rhs: NDArray) -> Array:
    """Solve the symmetric normal equations, refusing a singular system."""
    nxp = default_namespace()
    spectrum = nxp.linalg.eigvalsh(matrix)
    top = float(nxp.max(nxp.abs(spectrum)))
    floor = matrix.shape[0] * nxp.finfo(nxp.float64).eps * top
    if top == 0.0 or float(nxp.min(spectrum)) <= floor:
        raise ValueError(
            "the regression's normal equations are singular: the covariates are "
            "collinear or the coefficients are not identified; add a penalty "
            "(lam > 0) or drop a covariate"
        )
    return nxp.asarray(_linalg.solve_spd(matrix, rhs))


def _functional_response(
    y: FData,
    covariates: list[FData | NDArray],
    bases: list[Basis],
    lams: list[float],
    ops: list[Any],
    points: NDArray,
    orders: list[NDArray],
) -> tuple[NDArray, NDArray]:
    """Pointwise F statistics of the concurrent model for each response order.

    The model is ``y_i(t) = sum_j x_ij(t) beta_j(t)`` with ``beta_j`` expanded in
    ``bases[j]`` and penalised by ``lams[j] * int (L_j beta_j)^2``; a scalar
    covariate is constant in ``t``.  The normal equations

    ``sum_k [sum_i int x_ij x_ik phi_j phi_k'] b_k + lam_j R_j b_j
    = sum_i int x_ij y_i phi_j``

    are integrated by Gauss--Legendre quadrature on the joint break points of
    every basis involved (exact for splines), and the left-hand side, which
    does not depend on the order of the responses, is shared by all of them.
    """
    nxp = default_namespace()
    functional = [c.basis for c in covariates if isinstance(c, FData)]
    nodes, weights = _quadrature(y.basis, *bases, *functional)

    def at(values: FData | NDArray, grid: NDArray) -> NDArray:
        if isinstance(values, FData):
            return to_numpy(values(grid))
        return values[None, :]

    x_nodes = [at(c, nodes) for c in covariates]
    x_points = [at(c, points) for c in covariates]
    phi_nodes = [to_numpy(b(nodes)) for b in bases]
    phi_points = [to_numpy(b(points)) for b in bases]
    blocks = []
    for j in range(len(bases)):
        row = []
        for k in range(len(bases)):
            density = weights * nxp.sum(x_nodes[j] * x_nodes[k], axis=1)
            row.append(
                nxp.matmul(nxp.matrix_transpose(phi_nodes[j]), density[:, None] * phi_nodes[k])
            )
        blocks.append(nxp.concat(row, axis=1))
    normal = nxp.concat(blocks, axis=0) + _penalty_block(bases, lams, ops)
    y_nodes = to_numpy(y(nodes))
    y_points = to_numpy(y(points))
    rhs = nxp.stack(
        [
            nxp.concat(
                [
                    nxp.matmul(
                        nxp.matrix_transpose(phi),
                        weights * nxp.sum(xj * y_nodes[:, order], axis=1),
                    )
                    for phi, xj in zip(phi_nodes, x_nodes, strict=True)
                ]
            )
            for order in orders
        ],
        axis=1,
    )
    solution = _solve_normal(normal, rhs)
    sizes = [b.n_basis for b in bases]
    stats = nxp.empty((len(orders), points.shape[0]), dtype=nxp.float64)
    for column, order in enumerate(orders):
        fitted = nxp.zeros(y_points.shape, dtype=nxp.float64)
        start = 0
        for phi, xj, size in zip(phi_points, x_points, sizes, strict=True):
            beta = nxp.matmul(phi, solution[start : start + size, column])
            fitted = fitted + xj * beta[:, None]
            start += size
        stats[column] = _f_ratio(y_points[:, order], fitted)
    return stats[0], stats[1:]


def _f_ratio(observed: NDArray, fitted: NDArray) -> Array:
    """Pointwise ``var(fitted) / mean((observed - fitted)^2)`` across axis 1."""
    nxp = default_namespace()
    n = fitted.shape[1]
    centred = fitted - nxp.mean(fitted, axis=1, keepdims=True)
    explained = nxp.sum(centred * centred, axis=1) / (n - 1)
    residual = observed - fitted
    return explained / nxp.mean(residual * residual, axis=1)


def _scalar_response(
    y: NDArray,
    covariates: list[FData | NDArray],
    bases: list[Basis | None],
    lams: list[float],
    ops: list[Any],
    orders: list[NDArray],
) -> tuple[NDArray, NDArray]:
    """F statistics of ``y_i = sum_j int x_ij beta_j + sum_j x_ij b_j``."""
    nxp = default_namespace()
    columns: list[NDArray] = []
    for covariate, basis in zip(covariates, bases, strict=True):
        if isinstance(covariate, FData):
            target = covariate.basis if basis is None else basis
            columns.append(to_numpy(inprod(covariate, target)))
        else:
            columns.append(covariate[:, None])
    design = nxp.concat(columns, axis=1)
    normal = nxp.matmul(nxp.matrix_transpose(design), design) + _penalty_block(bases, lams, ops)
    responses = nxp.stack([y[order] for order in orders], axis=1)
    solution = _solve_normal(normal, nxp.matmul(nxp.matrix_transpose(design), responses))
    fitted = nxp.matmul(design, solution)
    stats = _f_ratio(nxp.matrix_transpose(responses), nxp.matrix_transpose(fitted))
    return stats[:1], stats[1:, None]


def f_test(
    y: FData | Any,
    x: Any,
    *,
    basis: Basis | Sequence[Basis | None] | None = None,
    lam: float | Sequence[float] = 0.0,
    penalty: int | LDO | Sequence[int | LDO] = 2,
    n_perm: int = 200,
    q: float = 0.05,
    t: Any = None,
    random_state: RandomState = None,
) -> PermutationTestResult:
    r"""Permutation F test of no effect in a functional linear regression.

    Replaces R's ``Fperm.fd``.  The regression is fitted by penalised least
    squares, then

    ``F(t) = var_i(yhat_i(t)) / mean_i((y_i(t) - yhat_i(t))^2)``

    (the variance with denominator ``n - 1``) is evaluated pointwise and the
    test statistic is its maximum.  The null distribution refits the model
    after randomly permuting the responses against fixed covariates.  No
    intercept is added: include a vector of ones in ``x`` for one, as in R.

    Two models are supported:

    * functional response (:class:`FData`): the concurrent model
      ``y_i(t) = sum_j x_ij(t) beta_j(t)``, where a scalar covariate is
      constant in ``t``, every ``beta_j`` is a curve in ``basis[j]`` (default
      ``y.basis``) and the penalty is ``lam[j] * int (L_j beta_j)^2``;
    * scalar response (1-D array): ``y_i = sum_j int x_ij(t) beta_j(t) dt`` for
      functional covariates (``beta_j`` in ``basis[j]``, default the
      covariate's basis) plus ordinary coefficients for scalar covariates, whose
      ``basis``, ``lam`` and ``penalty`` entries are ignored.  ``F`` is then a
      single number.

    Parameters
    ----------
    y : FData or array_like
        Response: ``n`` curves, or ``n`` numbers.
    x : array_like, FData or sequence of them
        Covariates, each a length-``n`` vector or an :class:`FData` of ``n``
        curves on the response domain.  A single covariate need not be wrapped
        in a list.
    basis : Basis or list of Basis, optional
        Basis of each coefficient function: one shared value, or a list (or
        tuple) with one entry per covariate, ``None`` meaning the default.
    lam : float or list of float, optional
        Smoothing parameter of each coefficient function, shared or one per
        covariate.  Default ``0.0``.
    penalty : int, LDO or list of them, optional
        Roughness operator of each coefficient function, shared or one per
        covariate.  Default ``2``.
    n_perm : int, optional
        Number of permutations.  Default ``200``.
    q : float, optional
        Upper-tail probability: ``critical_value`` is the ``1 - q`` quantile of
        the null distribution.  Default ``0.05``.
    t : array_like, optional
        Evaluation points for a functional response; 101 equally spaced points
        over the domain by default.
    random_state : int, numpy.random.Generator or None, optional
        Seed or generator; each permutation is ``permutation(n)``.

    Returns
    -------
    PermutationTestResult
        Observed statistic, null distribution, p-value and pointwise detail
        (``t`` is ``None`` for a scalar response).

    Raises
    ------
    ValueError
        If the covariates, bases or settings do not match the response, the
        normal equations are singular, or ``n_perm`` or ``q`` is out of range.
    TypeError
        If ``random_state`` is not a seed or a generator.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import FData, Fourier
    >>> from fabel.stats import f_test
    >>> rng = np.random.default_rng(0)
    >>> group = np.repeat([0.0, 1.0], 6)
    >>> y = FData(rng.normal(size=(5, 12)) + 2.0 * group, Fourier(n_basis=5))
    >>> result = f_test(y, [np.ones(12), group], lam=1e-4, n_perm=50, random_state=0)
    >>> result.pvalue < 0.05
    True
    """
    _check_test_args(n_perm, q)
    nxp = default_namespace()
    if isinstance(y, FData):
        if y.n_vars != 1:
            raise ValueError("f_test needs a response with a single variable")
        n = y.n_curves
        domain: tuple[float, float] | None = y.domain
    else:
        response = nxp.asarray(to_numpy(y), dtype=nxp.float64)
        if len(response.shape) != 1:
            raise ValueError(f"a scalar response must be one-dimensional, got {response.shape}")
        n = response.shape[0]
        domain = None
    covariates = _covariates(x, n, domain)
    n_cov = len(covariates)
    chosen = _per_covariate(basis, n_cov, "basis")
    lams = [float(v) for v in _per_covariate(lam, n_cov, "lam")]
    ops = _per_covariate(penalty, n_cov, "penalty")
    if any(v < 0.0 for v in lams):
        raise ValueError(f"lam must be non-negative, got {lams}")
    rng = _generator(random_state)
    orders = [nxp.arange(n)] + [nxp.asarray(rng.permutation(n)) for _ in range(n_perm)]
    if isinstance(y, FData):
        bases: list[Basis] = [y.basis if b is None else b for b in chosen]
        for b in bases:
            if not _same_domain(b.domain, y.domain):
                raise ValueError(
                    f"coefficient bases must share the response domain {y.domain}, got {b.domain}"
                )
        points = _default_grid(y.domain) if t is None else to_numpy(t)
        observed, null = _functional_response(y, covariates, bases, lams, ops, points, orders)
        return _summarise(observed, null, q, points)
    scalar_bases: list[Basis | None] = []
    for covariate, b in zip(covariates, chosen, strict=True):
        if isinstance(covariate, FData):
            target = covariate.basis if b is None else b
            if not _same_domain(target.domain, covariate.domain):
                raise ValueError(
                    f"coefficient bases must share the covariate domain {covariate.domain}, "
                    f"got {target.domain}"
                )
            scalar_bases.append(target)
        else:
            scalar_bases.append(None)
    observed, null = _scalar_response(response, covariates, scalar_bases, lams, ops, orders)
    return _summarise(observed, null, q, None)
