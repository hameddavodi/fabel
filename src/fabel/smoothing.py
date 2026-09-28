r"""Penalised smoothing: turn discrete observations into functional data.

A single entry point, :func:`smooth`, replaces the twelve R ``fda`` smoothing
functions (``smooth.basis``, ``smooth.basisPar``, ``Data2fd``,
``smooth.monotone``, ``smooth.pos``, ``smooth.morph``, ``lambda2df``,
``lambda2gcv``, ``df2lambda`` and friends).  The unconstrained problem is the
penalised least-squares fit of Ramsay & Silverman chapter 5,

.. math::

    (\Phi^{T} W \Phi + \lambda R)\,c = \Phi^{T} W y ,
    \qquad R_{ij} = \int (L\varphi_i)(L\varphi_j)\,dt ,

whose hat matrix ``H = Φ (ΦᵀWΦ + λR)⁻¹ ΦᵀW`` gives the equivalent degrees of
freedom ``df = tr H`` and the generalised cross-validation criterion
``gcv_i = (SSE_i / n) / (1 - df/n)²`` for curve ``i``.

Constrained fits reuse the same penalty on a latent function ``W``:
``constraint="positive"`` fits ``x = exp W``, ``constraint="monotone"`` fits
``x(t) = β₀ + β₁ ∫ₐᵗ exp W(u) du``, and ``constraint="morph"`` is the monotone
fit with ``β`` pinned so that the fitted curve maps the domain onto itself.
Their derivatives of any order are exact: with ``w_j = D^j W``, Faà di Bruno's
formula gives ``Dⁿ exp W = exp(W) Bₙ(w_1, …, w_n)`` for the complete Bell
polynomial ``Bₙ`` (``B₀ = 1``, ``B_{n+1} = Σ_{i=0}^{n} C(n, i) B_{n-i} w_{i+1}``),
and ``Dⁿ x = β₁ D^{n-1} exp W`` for a monotone fit.

Examples
--------
>>> import numpy as np
>>> from fabel import BSpline
>>> from fabel.smoothing import smooth
>>> t = np.linspace(0.0, 1.0, 21)
>>> y = np.sin(2 * np.pi * t)
>>> result = smooth(y, t, basis=BSpline(domain=(0.0, 1.0), n_basis=12), lam=1e-8)
>>> bool(np.max(np.abs(result.fd(t)[:, 0] - y)) < 1e-3)
True
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import comb, inf, isfinite, sqrt
from types import ModuleType
from typing import Any

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted, validate_data

from fabel import _linalg
from fabel._backend import array_namespace, asarray, default_namespace, to_numpy
from fabel._operator import LDO
from fabel.basis import Basis, BSpline
from fabel.core import FData

__all__ = [
    "SmoothResult",
    "Smoother",
    "df_to_lambda",
    "gcv_curve",
    "lambda_to_df",
    "smooth",
]

Array = Any

#: Bounds of the ``log10(lambda)`` search interval used by ``lam="gcv"``.
LOG10_LAMBDA_RANGE = (-8.0, 8.0)

#: Widest ``log10(lambda)`` bracket searched when solving ``df(lambda) = target``.
_LOG10_LAMBDA_BOUNDS = (-15.0, 15.0)

_GOLDEN_RATIO = 0.5 * (sqrt(5.0) - 1.0)
_COARSE_GRID_STEP = 0.5
_GOLDEN_TOL = 1e-8
_BISECTION_ITERATIONS = 200

#: Auto-basis size: ``min(n_t, _AUTO_BASIS_CAP) + _AUTO_BASIS_EXTRA`` cubic
#: B-spline functions on the observed range.  Deterministic and independent of
#: the data values, so two calls on the same design always agree.
_AUTO_BASIS_CAP = 40
_AUTO_BASIS_EXTRA = 2
_AUTO_BASIS_ORDER = 4

#: Gauss-Legendre nodes per panel for the monotone fit's cumulative integrals.
_MONOTONE_QUAD_DEGREE = 12
_MAX_GAUSS_NEWTON_STEPS = 100
_MAX_HALVINGS = 20
_GAUSS_NEWTON_TOL = 1e-12

_CONSTRAINTS = ("positive", "monotone", "morph")

#: Relative gap ``(n - df) / n`` below which a fit counts as interpolating and
#: its GCV score as undefined.  A genuine near-interpolant keeps a gap of
#: 1.7e-8 at ``lambda = 1e-12`` (8 splines on 8 points); an exact one 4e-16.
_INTERPOLATION_TOL = 1e-10


# --------------------------------------------------------------------------- #
# result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SmoothResult:
    """Outcome of a :func:`smooth` call.

    Attributes
    ----------
    fd : FData
        The fitted curves.  For a constrained fit this holds the *latent*
        function ``W``; call the result itself (``result(t, deriv)``) to
        evaluate the constrained curve ``exp W`` or ``β₀ + β₁ ∫ exp W`` and
        its derivatives of any order.
    df : float
        Equivalent degrees of freedom, ``tr H``.  For a constrained fit this is
        the trace of the hat matrix of the final Gauss-Newton linearisation, and
        for per-curve irregular designs it is the sum over the independent fits.
    gcv : array
        Generalised cross-validation score per curve, shaped like the curve axes
        of ``y`` (a scalar for a single unnamed curve, ``(n_curves,)`` for a
        matrix, ``(n_curves, n_vars)`` for a multivariate one).
    sse : float
        Residual sum of squares of the fit, summed over every point and curve.
        Unweighted even when ``weights`` is given, matching R's ``smooth.basis``.
    penalty_matrix : array
        The roughness penalty ``R``.
    lam : float
        The smoothing parameter actually used -- the value selected by GCV or by
        the degrees-of-freedom target when one of those was requested.
    y2c_map : array or tuple of arrays
        The linear map ``(ΦᵀWΦ + λR)⁻¹ ΦᵀW`` from data to coefficients, of shape
        ``(n_basis, n_t)``.  A tuple with one entry per curve for an irregular
        design, and ``None`` for a constrained fit, which is not linear in ``y``.
    beta : array or None
        ``(2, n_curves)`` location/scale pair of a monotone or morph fit; ``None``
        otherwise.
    constraint : str or None
        The constraint that was applied, if any.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.smoothing import smooth
    >>> t = np.linspace(0.0, 1.0, 15)
    >>> result = smooth(t**2, t, lam=1e-8)
    >>> result.constraint is None and result.beta is None
    True
    """

    fd: FData
    df: float
    gcv: Array
    sse: float
    penalty_matrix: Array
    lam: float
    y2c_map: Array
    beta: Array | None = None
    constraint: str | None = None

    def __call__(self, t: Any, deriv: int = 0) -> Array:
        """Evaluate the fitted curves or one of their derivatives.

        For a constrained fit the derivatives are exact to rounding error, for
        any order: ``Dⁿ exp W = exp(W) Bₙ(DW, …, DⁿW)`` with the complete Bell
        polynomial ``Bₙ`` (Faà di Bruno's formula), and a monotone curve has
        ``Dⁿ x = β₁ D^{n-1} exp W`` for ``n >= 1``.  The value itself
        (``deriv=0``) of a monotone curve is ``β₀ + β₁ ∫ₐᵗ exp W``, integrated by
        Gauss-Legendre quadrature between the basis break points.

        Parameters
        ----------
        t : array
            Points at which to evaluate.
        deriv : int, optional
            Derivative order, ``0`` (the default) or more.

        Returns
        -------
        array
            Fitted values with the shape of ``fd(t)``.

        Raises
        ------
        ValueError
            If ``deriv`` is negative, or a monotone fit carries no ``beta``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.smoothing import smooth
        >>> t = np.linspace(0.0, 1.0, 40)
        >>> fit = smooth(np.exp(2 * t), t, constraint="positive", lam=1e-8)
        >>> np.round(fit(np.array([0.5]), 3)[:, 0] / np.exp(1.0), 2)  # D³ e^{2t} = 8 e^{2t}
        array([8.])
        """
        if self.constraint is None:
            return self.fd(t, deriv)
        if deriv < 0:
            raise ValueError(f"deriv must be non-negative, got {deriv}")
        xp = array_namespace(self.fd.coefs)
        points = asarray(t, xp)
        if self.constraint == "positive":
            return _exp_derivative(self.fd, points, deriv, xp)
        if self.beta is None:
            raise ValueError("a monotone fit must carry its beta coefficients")
        beta: Array = self.beta
        if deriv:
            return _broadcast_curve(_exp_derivative(self.fd, points, deriv - 1, xp), beta[1], xp)
        values, _ = _cumulative_exp(self.fd.basis, self.fd.coefs, points, derivative=False)
        return _broadcast_curve(values, beta[1], xp) + _broadcast_curve(
            xp.ones_like(values), beta[0], xp
        )


def _exp_derivative(latent: FData, points: Array, order: int, xp: ModuleType) -> Array:
    """Return ``Dⁿ exp W`` at ``points`` for ``n = order``, exactly.

    Faà di Bruno's formula for the exponential: ``Dⁿ exp W = exp(W) Bₙ`` where
    the complete Bell polynomials in ``w_j = D^j W`` satisfy ``B₀ = 1`` and
    ``B_{k+1} = Σ_{i=0}^{k} C(k, i) B_{k-i} w_{i+1}``.
    """
    values = xp.exp(latent(points))
    if order == 0:
        return values
    slopes = [latent(points, j) for j in range(1, order + 1)]
    bell = [xp.ones_like(values)]
    for k in range(order):
        total = xp.zeros_like(values)
        for i in range(k + 1):
            total = total + comb(k, i) * bell[k - i] * slopes[i]
        bell.append(total)
    return values * bell[order]


# --------------------------------------------------------------------------- #
# shape handling
# --------------------------------------------------------------------------- #


def _is_irregular(t: Any) -> bool:
    """Return ``True`` when ``t`` is a sequence of per-curve argument vectors."""
    if not isinstance(t, (list, tuple)):
        return False
    return len(t) > 0 and hasattr(t[0], "__len__")


def _flatten_curves(y: Array, xp: ModuleType) -> tuple[Array, tuple[int, ...]]:
    """Reshape ``y`` to ``(n_t, n_columns)`` and return the trailing curve shape."""
    curve_shape = tuple(y.shape[1:])
    return xp.reshape(y, (y.shape[0], -1)), curve_shape


def _broadcast_curve(values: Array, per_curve: Array, xp: ModuleType) -> Array:
    """Multiply ``(n_t, ...)`` values by a per-curve vector broadcast over time."""
    if values.ndim == 1:
        return values * xp.reshape(per_curve, ())
    return values * per_curve


def _auto_basis(t: Array) -> BSpline:
    """Build the default cubic B-spline basis for an observation grid."""
    lower, upper = float(min(t)), float(max(t))
    if not upper > lower:
        raise ValueError("t must span a non-degenerate interval to build a default basis")
    n_points = len(t)
    n_basis = min(n_points, _AUTO_BASIS_CAP) + _AUTO_BASIS_EXTRA
    order = min(_AUTO_BASIS_ORDER, n_basis)
    return BSpline(domain=(lower, upper), n_basis=n_basis, order=order)


def _as_operator(penalty: int | LDO) -> LDO:
    """Coerce a penalty specification into a linear differential operator."""
    return penalty if isinstance(penalty, LDO) else LDO(int(penalty))


def _weight_vector(weights: Any, n_points: int, xp: ModuleType) -> Array:
    """Return the observation weights as a validated ``(n_t,)`` vector."""
    if weights is None:
        return xp.ones(n_points, dtype=xp.float64)
    w = asarray(weights, xp)
    if w.ndim != 1 or w.shape[0] != n_points:
        raise ValueError(f"weights must have shape ({n_points},), got {tuple(w.shape)}")
    return w


# --------------------------------------------------------------------------- #
# the linear fit
# --------------------------------------------------------------------------- #


def _normal_matrix(phi: Array, w: Array, xp: ModuleType) -> tuple[Array, Array]:
    """Return ``(ΦᵀWΦ, ΦᵀW)`` for a design matrix and weight vector."""
    weighted = w[:, None] * phi
    return xp.matmul(xp.matrix_transpose(phi), weighted), xp.matrix_transpose(weighted)


def _gcv_denominator(df: float, n_points: int) -> float:
    """Return ``(1 - df/n)²``, or ``0`` when the fit interpolates the data.

    At ``df = n`` the criterion is ``0 / 0``: the residuals and the denominator
    are both pure rounding, and their ratio is noise (a measured 0.27 for an
    exact interpolant).  R's ``smooth.basis`` returns no GCV there; Fabel
    reports ``inf`` once ``n - df`` is within :data:`_INTERPOLATION_TOL` of
    zero relative to ``n``.
    """
    if n_points - df <= _INTERPOLATION_TOL * n_points:
        return 0.0
    return (1.0 - df / n_points) ** 2


def _linear_fit(
    phi: Array, ymat: Array, w: Array, pen: Array, lam: float, xp: ModuleType
) -> tuple[Array, Array, float, Array, float]:
    """Solve the penalised normal equations.

    Returns ``(coefs, y2c_map, df, gcv_per_column, sse)`` for ``ymat`` shaped
    ``(n_t, n_columns)``.
    """
    gram, phi_w = _normal_matrix(phi, w, xp)
    y2c = _linalg.solve_spd(gram + lam * pen, phi_w)
    coefs = xp.matmul(y2c, ymat)
    fitted = xp.matmul(phi, coefs)
    n_points = phi.shape[0]
    df = float(xp.sum(phi * xp.matrix_transpose(y2c)))
    resid = ymat - fitted
    sse_columns = xp.sum(resid * resid, axis=0)
    denom = _gcv_denominator(df, n_points)
    gcv = (sse_columns / n_points) / denom if denom > 0 else xp.full_like(sse_columns, inf)
    return coefs, y2c, df, gcv, float(to_numpy(xp.sum(sse_columns)))


# --------------------------------------------------------------------------- #
# lambda selection -- one factorisation for the whole search
# --------------------------------------------------------------------------- #


class _Pencil:
    """Diagonalised smoothing pencil, giving ``df`` and ``gcv`` at O(K) per lambda.

    With ``μ, V`` from :func:`fabel._linalg.pencil_eigh` applied to
    ``(ΦᵀWΦ, ΦᵀWΦ + R)`` the inverse ``(ΦᵀWΦ + λR)⁻¹`` is
    ``V diag((μ + λ(1-μ))⁻¹) Vᵀ``, so a whole lambda grid reuses a single
    eigendecomposition.
    """

    def __init__(self, phi: Array, w: Array, pen: Array, ymat: Array | None) -> None:
        xp = default_namespace()
        gram, phi_w = _normal_matrix(phi, w, xp)
        mu, vec = _linalg.pencil_eigh(gram, gram + pen)
        self._mu = xp.clip(mu, 0.0, 1.0)
        self._n = int(phi.shape[0])
        self._projected = xp.matmul(phi, vec)
        self._rhs = None if ymat is None else xp.matmul(xp.matrix_transpose(vec), phi_w @ ymat)
        self._xp = xp

    def _scaling(self, lam: float) -> Array:
        """Return the diagonal ``(μ + λ(1-μ))⁻¹``."""
        denom = self._mu + lam * (1.0 - self._mu)
        return 1.0 / self._xp.clip(denom, 1e-300, None)

    def df(self, lam: float) -> float:
        """Equivalent degrees of freedom at ``lam``."""
        return float(self._xp.sum(self._mu * self._scaling(lam)))

    def gcv(self, lam: float) -> Array:
        """Per-column GCV score at ``lam``."""
        xp = self._xp
        scale = self._scaling(lam)
        df = float(xp.sum(self._mu * scale))
        resid = self._ymat - xp.matmul(self._projected, scale[:, None] * self._rhs)
        sse = xp.sum(resid * resid, axis=0)
        denom = _gcv_denominator(df, self._n)
        if denom <= 0.0:
            return xp.full_like(sse, inf)
        return (sse / self._n) / denom

    def attach(self, ymat: Array) -> None:
        """Record the response matrix used by :meth:`gcv`."""
        self._ymat = ymat


def _pencil_for(phi: Array, w: Array, pen: Array, ymat: Array | None) -> _Pencil:
    """Build a :class:`_Pencil` from NumPy copies of the fit inputs."""
    target = None if ymat is None else to_numpy(ymat)
    pencil = _Pencil(to_numpy(phi), to_numpy(w), to_numpy(pen), target)
    if ymat is not None:
        pencil.attach(to_numpy(ymat))
    return pencil


def _mean_gcv(pencil: _Pencil, log10_lam: float) -> float:
    """Mean GCV over curves at ``10 ** log10_lam``."""
    value = float(pencil._xp.mean(pencil.gcv(10.0**log10_lam)))
    return value if isfinite(value) else inf


def _minimise_gcv(pencil: _Pencil) -> float:
    """Minimise mean GCV over ``log10(lambda)``: coarse grid, then golden section."""
    lower, upper = LOG10_LAMBDA_RANGE
    n_grid = round((upper - lower) / _COARSE_GRID_STEP) + 1
    grid = [lower + i * _COARSE_GRID_STEP for i in range(n_grid)]
    scores = [_mean_gcv(pencil, g) for g in grid]
    best = min(range(n_grid), key=scores.__getitem__)
    left = grid[max(best - 1, 0)]
    right = grid[min(best + 1, n_grid - 1)]
    return _golden_section(lambda g: _mean_gcv(pencil, g), left, right)


def _golden_section(objective: Any, left: float, right: float) -> float:
    """Golden-section minimiser of a unimodal objective on ``[left, right]``."""
    if right - left < _GOLDEN_TOL:
        return 0.5 * (left + right)
    a, b = left, right
    c = b - _GOLDEN_RATIO * (b - a)
    d = a + _GOLDEN_RATIO * (b - a)
    fc, fd = objective(c), objective(d)
    while b - a > _GOLDEN_TOL:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - _GOLDEN_RATIO * (b - a)
            fc = objective(c)
        else:
            a, c, fc = c, d, fd
            d = a + _GOLDEN_RATIO * (b - a)
            fd = objective(d)
    return 0.5 * (a + b)


def _lambda_for_df(pencil: _Pencil, target: float) -> float:
    """Invert the (strictly decreasing) map ``lambda -> df`` by bisection."""
    low, high = _LOG10_LAMBDA_BOUNDS
    if pencil.df(10.0**low) < target:
        raise ValueError(f"target df {target} exceeds the {pencil.df(10.0**low):.6g} available")
    if pencil.df(10.0**high) > target:
        raise ValueError(f"target df {target} is below the achievable minimum")
    for _ in range(_BISECTION_ITERATIONS):
        mid = 0.5 * (low + high)
        if pencil.df(10.0**mid) > target:
            low = mid
        else:
            high = mid
    return float(10.0 ** (0.5 * (low + high)))


def _resolve_lambda(
    lam: float | str, df: float | None, phi: Array, w: Array, pen: Array, ymat: Array
) -> float:
    """Return the smoothing parameter requested by ``lam``/``df``."""
    if df is not None:
        return _lambda_for_df(_pencil_for(phi, w, pen, None), float(df))
    if isinstance(lam, str):
        if lam == "gcv":
            return float(10.0 ** _minimise_gcv(_pencil_for(phi, w, pen, ymat)))
        if lam.startswith("df="):
            return _lambda_for_df(_pencil_for(phi, w, pen, None), float(lam[3:]))
        raise ValueError(f"lam must be a number, 'gcv' or 'df=<value>', got {lam!r}")
    value = float(lam)
    if value < 0.0:
        raise ValueError(f"lam must be non-negative, got {value}")
    return value


# --------------------------------------------------------------------------- #
# public helpers
# --------------------------------------------------------------------------- #


def lambda_to_df(
    t: Any, basis: Basis, lam: float, *, penalty: int | LDO = 2, weights: Any = None
) -> float:
    """Equivalent degrees of freedom of a smooth at a given ``lam``.

    Parameters
    ----------
    t : array
        Observation points.
    basis : Basis
        Expansion basis.
    lam : float
        Smoothing parameter, on the linear scale.
    penalty : int or LDO, optional
        Roughness operator; default the curvature penalty ``D²``.
    weights : array, optional
        Observation weights.

    Returns
    -------
    float
        ``tr H``, the trace of the hat matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.smoothing import lambda_to_df
    >>> t = np.linspace(0.0, 1.0, 30)
    >>> round(lambda_to_df(t, BSpline(domain=(0.0, 1.0), n_basis=10), 0.0), 6)
    10.0
    """
    xp = default_namespace()
    points = asarray(t, xp)
    phi = basis(points)
    w = _weight_vector(weights, phi.shape[0], xp)
    return _pencil_for(phi, w, basis.penalty(_as_operator(penalty), xp=xp), None).df(float(lam))


def df_to_lambda(
    t: Any, basis: Basis, df: float, *, penalty: int | LDO = 2, weights: Any = None
) -> float:
    """Smoothing parameter whose smooth has ``df`` degrees of freedom.

    Parameters
    ----------
    t : array
        Observation points.
    basis : Basis
        Expansion basis.
    df : float
        Target equivalent degrees of freedom.
    penalty : int or LDO, optional
        Roughness operator; default the curvature penalty ``D²``.
    weights : array, optional
        Observation weights.

    Returns
    -------
    float
        The ``lam`` solving ``lambda_to_df(t, basis, lam) == df``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.smoothing import df_to_lambda, lambda_to_df
    >>> t = np.linspace(0.0, 1.0, 40)
    >>> basis = BSpline(domain=(0.0, 1.0), n_basis=12)
    >>> lam = df_to_lambda(t, basis, 6.0)
    >>> round(lambda_to_df(t, basis, lam), 8)
    6.0
    """
    xp = default_namespace()
    points = asarray(t, xp)
    phi = basis(points)
    w = _weight_vector(weights, phi.shape[0], xp)
    pen = basis.penalty(_as_operator(penalty), xp=xp)
    return _lambda_for_df(_pencil_for(phi, w, pen, None), float(df))


def gcv_curve(
    y: Any,
    t: Any,
    basis: Basis,
    lambdas: Any,
    *,
    penalty: int | LDO = 2,
    weights: Any = None,
) -> Array:
    """GCV score of every curve over a grid of smoothing parameters.

    One eigendecomposition serves the whole grid, so the cost is
    ``O(n_lambda · n_t · n_basis)`` rather than a factorisation per lambda.

    Parameters
    ----------
    y : array
        Observations of shape ``(n_t,)``, ``(n_t, n_curves)`` or
        ``(n_t, n_curves, n_vars)``.
    t : array
        Observation points, shape ``(n_t,)``.
    basis : Basis
        Expansion basis.
    lambdas : array
        Smoothing parameters on the linear scale.
    penalty : int or LDO, optional
        Roughness operator; default the curvature penalty ``D²``.
    weights : array, optional
        Observation weights.

    Returns
    -------
    array
        Shape ``(n_lambda, *curve_shape)``: one GCV score per lambda per curve.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.smoothing import gcv_curve
    >>> t = np.linspace(0.0, 1.0, 40)
    >>> y = np.sin(2 * np.pi * t)
    >>> scores = gcv_curve(y, t, BSpline(domain=(0.0, 1.0), n_basis=10), [1e-6, 1e-2])
    >>> scores.shape
    (2,)
    """
    xp = default_namespace()
    points = asarray(t, xp)
    values = asarray(y, xp)
    phi = basis(points)
    ymat, curve_shape = _flatten_curves(values, xp)
    w = _weight_vector(weights, phi.shape[0], xp)
    pencil = _pencil_for(phi, w, basis.penalty(_as_operator(penalty), xp=xp), ymat)
    grid = asarray(lambdas, xp)
    scores = xp.stack([pencil.gcv(float(lam)) for lam in grid])
    return xp.reshape(scores, (scores.shape[0], *curve_shape))


# --------------------------------------------------------------------------- #
# constrained fits
# --------------------------------------------------------------------------- #


def _integration_panels(basis: Basis, t: Array, xp: ModuleType) -> Array:
    """Panel boundaries for the cumulative integral: basis breaks plus the data."""
    lower = basis.domain[0]
    edges = xp.concat(
        [
            asarray(list(basis._natural_breaks()), xp),
            xp.reshape(asarray(t, xp), (-1,)),
            asarray([lower], xp),
        ]
    )
    edges = xp.unique_values(edges)
    upper = float(xp.max(asarray(t, xp)))
    return edges[(edges >= lower) & (edges <= upper)]


def _cumulative_exp(
    basis: Basis, coefs: Array, t: Array, *, derivative: bool
) -> tuple[Array, Array | None]:
    """Cumulative integrals ``h(t) = ∫ₐᵗ exp(W)`` and, optionally, ``∂h/∂c``.

    ``coefs`` may be a single ``(n_basis,)`` vector or a ``(n_basis, n_curves)``
    matrix; ``h`` comes back with the matching shape ``(n_t,)`` / ``(n_t, n_curves)``.
    The Jacobian is only defined for a single curve and is ``(n_t, n_basis)``.
    """
    xp = array_namespace(coefs)
    points = asarray(t, xp)
    panels = _integration_panels(basis, points, xp)
    nodes, weights = _linalg.composite_gauss_legendre(to_numpy(panels), _MONOTONE_QUAD_DEGREE)
    n_panels = panels.shape[0] - 1
    quad_nodes = asarray(nodes, xp)
    quad_weights = xp.reshape(asarray(weights, xp), (n_panels, _MONOTONE_QUAD_DEGREE))
    design = basis(quad_nodes)
    latent = xp.exp(xp.matmul(design, coefs))
    blocks = xp.reshape(latent, (n_panels, _MONOTONE_QUAD_DEGREE, *latent.shape[1:]))
    if latent.ndim == 1:
        panel_totals = xp.sum(blocks * quad_weights, axis=1)
    else:
        panel_totals = xp.sum(blocks * quad_weights[:, :, None], axis=1)
    zero = xp.zeros((1, *panel_totals.shape[1:]), dtype=panel_totals.dtype)
    cumulative = xp.concat([zero, xp.cumulative_sum(panel_totals, axis=0)])
    index = xp.searchsorted(panels, points)
    values = xp.take(cumulative, index, axis=0)
    if not derivative:
        return values, None
    weighted = latent * xp.reshape(quad_weights, (-1,))
    scaled = xp.reshape(weighted, (n_panels, _MONOTONE_QUAD_DEGREE, 1))
    per_panel = xp.sum(scaled * xp.reshape(design, (n_panels, _MONOTONE_QUAD_DEGREE, -1)), axis=1)
    zero_row = xp.zeros((1, per_panel.shape[1]), dtype=per_panel.dtype)
    jac = xp.concat([zero_row, xp.cumulative_sum(per_panel, axis=0)])
    return values, xp.take(jac, index, axis=0)


def _weighted_lstsq(design: Array, target: Array, w: Array, xp: ModuleType) -> Array:
    """Solve a small weighted least-squares problem ``(XᵀWX)β = XᵀWy``."""
    gram, cross = _normal_matrix(design, w, xp)
    return _linalg.solve_spd(gram, xp.matmul(cross, target))


def _gauss_newton(coefs: Array, evaluate: Any, step_of: Any, xp: ModuleType) -> tuple[Array, Any]:
    """Run damped Gauss-Newton from ``coefs`` until the step stops moving.

    ``evaluate(c)`` returns ``(criterion, parts)`` and ``step_of(c, parts)``
    returns the undamped Gauss-Newton step.  The step is halved until the
    criterion no longer increases; iteration stops when an accepted step is
    smaller than :data:`_GAUSS_NEWTON_TOL` relative to the coefficients, or when
    no damping of the step improves the criterion.
    """
    value, parts = evaluate(coefs)
    for _ in range(_MAX_GAUSS_NEWTON_STEPS):
        step = step_of(coefs, parts)
        scale = 1.0
        candidate = coefs
        improved = False
        for _ in range(_MAX_HALVINGS):
            candidate = coefs + scale * step
            new_value, new_parts = evaluate(candidate)
            if new_value <= value:
                improved = True
                break
            scale *= 0.5
        if not improved:
            break
        coefs, value, parts = candidate, new_value, new_parts
        moved = scale * float(xp.max(xp.abs(step)))
        if moved <= _GAUSS_NEWTON_TOL * (1.0 + float(xp.max(xp.abs(coefs)))):
            break
    return coefs, parts


def _fit_positive_curve(
    phi: Array, y: Array, w: Array, pen: Array, lam: float, xp: ModuleType
) -> tuple[Array, Array]:
    """Gauss-Newton fit of ``y ≈ exp(Φ c)`` with a roughness penalty on ``c``."""

    def evaluate(candidate: Array) -> tuple[float, Array]:
        latent = xp.exp(xp.matmul(phi, candidate))
        resid = y - latent
        penalty = lam * float(xp.matmul(candidate, xp.matmul(pen, candidate)))
        return float(xp.sum(w * resid * resid)) + penalty, latent

    def step_of(candidate: Array, latent: Array) -> Array:
        jac = latent[:, None] * phi
        hessian = xp.matmul(xp.matrix_transpose(jac), w[:, None] * jac) + lam * pen
        gradient = xp.matmul(xp.matrix_transpose(jac), w * (y - latent)) - lam * xp.matmul(
            pen, candidate
        )
        return _linalg.solve_spd(hessian, gradient)

    coefs, latent = _gauss_newton(xp.zeros(phi.shape[1], dtype=xp.float64), evaluate, step_of, xp)
    return coefs, latent[:, None] * phi


def _monotone_beta(
    values: Array, y: Array, w: Array, fixed: tuple[float, float] | None, xp: ModuleType
) -> Array:
    """Location/scale pair of a monotone fit: profiled out, or pinned for a morph."""
    if fixed is not None:
        total = float(values[-1])
        span = (fixed[1] - fixed[0]) / total if total > 0.0 else 1.0
        return asarray([fixed[0], span], xp)
    design = xp.stack([xp.ones_like(values), values], axis=1)
    return _weighted_lstsq(design, y, w, xp)


def _fit_monotone_curve(
    basis: Basis,
    y: Array,
    t: Array,
    w: Array,
    pen: Array,
    lam: float,
    fixed: tuple[float, float] | None,
    xp: ModuleType,
) -> tuple[Array, Array, Array]:
    """Gauss-Newton fit of ``y ≈ β₀ + β₁ ∫ₐᵗ exp(W)``.

    The criterion is exactly invariant under ``W -> W + s`` with
    ``β₁ -> β₁ e⁻ˢ`` whenever the penalty operator annihilates constants (every
    ``Dᵐ`` with ``m >= 1`` does), so the minimiser is a one-parameter family and
    the point reached depends on the iteration path.  The fitted curve, ``β₀``
    and the criterion itself are the same for every member.
    """

    def evaluate(candidate: Array) -> tuple[float, Any]:
        values, jac = _cumulative_exp(basis, candidate, t, derivative=True)
        beta = _monotone_beta(values, y, w, fixed, xp)
        resid = y - beta[0] - beta[1] * values
        penalty = lam * float(xp.matmul(candidate, xp.matmul(pen, candidate)))
        return float(xp.sum(w * resid * resid)) + penalty, (resid, beta, jac)

    def step_of(candidate: Array, parts: Any) -> Array:
        resid, beta, jac = parts
        design = beta[1] * jac
        hessian = xp.matmul(xp.matrix_transpose(design), w[:, None] * design) + lam * pen
        gradient = xp.matmul(xp.matrix_transpose(design), w * resid) - lam * xp.matmul(
            pen, candidate
        )
        return _linalg.solve_spd(hessian, gradient)

    coefs, parts = _gauss_newton(xp.zeros(basis.n_basis, dtype=xp.float64), evaluate, step_of, xp)
    _resid, beta, jac = parts
    return coefs, beta, beta[1] * jac


def _constrained_fit(
    basis: Basis,
    ymat: Array,
    t: Array,
    w: Array,
    pen: Array,
    lam: float,
    constraint: str,
    curve_shape: tuple[int, ...],
) -> SmoothResult:
    """Fit every column of ``ymat`` under a positivity or monotonicity constraint."""
    xp = default_namespace()
    n_points, n_columns = int(ymat.shape[0]), int(ymat.shape[1])
    phi = basis(t)
    fixed = basis.domain if constraint == "morph" else None
    # R's smooth.pos/smooth.monotone minimise mean(w r^2) + lam * c'Rc, where
    # smooth.basis minimises sum(w r^2) + lam * c'Rc.  Fabel keeps one criterion
    # (the sum) and rescales the constrained penalty, so `lam` means the same
    # thing here as in R.  Verified against the goldens: the stationarity
    # condition of R's reported coefficients holds at exactly n * lam.
    scaled = lam * n_points
    coef_columns: list[Array] = []
    beta_columns: list[Array] = []
    total_df = 0.0
    sse_columns: list[float] = []
    for column in range(n_columns):
        y = ymat[:, column]
        if constraint == "positive":
            coefs, design = _fit_positive_curve(phi, y, w, pen, scaled, xp)
            beta = None
            fitted = xp.exp(xp.matmul(phi, coefs))
        else:
            coefs, beta, design = _fit_monotone_curve(basis, y, t, w, pen, scaled, fixed, xp)
            values, _ = _cumulative_exp(basis, coefs, t, derivative=False)
            fitted = beta[0] + beta[1] * values
        gram, cross = _normal_matrix(design, w, xp)
        hat = _linalg.solve_spd(gram + scaled * pen, cross)
        total_df += float(xp.sum(design * xp.matrix_transpose(hat)))
        resid = y - fitted
        sse_columns.append(float(xp.sum(resid * resid)))
        coef_columns.append(coefs)
        if beta is not None:
            beta_columns.append(beta)
    coefs_all = xp.stack(coef_columns, axis=1)
    sse_array = asarray(sse_columns, xp)
    denom = _gcv_denominator(total_df / n_columns, n_points)
    gcv = (sse_array / n_points) / denom if denom > 0 else xp.full_like(sse_array, inf)
    return SmoothResult(
        fd=FData(xp.reshape(coefs_all, (basis.n_basis, *curve_shape)), basis),
        df=total_df,
        gcv=xp.reshape(gcv, curve_shape),
        sse=float(xp.sum(sse_array)),
        penalty_matrix=pen,
        lam=lam,
        y2c_map=None,
        beta=None if not beta_columns else xp.stack(beta_columns, axis=1),
        constraint=constraint,
    )


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #


def smooth(
    y: Any,
    t: Any,
    *,
    basis: Basis | None = None,
    lam: float | str = "gcv",
    penalty: int | LDO = 2,
    constraint: str | None = None,
    weights: Any = None,
    df: float | None = None,
) -> SmoothResult:
    r"""Fit smooth curves to discrete observations by penalised least squares.

    Parameters
    ----------
    y : array or sequence of arrays
        Observations of shape ``(n_t,)``, ``(n_t, n_curves)`` or
        ``(n_t, n_curves, n_vars)``.  For a per-curve irregular design, a
        sequence of 1-D arrays, one per curve.
    t : array or sequence of arrays
        Observation points shared by every curve, or one vector per curve for an
        irregular design (in which case ``y`` must be a matching sequence).
    basis : Basis, optional
        Expansion basis.  Defaults to a cubic B-spline on ``[min t, max t]`` with
        ``min(n_t, 40) + 2`` functions -- deterministic, and independent of the
        observed values.
    lam : float or str, optional
        Smoothing parameter.  A number is used directly; ``"gcv"`` (the default)
        minimises the mean GCV score over ``log10 λ ∈ [-8, 8]``; ``"df=12"``
        solves for the ``λ`` with 12 degrees of freedom.
    penalty : int or LDO, optional
        Roughness operator ``L``.  An integer means ``D^penalty``; the default
        ``2`` is the curvature penalty.
    constraint : {None, "positive", "monotone", "morph"}, optional
        Shape constraint.  ``"positive"`` fits ``x = exp W``, ``"monotone"`` fits
        ``x = β₀ + β₁ ∫ₐᵗ exp W``, and ``"morph"`` is the monotone fit with ``β``
        pinned so that ``x`` maps the domain onto itself.
    weights : array, optional
        Observation weights of shape ``(n_t,)``.
    df : float, optional
        Target degrees of freedom; equivalent to ``lam="df=<value>"`` and
        overriding ``lam``.

    Returns
    -------
    SmoothResult
        The fit, its degrees of freedom, GCV scores and data-to-coefficient map.

    Raises
    ------
    ValueError
        If the shapes of ``y`` and ``t`` disagree, if ``constraint`` is unknown,
        or if a constrained fit is asked for on an irregular design.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.smoothing import smooth
    >>> t = np.linspace(0.0, 1.0, 60)
    >>> y = np.exp(np.sin(4 * t))
    >>> fit = smooth(y, t, constraint="positive", lam=1e-6)
    >>> bool(np.min(fit(t)) > 0.0)
    True
    """
    if constraint is not None and constraint not in _CONSTRAINTS:
        raise ValueError(f"constraint must be one of {_CONSTRAINTS} or None, got {constraint!r}")
    if _is_irregular(t):
        if constraint is not None:
            raise ValueError("constrained fits require a single shared argument vector")
        return _smooth_irregular(y, t, basis, lam, penalty, weights, df)
    return _smooth_regular(y, t, basis, lam, penalty, constraint, weights, df)


def _smooth_regular(
    y: Any,
    t: Any,
    basis: Basis | None,
    lam: float | str,
    penalty: int | LDO,
    constraint: str | None,
    weights: Any,
    df: float | None,
) -> SmoothResult:
    """Fit curves observed on one shared argument vector."""
    xp = array_namespace(y, t)
    points = asarray(t, xp)
    values = asarray(y, xp)
    if points.ndim != 1:
        raise ValueError(f"t must be one-dimensional, got shape {tuple(points.shape)}")
    if values.shape[0] != points.shape[0]:
        raise ValueError(f"y has {values.shape[0]} rows but t has {points.shape[0]} points")
    used = _auto_basis(to_numpy(points)) if basis is None else basis
    operator = _as_operator(penalty)
    pen = used.penalty(operator, xp=xp)
    phi = used(points)
    w = _weight_vector(weights, points.shape[0], xp)
    ymat, curve_shape = _flatten_curves(values, xp)
    if constraint is not None:
        return _constrained_fit(
            used,
            to_numpy(ymat),
            to_numpy(points),
            to_numpy(w),
            to_numpy(pen),
            _resolve_lambda(lam, df, to_numpy(phi), to_numpy(w), to_numpy(pen), to_numpy(ymat)),
            constraint,
            curve_shape,
        )
    value = _resolve_lambda(lam, df, phi, w, pen, ymat)
    coefs, y2c, trace, gcv, sse = _linear_fit(phi, ymat, w, pen, value, xp)
    return SmoothResult(
        fd=FData(xp.reshape(coefs, (used.n_basis, *curve_shape)), used),
        df=trace,
        gcv=xp.reshape(gcv, curve_shape),
        sse=sse,
        penalty_matrix=pen,
        lam=value,
        y2c_map=y2c,
    )


def _smooth_irregular(
    y: Any,
    t: Sequence[Any],
    basis: Basis | None,
    lam: float | str,
    penalty: int | LDO,
    weights: Any,
    df: float | None,
) -> SmoothResult:
    """Fit one curve per argument vector, sharing a basis but not a design."""
    if not isinstance(y, (list, tuple)) or len(y) != len(t):
        raise ValueError("an irregular design needs y as a sequence matching t")
    xp = default_namespace()
    grids = [xp.reshape(asarray(points, xp), (-1,)) for points in t]
    responses = [xp.reshape(asarray(values, xp), (-1,)) for values in y]
    if basis is None:
        basis = _auto_basis(xp.concat(grids))
    operator = _as_operator(penalty)
    pen = basis.penalty(operator, xp=xp)
    if weights is None:
        weight_list: list[Any] = [None] * len(grids)
    elif isinstance(weights, (list, tuple)):
        weight_list = list(weights)
    else:
        weight_list = [weights] * len(grids)
    coef_columns: list[Array] = []
    maps: list[Array] = []
    gcv_values: list[float] = []
    total_df = 0.0
    total_sse = 0.0
    for grid, response, weight in zip(grids, responses, weight_list, strict=True):
        if grid.shape[0] != response.shape[0]:
            raise ValueError("each y entry must have as many points as its t entry")
        phi = basis(grid)
        w = _weight_vector(weight, grid.shape[0], xp)
        value = _resolve_lambda(lam, df, phi, w, pen, response[:, None])
        coefs, y2c, trace, gcv, sse = _linear_fit(phi, response[:, None], w, pen, value, xp)
        coef_columns.append(coefs[:, 0])
        maps.append(y2c)
        gcv_values.append(float(gcv[0]))
        total_df += trace
        total_sse += sse
        last_lambda = value
    coefs_all = xp.stack(coef_columns, axis=1)
    if len(coef_columns) == 1:
        coefs_all = coefs_all[:, 0]
    return SmoothResult(
        fd=FData(coefs_all, basis),
        df=total_df,
        gcv=asarray(gcv_values, xp) if len(gcv_values) > 1 else asarray(gcv_values[0], xp),
        sse=total_sse,
        penalty_matrix=pen,
        lam=last_lambda,
        y2c_map=maps[0] if len(maps) == 1 else tuple(maps),
    )


# --------------------------------------------------------------------------- #
# scikit-learn estimator
# --------------------------------------------------------------------------- #


class Smoother(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    """Penalised smoothing as a scikit-learn transformer.

    Each row of ``X`` is one curve sampled on the shared grid ``t``; ``transform``
    returns the basis coefficients, so a smoother can head a pipeline whose later
    steps work on the coefficient representation.

    Parameters
    ----------
    basis : Basis, optional
        Expansion basis.  ``None`` builds the default cubic B-spline of
        :func:`smooth` from the fitting grid.
    t : array, optional
        Sampling grid of length ``n_features``.  When omitted, ``fit`` uses
        ``0, 1, ..., n_features - 1``; it may also be supplied per call as the
        ``t`` fit parameter.
    lam : float or str, optional
        Smoothing parameter, as in :func:`smooth`.
    penalty : int or LDO, optional
        Roughness operator, as in :func:`smooth`.
    weights : array, optional
        Observation weights of length ``n_features``.

    Attributes
    ----------
    fd_ : FData
        The fitted curves, one per row of the training ``X``.
    result_ : SmoothResult
        The full fit, including ``df``, ``gcv`` and the data-to-coefficient map.
    t_ : numpy.ndarray
        The grid actually used.
    n_features_in_ : int
        Number of sampling points seen during ``fit``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.smoothing import Smoother
    >>> rng = np.random.default_rng(0)
    >>> t = np.linspace(0.0, 1.0, 25)
    >>> X = np.sin(2 * np.pi * t) + 0.05 * rng.standard_normal((6, 25))
    >>> coefs = Smoother(t=t, lam=1e-4).fit_transform(X)
    >>> coefs.shape[0]
    6
    """

    def __init__(
        self,
        basis: Basis | None = None,
        *,
        t: Any = None,
        lam: float | str = "gcv",
        penalty: int | LDO = 2,
        weights: Any = None,
    ) -> None:
        self.basis = basis
        self.t = t
        self.lam = lam
        self.penalty = penalty
        self.weights = weights

    def __sklearn_tags__(self) -> Any:
        """Declare that this transformer needs no ``y`` and takes dense 2-D input."""
        tags = super().__sklearn_tags__()
        tags.target_tags.required = False
        tags.input_tags.sparse = False
        return tags

    def _grid(self, n_features: int, override: Any) -> Array:
        """Return the sampling grid for ``n_features`` columns."""
        xp = default_namespace()
        chosen = self.t if override is None else override
        if chosen is None:
            return xp.arange(n_features, dtype=xp.float64)
        grid = xp.reshape(asarray(chosen, xp), (-1,))
        if grid.shape[0] != n_features:
            raise ValueError(f"t has {grid.shape[0]} points but X has {n_features} columns")
        return grid

    def fit(self, X: Any, y: Any = None, *, t: Any = None) -> Smoother:  # noqa: N803
        """Smooth every row of ``X`` and store the resulting curves.

        Parameters
        ----------
        X : array of shape (n_samples, n_features)
            One sampled curve per row.
        y : ignored
            Present for API compatibility.
        t : array, optional
            Sampling grid for this call, overriding the one given at construction.

        Returns
        -------
        Smoother
            The fitted estimator.
        """
        data = validate_data(self, X, ensure_min_features=2, ensure_min_samples=1)
        grid = self._grid(data.shape[1], t)
        self.t_ = grid
        self.result_ = smooth(
            data.T,
            grid,
            basis=self.basis,
            lam=self.lam,
            penalty=self.penalty,
            weights=self.weights,
        )
        self.fd_ = self.result_.fd
        return self

    def transform(self, X: Any) -> Array:  # noqa: N803
        """Return the basis coefficients of every row of ``X``.

        Parameters
        ----------
        X : array of shape (n_samples, n_features)
            Curves sampled on the grid seen during ``fit``.

        Returns
        -------
        numpy.ndarray
            Coefficients of shape ``(n_samples, n_basis)``.
        """
        check_is_fitted(self)
        data = validate_data(self, X, reset=False)
        xp = default_namespace()
        result = smooth(
            data.T,
            self.t_,
            basis=self.fd_.basis,
            lam=self.result_.lam,
            penalty=self.penalty,
            weights=self.weights,
        )
        return xp.matrix_transpose(xp.reshape(result.fd.coefs, (self.fd_.basis.n_basis, -1)))

    def get_feature_names_out(self, input_features: Any = None) -> Array:
        """Return the basis-function names produced by :meth:`transform`.

        Parameters
        ----------
        input_features : ignored
            Present for API compatibility; the output names come from the basis.

        Returns
        -------
        numpy.ndarray
            One name per basis function.
        """
        check_is_fitted(self)
        xp = default_namespace()
        return xp.asarray(self.fd_.basis.names, dtype=object)
