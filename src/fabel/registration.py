r"""Curve registration: remove phase variation by warping the time axis.

A single entry point, :func:`register`, replaces R ``fda``'s ``register.fd``
(continuous registration), ``landmarkreg`` (landmark registration),
``register.newfd`` and ``AmpPhaseDecomp`` (via
:meth:`RegistrationResult.decompose`).

Every warping function is strictly increasing by construction.  It is built
from an unconstrained latent function ``W`` exactly as a monotone smooth is
(see :func:`fabel.smoothing.smooth` with ``constraint="morph"``):

.. math::

    h(t) = a + (b - a)\,\frac{\int_a^t e^{W(u)}\,du}{\int_a^b e^{W(u)}\,du},
    \qquad W(u) = \sum_k c_k \varphi_k(u),

so ``h(a) = a``, ``h(b) = b`` and ``h' > 0``.  The registered curve is
``x(h(t))`` (``x(h(t) + δ)`` with a per-curve shift ``δ`` for periodic data).

Continuous registration chooses ``c`` for each curve by minimising, on the
uniform grid ``t_1 .. t_n`` with ``n = max(201, 10 K + 1)`` points (``K`` the
size of the curves' basis), either the least-squares criterion

.. math:: F_1 = \frac1n \sum_j \big(x_0(t_j) - x(h(t_j))\big)^2

or the minimum-eigenvalue criterion of Ramsay & Silverman (2005, §7.6),

.. math::

    F_2 = 2\,\mu_{\min}\!\begin{pmatrix}
        \overline{x_0^2} & \overline{x_0\,x\circ h} \\
        \overline{x_0\,x\circ h} & \overline{(x\circ h)^2}
    \end{pmatrix},

where the bar is the grid mean, plus the roughness penalty ``λ cᵀRc``.  These are
the discretisations R ``fda`` 6.3.0 uses, so ``lam`` means the same thing in
both.

Multivariate curves (``n_vars > 1``, such as the hip and knee angles of the
gait data, or the ``x`` and ``y`` coordinates of handwriting) get one warp per
curve, shared by all variables: the fitting criterion is the weighted sum
``Σ_v w_v F(x_{0v}, x_v∘h)`` of the criteria of the variables (Ramsay &
Silverman 2005, §7.6), with equal weights by default.  R ``fda`` 6.3.0's
``register.fd`` accepts multivariate curves but fits their warps to the first
variable alone (measured, black box); ``var_weights=[1, 0, ...]`` reproduces
that.  The first coefficient is pinned to zero, which removes the one exactly
flat direction (``W -> W + s`` leaves ``h`` unchanged).  The minimisation is a
safeguarded Newton iteration with the exact Hessian: the first and second
derivatives of ``h`` with respect to ``c`` are cumulative integrals of
``φ_k e^W`` and ``φ_k φ_l e^W``, computed by Gauss-Legendre quadrature.

Examples
--------
>>> import numpy as np
>>> from fabel import BSpline, FData
>>> from fabel.registration import register
>>> basis = BSpline(domain=(0.0, 1.0), n_basis=15)
>>> t = np.linspace(0.0, 1.0, 400)
>>> shifts = [0.42, 0.5, 0.58]
>>> curves = np.stack([np.exp(-((t - s) / 0.1) ** 2) for s in shifts], axis=1)
>>> fd = FData(np.linalg.lstsq(basis(t), curves, rcond=None)[0], basis)
>>> wbasis = BSpline(domain=(0.0, 1.0), n_basis=5)
>>> res = register(fd, criterion="least_squares", warp_basis=wbasis, lam=1e-3)
>>> peaks = t[np.argmax(res.registered(t), axis=0)]
>>> bool(np.ptp(peaks) < 0.5 * np.ptp(shifts))
True
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass
from math import atan2, cos, isfinite, sin, sqrt
from typing import Any, NamedTuple

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted, validate_data

from fabel import _linalg
from fabel._backend import asarray, default_namespace, is_torch, to_numpy
from fabel._operator import LDO
from fabel.basis import Basis, BSpline
from fabel.core import FData
from fabel.smoothing import smooth

__all__ = [
    "AmpPhaseDecomposition",
    "RegistrationResult",
    "Registrator",
    "landmark_register",
    "register",
]

Array = Any

#: The criterion, gradient and Hessian of one curve as a function of its
#: parameters, and a factory building it from a :class:`_CurveProblem`.
Objective = Callable[[Array], tuple[float, Array, Array]]
ObjectiveFactory = Callable[["_CurveProblem"], Objective]

#: The registration grid has ``max(_MIN_FINE_POINTS, _POINTS_PER_BASIS * K + 1)``
#: points, ``K`` being the size of the curves' basis -- R ``fda``'s choice.
_MIN_FINE_POINTS = 201
_POINTS_PER_BASIS = 10

#: Gauss-Legendre nodes per grid interval for the cumulative integrals of ``e^W``.
#: ``W`` is a polynomial on every interval, so eight nodes integrate ``e^W`` to
#: rounding error for any warp that can occur in practice.
_QUAD_DEGREE = 8

#: Gauss-Legendre nodes per panel for the warps returned to the user (the same
#: rule :func:`fabel.smoothing.smooth` uses for monotone fits).
_WARP_QUAD_DEGREE = 12

_ARMIJO = 1e-4
_EIGEN_FLOOR = 1e-12
_VALUE_NOISE = 1e-12
_STEP_FLOOR = 1e-15
_GRADIENT_NOISE = 1e-6

#: Default roughness weight for landmark warps: the landmark problem has more
#: warp coefficients than landmarks, so some penalty is needed to make it
#: well-posed; this one is small enough that the landmarks are met closely.
_DEFAULT_LANDMARK_LAMBDA = 1e-4

#: Newton polishing steps for the inverse warp, started from linear interpolation.
_INVERSE_NEWTON_STEPS = 6
_INVERSE_OVERSAMPLE = 8

_CRITERIA = ("eigen", "least_squares")


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


class AmpPhaseDecomposition(NamedTuple):
    """Amplitude/phase split of the variation in a set of curves.

    Replaces the list returned by R's ``AmpPhaseDecomp``.  Being a named tuple it
    unpacks as ``amp_mse, phase_mse, rsq, c = res.decompose()``.

    Attributes
    ----------
    amp_mse : float
        Mean squared error due to amplitude variation (R ``MS.amp``).
    phase_mse : float
        Mean squared error due to phase variation (R ``MS.pha``).  It can be
        negative when the registration does not improve the alignment.
    rsq : float
        Share of the total variation that is due to phase,
        ``phase_mse / (amp_mse + phase_mse)`` (R ``RSQR``).
    c : float
        The constant of the decomposition (R ``C``); it equals one when the warp
        derivatives are uncorrelated with the squared registered curves.

    Examples
    --------
    >>> from fabel.registration import AmpPhaseDecomposition
    >>> AmpPhaseDecomposition(1.0, 3.0, 0.75, 1.0).rsq
    0.75
    """

    amp_mse: float
    phase_mse: float
    rsq: float
    c: float


@dataclass(frozen=True, eq=False)
class RegistrationResult:
    """Outcome of :func:`register`.

    Attributes
    ----------
    registered : FData
        The registered curves ``x_i(h_i(t))`` (R ``regfd``), in the basis of the
        input curves: the least-squares fit to the warped values on the
        registration grid.
    warp : FData
        The warping functions ``h_i`` (R ``warpfd``).  For continuous
        registration they are expressed in the warp basis and include the shift
        of a periodic registration; for landmark registration they are expressed
        in the basis of the curves.  Both are least-squares fits on the
        registration grid, as in R.  :meth:`warp_values` evaluates the exact
        warps.
    unregistered : FData
        The input curves (R ``yfd``).
    latent : FData
        The latent functions ``W_i`` that define the warps (R ``Wfd``).
    shift : array
        Per-curve time shift ``δ_i`` of a periodic registration; zeros otherwise.
    target : FData or None
        The target curve(s) of a continuous registration; ``None`` for landmarks.
    warp_inverse : FData or None
        The inverse warps ``h_i⁻¹`` of a landmark registration (R ``warpinvfd``),
        in the basis of the curves; ``None`` for continuous registration.
    criterion : array or None
        Final value of the fitting criterion per curve (continuous only).
    n_iter : array or None
        Newton iterations used per curve (continuous only).

    Notes
    -----
    For multivariate curves ``registered``, ``unregistered`` and ``target``
    keep their variables (coefficients ``(n_basis, n_curves, n_vars)``), while
    ``warp``, ``latent``, ``shift`` and ``warp_inverse`` hold one warp per
    curve, shared by its variables.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import register
    >>> basis = BSpline(domain=(0.0, 1.0), n_basis=8)
    >>> fd = FData(np.random.default_rng(0).standard_normal((8, 3)), basis)
    >>> res = register(fd, landmarks=[[0.4], [0.5], [0.6]])
    >>> res.registered.n_curves, res.warp_inverse.n_curves
    (3, 3)
    """

    registered: FData
    warp: FData
    unregistered: FData
    latent: FData
    shift: Array
    target: FData | None = None
    warp_inverse: FData | None = None
    criterion: Array | None = None
    n_iter: Array | None = None

    def warp_values(self, t: Any) -> Array:
        """Evaluate the exact warping functions ``h_i(t) + δ_i``.

        Parameters
        ----------
        t : array
            Points in the domain of the curves, shape ``(n_points,)``.

        Returns
        -------
        array
            Warp values of shape ``(n_points, n_curves)``: a tensor (constant,
            on the device of the result) for a registration of tensors or for
            tensor ``t``, a NumPy array otherwise.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import BSpline, FData
        >>> from fabel.registration import register
        >>> basis = BSpline(domain=(0.0, 1.0), n_basis=8)
        >>> fd = FData(np.random.default_rng(1).standard_normal((8, 2)), basis)
        >>> res = register(fd, landmarks=[[0.45], [0.55]])
        >>> values = res.warp_values(np.linspace(0.0, 1.0, 5))
        >>> bool(np.all(np.diff(values, axis=0) > 0))
        True
        """
        xp = default_namespace()
        points = xp.reshape(asarray(to_numpy(t), xp), (-1,))
        coefs = asarray(to_numpy(self.latent.coefs), xp)
        values = _warp_values(self.latent.basis, coefs, points)
        values = values + xp.reshape(asarray(to_numpy(self.shift), xp), (1, -1))
        like = self.latent.coefs if is_torch(self.latent.coefs) else t
        if not is_torch(like):
            return values
        from fabel._internal.registration_torch import as_tensor_like

        return as_tensor_like(values, like)

    def apply(self, fd: FData) -> FData:
        """Warp new curves with the warps of this registration.

        Replaces R's ``register.newfd``: curve ``i`` of ``fd`` becomes
        ``x_i(h_i(t) + δ_i)`` with the exact warp ``h_i`` and shift ``δ_i``
        estimated here (warped times wrap around the domain when the shifts of
        a periodic registration are not all zero).  The result is the
        least-squares fit of those values, on the registration grid of ``fd``'s
        basis, in ``fd``'s basis -- the same projection :func:`register` uses,
        so ``res.apply(res.unregistered)`` reproduces ``res.registered``.  A
        typical use is carrying warps estimated on one set of curves over to
        derived or companion curves (derivatives, other variables).

        Parameters
        ----------
        fd : FData
            Curves to warp, univariate or multivariate, with one curve per warp
            and on the domain of the warps.

        Returns
        -------
        FData
            The warped curves in ``fd``'s basis, with ``fd``'s variables.  With
            :class:`torch.Tensor` coefficients the result holds a tensor of the
            same dtype and device, differentiable with respect to
            ``fd.coefs`` (the warps are constants).

        Raises
        ------
        ValueError
            If ``fd`` has a different number of curves than there are warps,
            lies on another domain, or has NaN or infinite coefficients.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import BSpline, FData
        >>> from fabel.registration import register
        >>> basis = BSpline(domain=(0.0, 1.0), n_basis=8)
        >>> fd = FData(np.random.default_rng(5).standard_normal((8, 3)), basis)
        >>> res = register(fd, landmarks=[[0.4], [0.5], [0.6]])
        >>> bool(np.allclose(res.apply(fd).coefs, res.registered.coefs))
        True
        >>> res.apply(fd.derivative()).n_curves
        3
        """
        n_curves = self.latent.n_curves
        if fd.n_curves != n_curves:
            raise ValueError(f"fd must hold {n_curves} curves, one per warp, got {fd.n_curves}")
        if fd.domain != self.latent.domain:
            raise ValueError(
                f"fd domain {fd.domain} differs from the warp domain {self.latent.domain}"
            )
        data = _as_numpy_fdata(fd)
        _check_finite(data.coefs, "fd")
        xp = default_namespace()
        grid = _fine_grid(fd.domain, fd.basis.n_basis)
        coefs = asarray(to_numpy(self.latent.coefs), xp)
        shift = xp.reshape(asarray(to_numpy(self.shift), xp), (-1,))
        warps = _warp_values(self.latent.basis, coefs, grid) + shift[None, :]
        if bool(xp.any(shift != 0.0)):
            lower, upper = fd.domain
            warps = lower + xp.remainder(warps - lower, upper - lower)
        if not is_torch(fd.coefs):
            return _warp_curves(data, grid, warps)
        from fabel._internal.registration_torch import warp_curves_torch

        return warp_curves_torch(fd, grid, warps, fd.coefs)

    def decompose(self, domain: tuple[float, float] | None = None) -> AmpPhaseDecomposition:
        r"""Split the variation of the curves into amplitude and phase parts.

        Replaces R's ``AmpPhaseDecomp`` (Kneip & Ramsay 2008).  With ``x_i`` the
        unregistered curves, ``y_i`` the registered ones and ``h_i`` the warps,

        .. math::

            C = 1 + \frac{\frac{1}{N-1}\sum_i \int (Dh_i - \overline{Dh})
                          (y_i^2 - \overline{y^2})}{\frac1N \sum_i \int y_i^2},
            \quad
            \text{amp} = \frac{C}{N}\sum_i \int (y_i - \bar y)^2,
            \quad
            \text{phase} = C\int \bar y^2 - \int \bar x^2 ,

        and ``rsq = phase / (amp + phase)``.  The integrals use the trapezoidal
        rule on ``max(201, 10 K + 1)`` equally spaced points, ``K`` being the size
        of the basis of the unregistered curves, which is the rule R uses.  For
        multivariate curves the squares are squared Euclidean norms over the
        variables (every integral is summed over the variables); R's
        ``AmpPhaseDecomp`` rejects multivariate curves.

        Parameters
        ----------
        domain : tuple of float, optional
            Sub-interval over which to decompose.  Defaults to the whole domain.

        Returns
        -------
        AmpPhaseDecomposition
            ``(amp_mse, phase_mse, rsq, c)``.

        Raises
        ------
        ValueError
            If fewer than two curves are registered, or ``domain`` is not inside
            the domain of the curves.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import BSpline, FData
        >>> from fabel.registration import register
        >>> basis = BSpline(domain=(0.0, 1.0), n_basis=10)
        >>> fd = FData(np.random.default_rng(2).standard_normal((10, 4)), basis)
        >>> amp, phase, rsq, c = register(fd, landmarks=[[0.4], [0.45], [0.55], [0.6]]).decompose()
        >>> bool(abs(rsq - phase / (amp + phase)) < 1e-12)
        True
        """
        return _amp_phase(self.unregistered, self.registered, self.warp, domain)


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #


def _fine_grid(domain: tuple[float, float], n_basis: int) -> Array:
    """Return R ``fda``'s registration grid for curves with ``n_basis`` functions."""
    xp = default_namespace()
    count = max(_MIN_FINE_POINTS, _POINTS_PER_BASIS * n_basis + 1)
    return xp.linspace(domain[0], domain[1], count, dtype=xp.float64)


def _as_numpy_fdata(fd: FData) -> FData:
    """Return ``fd`` with NumPy coefficients, the form the optimiser works on."""
    xp = default_namespace()
    return FData(asarray(to_numpy(fd.coefs), xp), fd.basis)


def _variable_weights(weights: Any, n_vars: int) -> tuple[float, ...]:
    """Validate the per-variable criterion weights (all ones by default)."""
    if weights is None:
        return (1.0,) * n_vars
    xp = default_namespace()
    values = xp.reshape(asarray(to_numpy(weights), xp), (-1,))
    if values.shape[0] != n_vars:
        raise ValueError(f"var_weights must have {n_vars} entries, got {values.shape[0]}")
    _check_finite(values, "var_weights")
    if bool(xp.any(values < 0.0)) or not bool(xp.any(values > 0.0)):
        raise ValueError("var_weights must be non-negative with at least one positive entry")
    return tuple(float(value) for value in values)


def _by_variable(values: Array) -> Array:
    """Reshape one curve's values ``(n, 1)`` or ``(n, 1, n_vars)`` to ``(n, n_vars)``."""
    xp = default_namespace()
    return xp.reshape(values, (values.shape[0], -1))


def _relative_cumulative(basis: Basis, coefs: Array, t: Array) -> tuple[Array, Array, Array]:
    """Cumulative integrals of ``e^{W - m}`` at ``t`` and over the whole domain.

    ``m`` is, per column, the largest value of ``W`` at the quadrature nodes.
    The warp depends on ``W`` only up to a constant, and with this offset the
    integrand is at most one (``e^W`` cannot overflow however large unpenalised
    coefficients grow) and equals one at some node (the total cannot underflow).
    The panels are the breaks of ``basis`` refined by ``t``, so every panel
    sees one polynomial piece of ``W``.

    Returns
    -------
    tuple of array
        ``(values, total, offset)`` of shapes ``(n_t, n_curves)``,
        ``(1, n_curves)`` and ``(1, n_curves)``.
    """
    xp = default_namespace()
    lower, upper = basis.domain
    points = xp.reshape(t, (-1,))
    edges = xp.unique_values(
        xp.concat([asarray(list(basis._natural_breaks()), xp), points, asarray([lower, upper], xp)])
    )
    edges = edges[(edges >= lower) & (edges <= upper)]
    nodes, weights = _linalg.composite_gauss_legendre(to_numpy(edges), _WARP_QUAD_DEGREE)
    n_panels = int(edges.shape[0]) - 1
    latent = xp.matmul(basis(asarray(nodes, xp)), coefs)
    offset = xp.max(latent, axis=0, keepdims=True)
    weighted = xp.exp(latent - offset) * asarray(weights, xp)[:, None]
    panels = xp.sum(xp.reshape(weighted, (n_panels, _WARP_QUAD_DEGREE, -1)), axis=1)
    zero = xp.zeros((1, panels.shape[1]), dtype=panels.dtype)
    cumulative = xp.concat([zero, xp.cumulative_sum(panels, axis=0)])
    values = xp.take(cumulative, xp.searchsorted(edges, points), axis=0)
    return values, cumulative[-1:, :], offset


def _check_finite(values: Array, name: str) -> None:
    """Reject NaN or infinite entries, which no registration criterion survives."""
    xp = default_namespace()
    if not bool(xp.all(xp.isfinite(values))):
        raise ValueError(f"{name} must be finite (no NaN or infinite entries)")


def _warp_values(basis: Basis, coefs: Array, t: Array) -> Array:
    """Exact ``h(t) = a + (b - a) H(t) / H(b)`` for every column of ``coefs``."""
    xp = default_namespace()
    lower, upper = basis.domain
    values, total, _ = _relative_cumulative(basis, coefs, t)
    return xp.clip(lower + (upper - lower) * values / total, lower, upper)


def _warp_slopes(basis: Basis, coefs: Array, t: Array) -> Array:
    """Exact ``h'(t) = (b - a) e^{W(t)} / H(b)`` for every column of ``coefs``."""
    xp = default_namespace()
    lower, upper = basis.domain
    _, total, offset = _relative_cumulative(basis, coefs, asarray([upper], xp))
    return (upper - lower) * xp.exp(xp.matmul(basis(t), coefs) - offset) / total


def _project(basis: Basis, grid: Array, values: Array) -> FData:
    """Least-squares fit of ``values`` sampled on ``grid`` in ``basis``.

    ``values`` is ``(n_grid, n_curves)`` or ``(n_grid, n_curves, n_vars)``.
    """
    xp = default_namespace()
    if values.ndim == 3:
        n_grid, n_curves, n_vars = values.shape
        flat = _linalg.lstsq(basis(grid), xp.reshape(values, (n_grid, n_curves * n_vars)))
        return FData(xp.reshape(flat, (basis.n_basis, n_curves, n_vars)), basis)
    return FData(_linalg.lstsq(basis(grid), values), basis)


def _warp_curves(fd: FData, grid: Array, where: Array) -> FData:
    """Fit ``x_i(where[:, i])`` on ``grid`` in ``fd``'s basis, keeping the variables."""
    xp = default_namespace()
    warped = xp.stack([_by_variable(fd[i](where[:, i])) for i in range(fd.n_curves)], axis=1)
    if len(fd.coefs.shape) == 2:
        warped = warped[:, :, 0]
    return _project(fd.basis, grid, warped)


def _penalty_matrix(basis: Basis, penalty: int | LDO, lam: float) -> Array:
    """Return ``R`` for a positive ``lam``, a zero matrix otherwise."""
    xp = default_namespace()
    if lam < 0.0:
        raise ValueError(f"lam must be non-negative, got {lam}")
    if lam == 0.0:
        return xp.zeros((basis.n_basis, basis.n_basis), dtype=xp.float64)
    return asarray(to_numpy(basis.penalty(penalty)), xp)


# --------------------------------------------------------------------------- #
# continuous registration
# --------------------------------------------------------------------------- #


class _WarpQuadrature:
    """Cumulative integrals of ``e^W``, ``φ_k e^W`` and ``φ_k φ_l e^W`` on a grid.

    The quadrature panels are the grid intervals refined by the warp basis
    breaks, so every panel sees a single polynomial piece of ``W``.
    """

    def __init__(self, basis: Basis, grid: Array) -> None:
        xp = default_namespace()
        lower, upper = basis.domain
        edges = xp.unique_values(xp.concat([grid, asarray(list(basis._natural_breaks()), xp)]))
        edges = edges[(edges >= lower) & (edges <= upper)]
        nodes, weights = _linalg.composite_gauss_legendre(to_numpy(edges), _QUAD_DEGREE)
        self.n_panels = int(edges.shape[0]) - 1
        self.phi = basis(asarray(nodes, xp))
        self.weights = asarray(weights, xp)
        self.index = xp.searchsorted(edges, grid)
        self.span = (float(lower), float(upper))

    def moments(self, coefs: Array) -> tuple[Array, Array, Array, Array, Array, Array]:
        """Return ``H, H_k, H_kl`` at the grid points and their totals at ``b``.

        The integrals are of ``e^{W - m}`` with ``m`` the largest value of
        ``W`` at the quadrature nodes, so the integrand is at most one (no
        overflow) and equals one at some node (the totals cannot underflow).
        Every ratio built from the moments is the same as for ``e^W``.
        """
        xp = default_namespace()
        n_basis = self.phi.shape[1]
        latent = xp.matmul(self.phi, coefs)
        weighted = xp.exp(latent - xp.max(latent)) * self.weights
        order0 = xp.sum(xp.reshape(weighted, (self.n_panels, _QUAD_DEGREE)), axis=1)
        first = weighted[:, None] * self.phi
        order1 = xp.sum(xp.reshape(first, (self.n_panels, _QUAD_DEGREE, n_basis)), axis=1)
        second = first[:, :, None] * self.phi[:, None, :]
        order2 = xp.sum(xp.reshape(second, (self.n_panels, _QUAD_DEGREE, n_basis, n_basis)), axis=1)
        cum0 = xp.concat([xp.zeros(1, dtype=xp.float64), xp.cumulative_sum(order0)])
        cum1 = xp.concat(
            [xp.zeros((1, n_basis), dtype=xp.float64), xp.cumulative_sum(order1, axis=0)]
        )
        cum2 = xp.concat(
            [
                xp.zeros((1, n_basis, n_basis), dtype=xp.float64),
                xp.cumulative_sum(order2, axis=0),
            ]
        )
        idx = self.index
        return (
            xp.take(cum0, idx, axis=0),
            xp.take(cum1, idx, axis=0),
            xp.take(cum2, idx, axis=0),
            cum0[-1],
            cum1[-1, :],
            cum2[-1, :, :],
        )

    def warp(self, coefs: Array) -> tuple[Array, Array, Array]:
        """Return ``h`` and its first and second derivatives in ``coefs``."""
        xp = default_namespace()
        lower, upper = self.span
        width = upper - lower
        a0, a1, a2, t0, t1, t2 = self.moments(coefs)
        h = xp.clip(lower + width * a0 / t0, lower, upper)
        grad = width * (a1 / t0 - a0[:, None] * t1[None, :] / t0**2)
        cross = a1[:, :, None] * t1[None, None, :]
        hess = width * (
            a2 / t0
            - (cross + xp.permute_dims(cross, (0, 2, 1))) / t0**2
            - a0[:, None, None] * t2[None, :, :] / t0**2
            + 2.0 * a0[:, None, None] * (t1[:, None] * t1[None, :])[None, :, :] / t0**3
        )
        return h, grad, hess


@dataclass(frozen=True, eq=False)
class _CurveProblem:
    """The registration criterion of one curve as a function of its parameters."""

    curve: FData
    target: Array
    quadrature: _WarpQuadrature
    penalty: Array
    lam: float
    criterion: str
    periodic: bool
    has_curvature: bool
    weights: tuple[float, ...] = (1.0,)

    @property
    def n_coefs(self) -> int:
        """Number of warp coefficients, the pinned first one included."""
        return int(self.penalty.shape[0])

    def split(self, params: Array) -> tuple[Array, float]:
        """Return the full coefficient vector and the shift for ``params``."""
        xp = default_namespace()
        free = self.n_coefs - 1
        coefs = xp.concat([xp.zeros(1, dtype=xp.float64), params[:free]])
        shift = float(params[free]) if self.periodic else 0.0
        return coefs, shift

    def _position(self, values: Array) -> Array:
        """Map warped times back into the domain (periodically if required)."""
        xp = default_namespace()
        lower, upper = self.quadrature.span
        if self.periodic:
            return lower + xp.remainder(values - lower, upper - lower)
        return xp.clip(values, lower, upper)

    def evaluate(self, params: Array) -> tuple[float, Array, Array]:
        """Return the criterion, its gradient and its Hessian at ``params``.

        For multivariate curves the criterion is the ``weights``-weighted sum
        of the criteria of the variables, all warped by the same ``h``.
        """
        xp = default_namespace()
        coefs, shift = self.split(params)
        h, h_c, h_cc = self.quadrature.warp(coefs)
        n_points = h.shape[0]
        jac_h = h_c[:, 1:]
        hess_h = h_cc[:, 1:, 1:]
        if self.periodic:
            jac_h = xp.concat([jac_h, xp.ones((n_points, 1), dtype=xp.float64)], axis=1)
            width = jac_h.shape[1]
            padded = xp.zeros((n_points, width, width), dtype=xp.float64)
            padded[:, : width - 1, : width - 1] = hess_h
            hess_h = padded
        where = self._position(h + shift)
        values = _by_variable(self.curve(where))
        slopes = _by_variable(self.curve(where, 1))
        curvatures = _by_variable(self.curve(where, 2)) if self.has_curvature else None
        target = xp.reshape(self.target, (n_points, -1))
        part = _least_squares if self.criterion == "least_squares" else _min_eigenvalue
        size = jac_h.shape[1]
        fit = 0.0
        grad = xp.zeros(size, dtype=xp.float64)
        second = xp.zeros((size, size), dtype=xp.float64)
        for var, weight in enumerate(self.weights):
            if weight == 0.0:
                continue
            slope = slopes[:, var]
            jac = slope[:, None] * jac_h
            hess = slope[:, None, None] * hess_h
            if curvatures is not None:
                curvature = curvatures[:, var]
                hess = hess + curvature[:, None, None] * jac_h[:, :, None] * jac_h[:, None, :]
            var_fit, var_grad, var_second = part(target[:, var], values[:, var], jac, hess)
            fit = fit + weight * var_fit
            grad = grad + weight * var_grad
            second = second + weight * var_second
        free = self.n_coefs - 1
        pen = self.penalty[1:, 1:]
        rough = xp.matmul(pen, params[:free])
        total = fit + self.lam * float(xp.sum(params[:free] * rough))
        grad = xp.concat([grad[:free] + 2.0 * self.lam * rough, grad[free:]])
        second = second.copy()
        second[:free, :free] = second[:free, :free] + 2.0 * self.lam * pen
        return total, grad, second


def _least_squares(
    target: Array, value: Array, jac: Array, hess: Array
) -> tuple[float, Array, Array]:
    """Mean squared distance to the target and its derivatives."""
    xp = default_namespace()
    n_points = value.shape[0]
    resid = target - value
    fit = float(xp.sum(resid * resid)) / n_points
    grad = -2.0 * xp.matmul(xp.matrix_transpose(jac), resid) / n_points
    second = (
        2.0
        * (xp.matmul(xp.matrix_transpose(jac), jac) - xp.sum(resid[:, None, None] * hess, axis=0))
        / n_points
    )
    return fit, grad, second


def _min_eigenvalue(
    target: Array, value: Array, jac: Array, hess: Array
) -> tuple[float, Array, Array]:
    """Twice the smaller eigenvalue of the grid-mean cross-product matrix."""
    xp = default_namespace()
    n_points = value.shape[0]
    aa = float(xp.sum(target * target)) / n_points
    bb = float(xp.sum(target * value)) / n_points
    dd = float(xp.sum(value * value)) / n_points
    half_gap = sqrt(0.25 * (aa - dd) ** 2 + bb * bb)
    # Unit eigenvectors of the larger (u) and the smaller (v) eigenvalue from
    # the rotation angle.  The small eigenvalue is then the mean square of the
    # residual v1 x0 + v2 (x o h): a sum of squares, free of the cancellation
    # in (a + d)/2 - half_gap that would put a 1e-11 noise floor under Newton.
    angle = 0.5 * atan2(2.0 * bb, aa - dd)
    u1, u2 = cos(angle), sin(angle)
    v1, v2 = -u2, u1
    resid = v1 * target + v2 * value
    other = u1 * target + u2 * value
    low = float(xp.sum(resid * resid)) / n_points
    jac_t = xp.matrix_transpose(jac)
    grad = 2.0 * v2 * xp.matmul(jac_t, resid) / n_points
    second = (
        2.0
        * (v2 * v2 * xp.matmul(jac_t, jac) + v2 * xp.sum(resid[:, None, None] * hess, axis=0))
        / n_points
    )
    if half_gap > 0.0:
        mixed = (v2 * xp.matmul(jac_t, other) + u2 * xp.matmul(jac_t, resid)) / n_points
        second = second - (mixed[:, None] * mixed[None, :]) / half_gap
    return 2.0 * low, 2.0 * grad, 2.0 * second


def _newton_step(hess: Array, grad: Array) -> Array:
    """Newton step with the Hessian's non-positive curvature reflected upward.

    Falls back to the steepest-descent step ``-grad`` when the Hessian is not
    finite, is zero, or its eigendecomposition fails.
    """
    xp = default_namespace()
    sym = 0.5 * (hess + xp.matrix_transpose(hess))
    if not bool(xp.all(xp.isfinite(sym))):
        # A Hessian that overflowed carries no usable curvature.
        return -grad
    try:
        values, vectors = xp.linalg.eigh(sym)
    except ValueError:
        # LAPACK's eigensolver can fail to converge on badly scaled matrices
        # (numpy's LinAlgError is a ValueError): take a steepest-descent step.
        return -grad
    scale = float(xp.max(xp.abs(values))) if values.shape[0] else 0.0
    if scale == 0.0:
        # No curvature information at all: fall back to steepest descent.
        return -grad
    floor = _EIGEN_FLOOR * scale
    safe = xp.where(xp.abs(values) > floor, xp.abs(values), floor)
    return -xp.matmul(vectors, xp.matmul(xp.matrix_transpose(vectors), grad) / safe)


def _minimise(
    evaluate: Objective,
    start: Array,
    max_iter: int,
    tol: float,
) -> tuple[Array, float, int, bool]:
    """Safeguarded Newton descent with an Armijo backtracking line search.

    A trial point whose criterion, gradient or Hessian is not finite is
    rejected like one that does not decrease the criterion, so the line search
    shortens the step instead of accepting an overflowed point.

    Raises
    ------
    ValueError
        If the criterion, gradient or Hessian is not finite at ``start``.
    """
    xp = default_namespace()

    def finite(point_value: float, point_grad: Array, point_hess: Array) -> bool:
        return (
            isfinite(point_value)
            and bool(xp.all(xp.isfinite(point_grad)))
            and bool(xp.all(xp.isfinite(point_hess)))
        )

    def stationary(point_value: float, point_grad: Array) -> bool:
        if not point_grad.shape[0]:
            return True
        return float(xp.max(xp.abs(point_grad))) <= tol * max(1.0, abs(point_value))

    params = start
    value, grad, hess = evaluate(params)
    if not finite(value, grad, hess):
        raise ValueError("the registration criterion is not finite at the starting point")
    for iteration in range(max_iter + 1):
        if stationary(value, grad):
            return params, value, iteration, True
        if iteration == max_iter:
            break
        step = _newton_step(hess, grad)
        slope = float(xp.sum(grad * step))
        slack = _VALUE_NOISE * max(1.0, abs(value))
        reach = float(xp.max(xp.abs(step)))
        floor = _STEP_FLOOR * (1.0 + float(xp.max(xp.abs(params))))
        scale = 1.0
        accepted = False
        while scale * reach > floor:
            candidate = params + scale * step
            new_value, new_grad, new_hess = evaluate(candidate)
            if not finite(new_value, new_grad, new_hess):
                scale *= 0.5
                continue
            if new_value <= value + _ARMIJO * scale * slope:
                accepted = True
                break
            # Next to the minimum the decrease is below the rounding noise of the
            # criterion; a step that lands on a stationary point is still taken.
            if stationary(new_value, new_grad) and new_value <= value + slack:
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            # No step length decreases the criterion: the iterate sits on the
            # rounding floor of the criterion.  That is convergence only if the
            # gradient is down at the noise level too.
            noise = float(xp.max(xp.abs(grad))) <= _GRADIENT_NOISE * max(1.0, abs(value))
            return params, value, iteration, noise
        params, value, grad, hess = candidate, new_value, new_grad, new_hess
    return params, value, max_iter, False


def _initial_latent(init: Any, warp_basis: Basis, n_curves: int) -> tuple[Array, Basis]:
    """Resolve ``init`` into a ``(K, n_curves)`` coefficient matrix."""
    xp = default_namespace()
    if isinstance(init, FData):
        basis = init.basis
        coefs = asarray(to_numpy(init.coefs), xp)
    else:
        basis = warp_basis
        coefs = (
            xp.zeros((basis.n_basis, n_curves), dtype=xp.float64)
            if init is None
            else asarray(to_numpy(init), xp)
        )
    if coefs.ndim == 1:
        coefs = coefs[:, None]
    coefs = xp.reshape(coefs, (coefs.shape[0], -1))
    if coefs.shape[0] != basis.n_basis:
        raise ValueError(f"init must have {basis.n_basis} rows, got {coefs.shape[0]}")
    if coefs.shape[1] == 1:
        coefs = xp.tile(coefs, (1, n_curves))
    if coefs.shape[1] != n_curves:
        raise ValueError(f"init must have 1 or {n_curves} columns, got {coefs.shape[1]}")
    _check_finite(coefs, "init")
    # Adding a constant to W leaves h unchanged; pin the first coefficient to 0.
    return coefs - coefs[:1, :], basis


def _continuous(
    fd: FData,
    target: FData | None,
    warp_basis: Basis | None,
    lam: float,
    penalty: int | LDO,
    criterion: str,
    periodic: bool,
    init: Any,
    init_shift: Any,
    max_iter: int,
    tol: float,
    objective: ObjectiveFactory | None = None,
    weights: tuple[float, ...] = (1.0,),
) -> RegistrationResult:
    """Continuous registration of every curve of ``fd`` to ``target``.

    ``objective`` builds the function handed to the optimiser from each curve's
    problem; the default is the analytic :meth:`_CurveProblem.evaluate`.  The
    PyTorch path passes one that differentiates the criterion by autograd.
    """
    xp = default_namespace()
    if criterion not in _CRITERIA:
        raise ValueError(f"criterion must be one of {_CRITERIA}, got {criterion!r}")
    if max_iter < 0:
        raise ValueError(f"max_iter must be non-negative, got {max_iter}")
    n_curves = fd.n_curves
    goal = fd.mean() if target is None else _as_numpy_fdata(target)
    if goal.n_vars != fd.n_vars:
        raise ValueError(f"target must have {fd.n_vars} variables per curve, got {goal.n_vars}")
    _check_finite(goal.coefs, "target")
    if goal.n_curves not in (1, n_curves):
        raise ValueError(f"target must hold 1 or {n_curves} curves, got {goal.n_curves}")
    default_basis = BSpline(domain=fd.domain, n_basis=2, order=2)
    coefs, basis = _initial_latent(
        init, default_basis if warp_basis is None else warp_basis, n_curves
    )
    if not isinstance(basis, BSpline):
        raise ValueError("the warp basis must be a B-spline basis")
    if basis.domain != fd.domain:
        raise ValueError(f"warp basis domain {basis.domain} differs from {fd.domain}")
    shifts = (
        xp.zeros(n_curves, dtype=xp.float64)
        if init_shift is None
        else xp.reshape(asarray(to_numpy(init_shift), xp), (-1,)) * xp.ones(n_curves)
    )
    _check_finite(shifts, "init_shift")
    grid = _fine_grid(fd.domain, fd.basis.n_basis)
    goal_values = xp.reshape(goal(grid), (grid.shape[0], goal.n_curves, fd.n_vars))
    quadrature = _WarpQuadrature(basis, grid)
    pen = _penalty_matrix(basis, penalty, lam)
    curvature = not (isinstance(fd.basis, BSpline) and fd.basis.order < 3)
    coef_columns: list[Array] = []
    shift_values: list[float] = []
    criteria: list[float] = []
    iterations: list[int] = []
    unconverged: list[int] = []
    for i in range(n_curves):
        problem = _CurveProblem(
            curve=fd[i],
            target=goal_values[:, 0 if goal.n_curves == 1 else i, :],
            quadrature=quadrature,
            penalty=pen,
            lam=lam,
            criterion=criterion,
            periodic=periodic,
            has_curvature=curvature,
            weights=weights,
        )
        start = coefs[1:, i]
        if periodic:
            start = xp.concat([start, shifts[i : i + 1]])
        evaluate = problem.evaluate if objective is None else objective(problem)
        params, value, used, converged = _minimise(evaluate, start, max_iter, tol)
        if not converged:
            unconverged.append(i)
        full, shift = problem.split(params)
        coef_columns.append(full)
        shift_values.append(shift)
        criteria.append(value)
        iterations.append(used)
    if unconverged and max_iter > 0:
        warnings.warn(
            f"registration did not converge in {max_iter} iterations for curves {unconverged}",
            RuntimeWarning,
            stacklevel=3,
        )
    latent = xp.stack(coef_columns, axis=1)
    shift_array = asarray(shift_values, xp)
    lower, upper = fd.domain
    warps = _warp_values(basis, latent, grid) + shift_array[None, :]
    where = lower + xp.remainder(warps - lower, upper - lower) if periodic else warps
    return RegistrationResult(
        registered=_warp_curves(fd, grid, where),
        warp=_project(basis, grid, warps),
        unregistered=fd,
        latent=FData(latent, basis),
        shift=shift_array,
        target=goal,
        criterion=asarray(criteria, xp),
        n_iter=xp.asarray(iterations, dtype=xp.int64),
    )


# --------------------------------------------------------------------------- #
# landmark registration
# --------------------------------------------------------------------------- #


def _landmark_matrix(landmarks: Any, n_curves: int, domain: tuple[float, float]) -> Array:
    """Validate landmarks into an ``(n_curves, n_landmarks)`` matrix."""
    xp = default_namespace()
    marks = asarray(to_numpy(landmarks), xp)
    if marks.ndim == 1:
        marks = marks[:, None]
    if marks.ndim != 2 or marks.shape[0] != n_curves:
        raise ValueError(
            f"landmarks must have one row per curve ({n_curves}), got shape {tuple(marks.shape)}"
        )
    _check_finite(marks, "landmarks")
    lower, upper = domain
    if bool(xp.any((marks <= lower) | (marks >= upper))):
        raise ValueError(f"landmarks must lie strictly inside the domain {domain}")
    if marks.shape[1] > 1 and bool(xp.any(marks[:, 1:] <= marks[:, :-1])):
        raise ValueError("landmarks must be strictly increasing along each row")
    return marks


def _inverse_warp(basis: Basis, coefs: Array, grid: Array) -> Array:
    """Invert the warps ``h_i`` at the grid points to rounding error."""
    xp = default_namespace()
    lower, upper = basis.domain
    dense = xp.linspace(lower, upper, _INVERSE_OVERSAMPLE * (grid.shape[0] - 1) + 1)
    forward = _warp_values(basis, coefs, dense)
    columns: list[Array] = []
    for i in range(coefs.shape[1]):
        column = coefs[:, i : i + 1]
        guess = asarray(
            _linalg_interp(to_numpy(grid), to_numpy(forward[:, i]), to_numpy(dense)), xp
        )
        for _ in range(_INVERSE_NEWTON_STEPS):
            value = _warp_values(basis, column, guess)[:, 0]
            slope = _warp_slopes(basis, column, guess)[:, 0]
            guess = xp.clip(guess - (value - grid) / slope, lower, upper)
        columns.append(guess)
    return xp.stack(columns, axis=1)


def _linalg_interp(x: Array, xs: Array, ys: Array) -> Array:
    """Piecewise-linear interpolation ``y(x)`` through the increasing nodes ``xs``."""
    xp = default_namespace()
    points = asarray(x, xp)
    nodes = asarray(xs, xp)
    values = asarray(ys, xp)
    right = xp.clip(xp.searchsorted(nodes, points), 1, nodes.shape[0] - 1)
    left = right - 1
    x0 = xp.take(nodes, left)
    x1 = xp.take(nodes, right)
    y0 = xp.take(values, left)
    y1 = xp.take(values, right)
    frac = xp.where(x1 > x0, (points - x0) / xp.where(x1 > x0, x1 - x0, 1.0), 0.0)
    return y0 + frac * (y1 - y0)


def _landmark(
    fd: FData,
    landmarks: Any,
    target_landmarks: Any,
    warp_basis: Basis | None,
    lam: float,
    penalty: int | LDO,
) -> RegistrationResult:
    """Landmark registration: warp each curve so its landmarks meet the targets."""
    xp = default_namespace()
    n_curves = fd.n_curves
    lower, upper = fd.domain
    marks = _landmark_matrix(landmarks, n_curves, fd.domain)
    goal = (
        xp.mean(marks, axis=0)
        if target_landmarks is None
        else xp.reshape(asarray(to_numpy(target_landmarks), xp), (-1,))
    )
    if goal.shape[0] != marks.shape[1]:
        raise ValueError(
            f"target_landmarks must have {marks.shape[1]} entries, got {goal.shape[0]}"
        )
    _landmark_matrix(goal[None, :], 1, fd.domain)
    basis = (
        BSpline(domain=fd.domain, order=4, breaks=[lower, *to_numpy(goal).tolist(), upper])
        if warp_basis is None
        else warp_basis
    )
    if not isinstance(basis, BSpline) or basis.domain != fd.domain:
        raise ValueError("the warp basis must be a B-spline basis on the domain of the curves")
    knots_x = xp.concat([asarray([lower], xp), goal, asarray([upper], xp)])
    columns: list[Array] = []
    for i in range(n_curves):
        knots_y = xp.concat([asarray([lower], xp), marks[i, :], asarray([upper], xp)])
        fit = smooth(knots_y, knots_x, basis=basis, lam=lam, penalty=penalty, constraint="morph")
        coefs = asarray(to_numpy(fit.fd.coefs), xp)[:, 0]
        columns.append(coefs - xp.mean(coefs))
    latent = xp.stack(columns, axis=1)
    grid = _fine_grid(fd.domain, fd.basis.n_basis)
    warps = _warp_values(basis, latent, grid)
    inverse = _inverse_warp(basis, latent, grid)
    return RegistrationResult(
        registered=_warp_curves(fd, grid, warps),
        warp=_project(fd.basis, grid, warps),
        unregistered=fd,
        latent=FData(latent, basis),
        shift=xp.zeros(n_curves, dtype=xp.float64),
        warp_inverse=_project(fd.basis, grid, inverse),
    )


# --------------------------------------------------------------------------- #
# amplitude / phase decomposition
# --------------------------------------------------------------------------- #


def _amp_phase(
    unregistered: FData,
    registered: FData,
    warp: FData,
    domain: tuple[float, float] | None,
) -> AmpPhaseDecomposition:
    """Kneip-Ramsay amplitude/phase decomposition on a trapezoidal grid."""
    xp = default_namespace()
    n_curves = unregistered.n_curves
    if n_curves < 2:
        raise ValueError("the decomposition needs at least two curves")
    span = unregistered.domain if domain is None else (float(domain[0]), float(domain[1]))
    lower, upper = unregistered.domain
    if not lower <= span[0] < span[1] <= upper:
        raise ValueError(f"domain {span} must be an interval inside {unregistered.domain}")
    count = max(_MIN_FINE_POINTS, _POINTS_PER_BASIS * unregistered.basis.n_basis + 1)
    grid = xp.linspace(span[0], span[1], count, dtype=xp.float64)
    step = (span[1] - span[0]) / (count - 1)
    weights = xp.full(count, step, dtype=xp.float64)
    weights[0] = weights[-1] = 0.5 * step
    # Curves as (n_grid, n_curves, n_vars); every square is a squared norm
    # over the variables, which is the plain square for univariate curves.
    shape = (count, n_curves, unregistered.n_vars)
    x = xp.reshape(asarray(to_numpy(unregistered(grid)), xp), shape)
    y = xp.reshape(asarray(to_numpy(registered(grid)), xp), shape)
    slope = asarray(to_numpy(warp(grid, 1)), xp)
    x_mean = xp.mean(x, axis=1)
    y_mean = xp.mean(y, axis=1)
    y_sq = xp.sum(y * y, axis=2)
    slope_dev = slope - xp.mean(slope, axis=1, keepdims=True)
    y_sq_dev = y_sq - xp.mean(y_sq, axis=1, keepdims=True)
    covariance = float(xp.sum(weights[:, None] * slope_dev * y_sq_dev)) / (n_curves - 1)
    power = float(xp.sum(weights[:, None] * y_sq)) / n_curves
    const = 1.0 + covariance / power
    spread = float(xp.sum(weights[:, None, None] * (y - y_mean[:, None, :]) ** 2)) / n_curves
    amp = const * spread
    y_norm = xp.sum(y_mean**2, axis=1)
    x_norm = xp.sum(x_mean**2, axis=1)
    phase = const * float(xp.sum(weights * y_norm)) - float(xp.sum(weights * x_norm))
    return AmpPhaseDecomposition(amp, phase, phase / (amp + phase), const)


# --------------------------------------------------------------------------- #
# entry points
# --------------------------------------------------------------------------- #


def register(
    fd: FData,
    target: FData | None = None,
    *,
    landmarks: Any = None,
    target_landmarks: Any = None,
    warp_basis: Basis | None = None,
    lam: float | None = None,
    penalty: int | LDO = 2,
    criterion: str = "eigen",
    periodic: bool = False,
    init: Any = None,
    init_shift: Any = None,
    max_iter: int = 100,
    tol: float = 1e-10,
    var_weights: Any = None,
) -> RegistrationResult:
    r"""Register curves by warping their time axis.

    Without ``landmarks`` this is continuous registration (R ``register.fd``):
    each curve is warped to minimise its distance to ``target`` (the mean curve
    by default).  With ``landmarks`` it is landmark registration (R
    ``landmarkreg``): each curve is warped so that its landmark times move to
    ``target_landmarks`` (the mean landmark times by default).

    Multivariate curves get one warp per curve, shared by all their
    variables; the continuous criterion is then the sum over the variables,
    weighted by ``var_weights``.

    Parameters
    ----------
    fd : FData
        The curves to register, univariate or multivariate.
    target : FData, optional
        One target curve, or one per curve, with the variables of ``fd``.
        Defaults to ``fd.mean()``.  Ignored for landmark registration.
    landmarks : array, optional
        Landmark times, shape ``(n_curves, n_landmarks)`` (or ``(n_curves,)``
        for one landmark), strictly inside the domain and increasing per curve.
    target_landmarks : array, optional
        Where the landmarks should land, shape ``(n_landmarks,)``.  Defaults to
        the mean landmark times.
    warp_basis : Basis, optional
        B-spline basis for the latent functions ``W``.  Defaults to the linear
        two-function basis for continuous registration (R's default) and to a
        cubic spline with knots at the target landmarks otherwise.
    lam : float, optional
        Roughness weight on ``W``.  Defaults to ``0`` for continuous registration
        (R's default) and to ``1e-4`` for landmark registration, where the
        latent coefficients outnumber the landmarks.
    penalty : int or LDO, optional
        Roughness operator on ``W``.  Default ``2``.
    criterion : {"eigen", "least_squares"}, optional
        Continuous criterion: the minimum-eigenvalue criterion (R ``crit=2``,
        the default) or least squares (R ``crit=1``).
    periodic : bool, optional
        Treat the curves as periodic over their domain and estimate a time shift
        per curve as well as the warp (R ``periodic=TRUE``).
    init : FData or array, optional
        Starting latent functions: an FData on the warp basis (whose basis then
        is the warp basis), or a ``(K,)`` / ``(K, n_curves)`` coefficient array.
        Defaults to zero, the identity warp.
    init_shift : float or array, optional
        Starting shifts of a periodic registration.  Default zero.
    max_iter : int, optional
        Maximum Newton iterations per curve.  ``0`` returns the warps defined by
        ``init`` without optimising.  Default ``100``.
    tol : float, optional
        Convergence tolerance on the largest gradient component, relative to
        ``max(1, criterion)``.  Default ``1e-10``.
    var_weights : array, optional
        Non-negative weight of each variable in the continuous criterion of
        multivariate curves, shape ``(n_vars,)``, at least one positive.
        Defaults to ones: the criteria of the variables are summed.  R
        ``fda`` 6.3.0's ``register.fd`` fits multivariate warps to the first
        variable only, which ``var_weights=[1, 0, ...]`` reproduces.  Ignored
        for landmark registration.

    Returns
    -------
    RegistrationResult
        The registered curves, warps, latent functions and shifts.  With
        PyTorch input (see Notes) every array and coefficient matrix in it is a
        tensor of the input's dtype and device.

    Raises
    ------
    ValueError
        On a non-B-spline warp basis, badly shaped landmarks or targets, a
        target whose number of variables differs from ``fd``'s, invalid
        ``var_weights``, an unknown criterion, a negative ``lam``, or NaN or
        infinite values in ``fd``, ``target``, ``init`` or ``init_shift``.

    Warns
    -----
    RuntimeWarning
        If some curve has not converged after ``max_iter`` iterations.  This
        includes a criterion without a finite minimum: with ``lam=0`` the
        eigenvalue criterion can keep decreasing as the warp squeezes a curve
        into a vanishing stretch of the domain, and ``W`` grows without bound.
        The line search rejects steps at which the criterion or its
        derivatives overflow and a failed eigendecomposition of the Hessian
        falls back to steepest descent, so such curves stop at the last finite
        iterate (an increasing warp) with this warning rather than an
        exception.  A positive ``lam`` keeps ``W`` bounded.

    Notes
    -----
    **PyTorch input.**  When ``fd`` (or, for continuous registration, an FData
    ``target``) has :class:`torch.Tensor` coefficients, the result holds
    tensors in the dtype and on the device of that input, and ``registered``
    is differentiable with respect to the input coefficients.  Continuous
    registration then evaluates the criterion in PyTorch and takes its
    gradient and Hessian by automatic differentiation instead of the analytic
    formulas; the safeguarded Newton iteration, the line search and the
    stopping rules are the same, so both paths reach the same optimum.  The
    optimisation always runs in float64 on the CPU (a registration needs
    ``tol``-level accuracy, and some devices have no float64).  Gradients flow
    through the final evaluation ``x_i(h_i(t))`` with the optimal warps held
    fixed: the dependence of the optimum on the curves (the implicit
    gradient) is not propagated.  Landmark registration needs no optimisation;
    its warps depend only on the landmarks, so the gradient of ``registered``
    with respect to the coefficients is exact.  ``warp``, ``latent``,
    ``warp_inverse``, ``shift``, ``criterion`` and ``n_iter`` are constant
    tensors.  Multivariate tensor curves work the same way.  ``import
    fabel`` never imports torch; this path imports it on first use.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import register
    >>> basis = BSpline(domain=(0.0, 1.0), n_basis=12)
    >>> t = np.linspace(0.0, 1.0, 300)
    >>> curves = np.stack([np.sin(np.pi * t**p) for p in (0.8, 1.0, 1.25)], axis=1)
    >>> fd = FData(np.linalg.lstsq(basis(t), curves, rcond=None)[0], basis)
    >>> peaks = t[np.argmax(curves, axis=0)]
    >>> res = register(fd, landmarks=peaks)
    >>> bool(np.allclose(res.warp_values(np.full(1, peaks.mean())), peaks, atol=1e-3))
    True

    Multivariate curves (here two variables that share their phase) get one
    warp per curve:

    >>> both = FData(np.stack([fd.coefs, 2.0 * fd.coefs], axis=2), basis)
    >>> multi = register(both, criterion="least_squares", warp_basis=BSpline(n_basis=4))
    >>> multi.registered.coefs.shape, multi.latent.coefs.shape
    ((12, 3, 2), (4, 3))
    """
    data = _as_numpy_fdata(fd)
    _check_finite(data.coefs, "fd")
    like = _torch_reference(fd, None if landmarks is not None else target)
    if landmarks is not None:
        weight = _DEFAULT_LANDMARK_LAMBDA if lam is None else float(lam)
        if weight < 0.0:
            raise ValueError(f"lam must be non-negative, got {weight}")
        result = _landmark(data, landmarks, target_landmarks, warp_basis, weight, penalty)
    else:
        weights = _variable_weights(var_weights, data.n_vars)
        objective: ObjectiveFactory | None = None
        if like is not None:
            from fabel._internal.registration_torch import AutogradObjective

            objective = AutogradObjective
        result = _continuous(
            data,
            target,
            warp_basis,
            0.0 if lam is None else float(lam),
            penalty,
            criterion,
            periodic,
            init,
            init_shift,
            max_iter,
            tol,
            objective,
            weights,
        )
    if like is None:
        return result
    from fabel._internal.registration_torch import to_torch_result

    return to_torch_result(result, fd, target, like, periodic=landmarks is None and periodic)


def _torch_reference(fd: FData, target: Any) -> Array | None:
    """Return the tensor whose dtype and device a torch result takes, if any.

    That is ``fd``'s coefficients when they are a tensor, else the target's
    when the target is an FData with tensor coefficients, else ``None`` (the
    NumPy path).  Checking never imports torch.
    """
    if is_torch(fd.coefs):
        return fd.coefs
    if isinstance(target, FData) and is_torch(target.coefs):
        return target.coefs
    return None


def landmark_register(
    fd: FData,
    landmarks: Any,
    target_landmarks: Any = None,
    **kwargs: Any,
) -> RegistrationResult:
    """Landmark registration; shorthand for ``register(fd, landmarks=...)``.

    Parameters
    ----------
    fd : FData
        The curves to register, univariate or multivariate (one warp per
        curve, shared by its variables).
    landmarks : array
        Landmark times, shape ``(n_curves, n_landmarks)`` or ``(n_curves,)``.
    target_landmarks : array, optional
        Target landmark times; the mean landmark times by default.
    **kwargs
        ``warp_basis``, ``lam`` and ``penalty``, as in :func:`register`.

    Returns
    -------
    RegistrationResult
        As returned by :func:`register`.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import landmark_register
    >>> fd = FData(np.random.default_rng(3).standard_normal((6, 2)), BSpline(n_basis=6))
    >>> landmark_register(fd, [0.3, 0.5]).latent.n_curves
    2
    """
    return register(fd, landmarks=landmarks, target_landmarks=target_landmarks, **kwargs)


# --------------------------------------------------------------------------- #
# scikit-learn estimator
# --------------------------------------------------------------------------- #


class Registrator(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    """Continuous registration as a scikit-learn transformer.

    ``fit`` registers the training curves to their mean and keeps that mean as
    the target; ``transform`` registers new curves to the same target and returns
    the coefficients of the registered curves, one row per curve.

    Parameters
    ----------
    warp_basis : BSpline, optional
        Basis of the latent warp functions, as in :func:`register`.
    lam : float, optional
        Roughness weight on the warps.  Default ``0``.
    penalty : int or LDO, optional
        Roughness operator on the warps.  Default ``2``.
    criterion : {"eigen", "least_squares"}, optional
        Registration criterion.  Default ``"eigen"``.
    periodic : bool, optional
        Estimate a periodic time shift as well.  Default ``False``.
    max_iter : int, optional
        Newton iterations per curve.  Default ``100``.
    tol : float, optional
        Convergence tolerance.  Default ``1e-10``.
    basis : Basis, optional
        Basis to attach to a plain coefficient matrix; ignored for FData input.
        ``None`` builds a cubic B-spline on ``(0, 1)`` with one function per
        column.

    Notes
    -----
    Being a scikit-learn transformer it works on univariate curves (one row
    of coefficients per curve); multivariate FData raise ``ValueError``.

    Attributes
    ----------
    target_ : FData
        The registration target (mean of the training curves).
    result_ : RegistrationResult
        The registration of the training curves.
    n_features_in_ : int
        Number of basis coefficients per curve seen during ``fit``.
    n_iter_ : int
        Largest number of Newton iterations used by a training curve.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline, FData
    >>> from fabel.registration import Registrator
    >>> basis = BSpline(domain=(0.0, 1.0), n_basis=9)
    >>> fd = FData(np.random.default_rng(4).standard_normal((9, 5)), basis)
    >>> Registrator(criterion="least_squares").fit_transform(fd).shape
    (5, 9)
    """

    def __init__(
        self,
        warp_basis: BSpline | None = None,
        *,
        lam: float = 0.0,
        penalty: int | LDO = 2,
        criterion: str = "eigen",
        periodic: bool = False,
        max_iter: int = 100,
        tol: float = 1e-10,
        basis: Basis | None = None,
    ) -> None:
        self.warp_basis = warp_basis
        self.lam = lam
        self.penalty = penalty
        self.criterion = criterion
        self.periodic = periodic
        self.max_iter = max_iter
        self.tol = tol
        self.basis = basis

    def __sklearn_tags__(self) -> Any:
        """Declare a dense, unsupervised transformer."""
        tags = super().__sklearn_tags__()
        tags.target_tags.required = False
        tags.input_tags.sparse = False
        return tags

    def _as_fdata(self, X: Any, *, reset: bool) -> FData:  # noqa: N803
        """Interpret ``X`` as functional data, validating a plain array."""
        if isinstance(X, FData):
            if X.n_vars != 1:
                raise ValueError(
                    "Registrator handles univariate curves; register multivariate "
                    "curves with fabel.registration.register"
                )
            if reset:
                self.n_features_in_ = X.basis.n_basis
            elif X.basis.n_basis != self.n_features_in_:
                raise ValueError(
                    f"X has {X.basis.n_basis} basis functions, expected {self.n_features_in_}"
                )
            return X
        data = validate_data(self, X, reset=reset, ensure_min_samples=1)
        xp = default_namespace()
        basis = (
            self.basis
            if self.basis is not None
            else BSpline(domain=(0.0, 1.0), n_basis=data.shape[1], order=min(4, data.shape[1]))
        )
        if basis.n_basis != data.shape[1]:
            raise ValueError(
                f"X has {data.shape[1]} columns but the basis has {basis.n_basis} functions"
            )
        return FData(xp.matrix_transpose(asarray(data, xp)), basis)

    def _register(self, fd: FData, target: FData) -> RegistrationResult:
        return register(
            fd,
            target,
            warp_basis=self.warp_basis,
            lam=self.lam,
            penalty=self.penalty,
            criterion=self.criterion,
            periodic=self.periodic,
            max_iter=self.max_iter,
            tol=self.tol,
        )

    def fit(self, X: Any, y: Any = None) -> Registrator:  # noqa: N803
        """Register the training curves to their mean.

        Parameters
        ----------
        X : FData or array of shape (n_samples, n_basis)
            Curves, or their basis coefficients one curve per row.
        y : ignored
            Present for API compatibility.

        Returns
        -------
        Registrator
            The fitted estimator.
        """
        fd = self._as_fdata(X, reset=True)
        self.target_ = _as_numpy_fdata(fd).mean()
        self.result_ = self._register(fd, self.target_)
        # A continuous registration always records its per-curve iterations.
        self.n_iter_ = int(to_numpy(self.result_.n_iter).max())
        return self

    def transform(self, X: Any) -> Array:  # noqa: N803
        """Register curves to the fitted target and return their coefficients.

        Parameters
        ----------
        X : FData or array of shape (n_samples, n_basis)
            Curves in the basis seen during ``fit``.

        Returns
        -------
        numpy.ndarray
            Registered coefficients of shape ``(n_samples, n_basis)``.
        """
        check_is_fitted(self)
        fd = self._as_fdata(X, reset=False)
        xp = default_namespace()
        result = self._register(fd, self.target_)
        return xp.matrix_transpose(asarray(to_numpy(result.registered.coefs), xp))
