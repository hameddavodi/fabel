r"""Dynamics: principal differential analysis and phase-plane plots.

Principal differential analysis (PDA, Ramsay & Silverman chapter 19; Ramsay,
Hooker & Graves chapter 11) looks for the linear differential equation that a
sample of curves satisfies as closely as possible.  For a single variable and
order ``m`` it estimates weight functions :math:`\beta_0, \dots, \beta_{m-1}`
such that

.. math::

    L x_n = D^m x_n + \sum_{j=0}^{m-1} \beta_j(t)\, D^j x_n \approx 0
    \qquad (n = 1, \dots, N),

by minimising

.. math::

    \frac{1}{N} \sum_{n=1}^{N} \int (L x_n)(t)^2 \, dt
    + \sum_{j} \lambda_j \int \bigl(K \beta_j\bigr)(t)^2 \, dt ,

where :math:`K` is a roughness operator on the weights.  Expanding each weight
in a basis, :math:`\beta_j = \theta_j^{T} b_j`, makes this a linear least-squares
problem in the stacked coefficients :math:`b`:

.. math::

    \Bigl(A + \operatorname{diag}_j \lambda_j R_j\Bigr) b = r , \quad
    A_{(j,a),(l,c)} = \frac{1}{N} \sum_n \int \theta_{ja}\theta_{lc}\,
        D^j x_n\, D^l x_n , \quad
    r_{(j,a)} = -\frac{1}{N} \sum_n \int \theta_{ja}\, D^j x_n\, D^m x_n .

A system of ``d`` coupled variables is the same problem once per equation
``i``, with the unknowns :math:`\beta_{ikj}` multiplying :math:`D^j x_k` for
every variable ``k``.

Forcing functions :math:`u_1, \dots, u_K` (inputs from outside the system)
enter an equation through weight functions :math:`\alpha_1, \dots,
\alpha_K`,

.. math::

    D^m x = -\sum_{j} \beta_j D^j x + \sum_{k} \alpha_k u_k ,

so the residual becomes :math:`L x - \sum_k \alpha_k u_k` and each
:math:`\alpha_k u_k` is one more regressor, :math:`-u_k`, in the same
least-squares problem (fitted jointly with the :math:`\beta_j`).

The fitted equation, written as a first-order system :math:`Dz = A(t) z +
f(t)` in the state :math:`z = (x, Dx, \dots, D^{m-1} x)`, has a local
stability picture: the eigenvalues of :math:`A(t)` (negative real parts
decay, imaginary parts oscillate) and the equilibrium :math:`z^* = -A(t)^{-1}
f(t)` the system is drawn towards (:meth:`PDA.stability`).

Replaces R's ``pda.fd`` (:class:`PDA`, including ``awtlist``/``ufdlist``),
``eigen.pda`` (:meth:`PDA.stability`), ``pda.overlay``
(:meth:`PDA.plot_overlay`) and ``phaseplanePlot`` (:func:`phase_plane`); the
fitted equation is integrated with :func:`scipy.integrate.solve_ivp`
(:meth:`PDA.solve`).

Examples
--------
>>> import numpy as np
>>> import fdatools as fdt
>>> from fdatools.dynamics import PDA
>>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
>>> curves = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.5], [0.0, 2.0]]), basis)
>>> pda = PDA(order=2, n_grid=None).fit(curves)  # a sin t + b cos t solve D²x + x = 0
>>> [round(float(w.coefs[0, 0]), 10) + 0.0 for w in pda.weights_]
[1.0, 0.0]
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import TYPE_CHECKING, Any

from scipy.integrate import solve_ivp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from fdatools import _linalg
from fdatools._backend import asarray, default_namespace, to_numpy
from fdatools._operator import LDO
from fdatools.basis import Basis, Constant, _same_domain
from fdatools.core import FData, _quadrature

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

__all__ = ["PDA", "PDAStability", "phase_plane"]

Array = Any

#: Points on the equally spaced grid R's ``pda.fd`` integrates with the
#: trapezoidal rule; the default, so that fdatools reproduces ``pda.fd``.
_R_GRID = 501

#: Points used to draw a curve or a weight trajectory.
_PLOT_POINTS = 201

#: Tolerances handed to :func:`scipy.integrate.solve_ivp`.
_ODE_RTOL = 1e-10
_ODE_ATOL = 1e-12

#: Smallest eigenvalue, per unknown, of the unit-diagonal normal matrix below
#: which the PDA problem is reported as singular.
_SINGULAR_TOL = 1e-13

#: Relative slack when checking that evaluation points lie in the domain.
_DOMAIN_TOL = 1e-12


# --------------------------------------------------------------------------- #
# the discretised criterion
# --------------------------------------------------------------------------- #


def _trapezoid_rule(domain: tuple[float, float], n: int, xp: ModuleType) -> tuple[Array, Array]:
    """Return the nodes and weights of the ``n``-point trapezoidal rule on ``domain``."""
    nodes = xp.linspace(domain[0], domain[1], n, dtype=xp.float64)
    step = (domain[1] - domain[0]) / (n - 1)
    weights = xp.full((n,), step, dtype=xp.float64)
    weights[0] = weights[-1] = 0.5 * step
    return nodes, weights


def _derivatives(curves: FData, nodes: Array, order: int) -> list[Array]:
    """Return ``D^j x`` at ``nodes`` for ``j = 0..order``, each ``(n_nodes, N, d)``."""
    out = []
    for j in range(order + 1):
        values = asarray(to_numpy(curves(nodes, j)))
        out.append(values if len(values.shape) == 3 else values[:, :, None])
    return out


def _solve_normal_equations(matrix: Array, rhs: Array, xp: ModuleType) -> Array:
    """Solve the symmetric PDA normal equations, reporting a singular problem.

    The matrix is first scaled to unit diagonal, which removes the (often
    large) differences in scale between ``x`` and its derivatives; the problem
    is singular when the smallest eigenvalue of the scaled matrix is at the
    level of rounding error.
    """
    diagonal = xp.linalg.diagonal(matrix)
    if bool(xp.any(diagonal <= 0.0)):
        raise ValueError(
            "the PDA normal equations are singular: some weight coefficient multiplies a "
            "derivative that vanishes wherever its basis function is non-zero"
        )
    scale = 1.0 / xp.sqrt(diagonal)
    scaled = scale[:, None] * matrix * scale[None, :]
    smallest = float(xp.min(xp.linalg.eigvalsh(0.5 * (scaled + xp.matrix_transpose(scaled)))))
    if smallest < _SINGULAR_TOL * matrix.shape[0]:
        raise ValueError(
            "the PDA normal equations are singular: the derivatives D^j x_k are linearly "
            "dependent over the weight bases (add a roughness penalty or fewer weights)"
        )
    return scale * asarray(_linalg.solve_spd(scaled, scale * rhs))


class _Term:
    """One regressor of an equation: its values, its weight basis and penalty.

    ``values`` is ``(n_nodes, N)``: ``D^j x_k`` for a weight ``β``, ``-u`` for a
    forcing weight ``a``.  ``block`` is the weight basis at the nodes.
    """

    __slots__ = ("block", "lam", "penalty", "values")

    def __init__(self, values: Array, block: Array, penalty: Array | None, lam: float) -> None:
        self.values = values
        self.block = block
        self.penalty = penalty
        self.lam = lam


def _fit_equation(target: Array, terms: list[_Term], quad: Array, xp: ModuleType) -> list[Array]:
    """Estimate the weight coefficients of one equation.

    ``target`` is ``D^m x_i`` of shape ``(n_nodes, N)``.  The residual is
    ``target + Σ_p w_p * values_p`` with ``w_p = block_p @ coefs_p``; the
    coefficients minimise its mean integrated square plus the penalties.
    Returns one coefficient vector per term.
    """
    n_curves = target.shape[1]
    sizes = [term.block.shape[1] for term in terms]
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1] + size)
    matrix = xp.zeros((offsets[-1], offsets[-1]), dtype=xp.float64)
    rhs = xp.zeros((offsets[-1],), dtype=xp.float64)
    for p, term in enumerate(terms):
        left = term.values
        for q, other in enumerate(terms[: p + 1]):
            kernel = xp.sum(left * other.values, axis=1) * quad / n_curves
            block = xp.matmul(xp.matrix_transpose(term.block), kernel[:, None] * other.block)
            matrix[offsets[p] : offsets[p + 1], offsets[q] : offsets[q + 1]] = block
            matrix[offsets[q] : offsets[q + 1], offsets[p] : offsets[p + 1]] = xp.matrix_transpose(
                block
            )
        forcing = xp.sum(left * target, axis=1) * quad / n_curves
        rhs[offsets[p] : offsets[p + 1]] = -xp.matmul(xp.matrix_transpose(term.block), forcing)
        if term.penalty is not None:
            matrix[offsets[p] : offsets[p + 1], offsets[p] : offsets[p + 1]] += (
                term.lam * term.penalty
            )
    solution = _solve_normal_equations(matrix, rhs, xp)
    return [solution[offsets[p] : offsets[p + 1]] for p in range(len(terms))]


def _forcing_values(forcing: list[list[FData]], nodes: Array, n_curves: int) -> list[list[Array]]:
    """Evaluate every forcing function at ``nodes`` as ``(n_nodes, n_curves)``.

    A forcing function with a single curve is shared by all ``n_curves``
    curves.
    """
    xp = default_namespace()
    out: list[list[Array]] = []
    for functions in forcing:
        row = []
        for u in functions:
            values = asarray(to_numpy(u(nodes)))
            row.append(xp.broadcast_to(values, (values.shape[0], n_curves)))
        out.append(row)
    return out


# --------------------------------------------------------------------------- #
# plotting helpers
# --------------------------------------------------------------------------- #


def _axes(ax: Axes | None) -> Axes:
    """Return ``ax``, or the axes of a new figure."""
    if ax is not None:
        return ax
    import matplotlib.pyplot as plt

    _, new = plt.subplots()
    return new


def _annotate(ax: Axes, points: Mapping[float, str], x_of: Any, y_of: Any) -> None:
    """Write each label at the position the curve reaches at its time point."""
    for when, text in points.items():
        ax.text(float(x_of(when)), float(y_of(when)), str(text))


# --------------------------------------------------------------------------- #
# stability analysis
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PDAStability:
    r"""Local stability of a fitted differential equation over time.

    Replaces the value of R's ``eigen.pda``.  At each time ``t`` the fitted
    equations are written as the first-order system ``Dz = A(t) z + f(t)`` in
    the state ``z = (x_1, Dx_1, …, D^{m-1}x_1, x_2, …)``, where ``A(t)`` is the
    companion matrix of the weights ``β`` and ``f(t)`` holds the forcing terms
    ``Σ_k a_k(t) u_k(t)``.

    Attributes
    ----------
    t : array
        The ``(n_t,)`` time points.
    eigenvalues : array
        ``(n_t, n_vars * order)`` complex eigenvalues of ``A(t)``, each row in
        decreasing modulus (a complex pair keeps its positive imaginary part
        first).  A negative real part means local exponential decay, a
        positive one growth, and a non-zero imaginary part oscillation with
        angular frequency ``|Im λ|``.
    limits : array
        ``(n_t, n_vars * order)`` equilibrium states ``z* = -A(t)⁻¹ f(t)``, the
        point the system is drawn to (when it is stable) if the weights and
        forcing were frozen at their value at ``t``.  Zero for an unforced
        equation and ``nan`` where ``A(t)`` is exactly singular.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.dynamics import PDA
    >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
    >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
    >>> result = PDA(order=2, n_grid=None).fit(fd).stability(n_points=3)
    >>> np.round(result.eigenvalues[0], 10) + 0.0
    array([0.+1.j, 0.-1.j])
    """

    t: Array
    eigenvalues: Array
    limits: Array

    def plot(self, ax: Axes | None = None, **kwargs: Any) -> Axes:
        """Plot the real (solid) and imaginary (dashed) parts of the eigenvalues.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when omitted.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.plot` for every line.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on: one solid line per eigenvalue for its real part,
            then one dashed line per eigenvalue for its imaginary part.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> len(PDA(order=2).fit(fd).stability().plot().lines)
        5
        """
        axes = _axes(ax)
        times = to_numpy(self.t)
        values = to_numpy(self.eigenvalues)
        for column in range(values.shape[1]):
            axes.plot(times, values[:, column].real, **kwargs)
        for column in range(values.shape[1]):
            axes.plot(times, values[:, column].imag, linestyle="--", **kwargs)
        axes.axhline(0.0, color="grey", linewidth=0.5)
        axes.set_xlabel("t")
        axes.set_ylabel("eigenvalue (real solid, imaginary dashed)")
        return axes


def _is_basis_spec(value: Any) -> bool:
    return value is None or isinstance(value, Basis)


def _is_lambda_spec(value: Any) -> bool:
    return isinstance(value, (int, float))


def _per_forcing(spec: Any, counts: list[int], leaf: Any, name: str) -> list[list[Any]]:
    """Spread ``forcing_basis``/``forcing_lam`` over the forcing functions.

    A single leaf applies to every forcing function.  Otherwise ``spec``
    mirrors the forcing: one entry per forcing function for one equation; for
    a system one entry per equation, each a leaf or one entry per forcing
    function of that equation.
    """
    if leaf(spec):
        return [[spec] * count for count in counts]
    items = list(spec)
    if len(counts) == 1:
        groups = [items]
    else:
        if len(items) != len(counts):
            raise ValueError(f"{name} must hold {len(counts)} entries (one per equation)")
        groups = [
            [entry] * count if leaf(entry) else list(entry)
            for entry, count in zip(items, counts, strict=True)
        ]
    for group, count in zip(groups, counts, strict=True):
        if len(group) != count:
            raise ValueError(f"{name} must hold {count} entries for an equation, got {len(group)}")
        if not all(leaf(entry) for entry in group):
            raise TypeError(f"{name} holds an entry of the wrong type")
    return groups


# --------------------------------------------------------------------------- #
# the estimator
# --------------------------------------------------------------------------- #


class PDA(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    r"""Principal differential analysis: fit a linear ODE to a sample of curves.

    Replaces R's ``pda.fd``.  Estimates weight functions ``β_j`` so that
    ``D^m x + Σ_j β_j D^j x`` is as small as possible over the sample, in the
    mean integrated square, with an optional roughness penalty on each weight
    (see the module docstring for the criterion).  Curves with several
    variables (``n_vars = d > 1``) define a system of ``d`` coupled equations.
    Forcing functions ``u_k`` passed to :meth:`fit` add terms ``a_k u_k`` to
    the right-hand side, ``D^m x = -Σ_j β_j D^j x + Σ_k a_k u_k`` (R's
    ``awtlist``/``ufdlist``); their weights ``a_k`` are estimated jointly with
    the ``β_j``.

    Parameters
    ----------
    order : int, optional
        Order ``m`` of the equation.  Default ``2``.
    weight_basis : Basis or sequence of Basis, optional
        Basis for the weight functions; one per derivative order ``j`` when a
        sequence of length ``order``.  ``None`` (default) uses a constant
        basis, i.e. a constant-coefficient equation.
    lam : float or sequence of float, optional
        Roughness penalty on each weight function, per derivative order when a
        sequence.  Default ``0.0``.
    penalty : int or LDO, optional
        Roughness operator applied to the weights (``β`` and ``a``).  Default
        ``2``.
    n_grid : int or None, optional
        How the integrals are computed.  An integer is the number of equally
        spaced points of a trapezoidal rule; the default ``501`` is the rule
        R's ``pda.fd`` uses, so fdatools reproduces it.  (R raises its grid to
        five times the number of basis functions of the curves when that is
        larger; pass that number to reproduce R for curve bases of more than
        100 functions.)  The residual functions are then least-squares fits
        of the residual on the same grid.  ``None`` integrates exactly
        (composite Gauss-Legendre on the break points, exact for spline
        curves and weights) and projects the residuals in ``L²``; it differs
        from R's rule by about ``1e-5`` relative on realistic data.
    forcing_basis : Basis or sequence, optional
        Basis for the forcing weights ``a``.  ``None`` (default) is a constant
        basis.  A single basis applies to every forcing function; otherwise one
        entry per forcing function (for a system: one entry per equation, each
        a basis or one basis per forcing function of that equation).
    forcing_lam : float or sequence, optional
        Roughness penalty on each forcing weight, spread like
        ``forcing_basis``.  Default ``0.0``.

    Attributes
    ----------
    weights_ : tuple
        For one variable, ``weights_[j]`` is ``β_j`` as a one-curve
        :class:`~fdatools.core.FData`.  For ``d`` variables, ``weights_[i][k][j]``
        multiplies ``D^j x_k`` in equation ``i``.  This mirrors the nesting of
        R's ``bwtlist``.
    forcing_weights_ : tuple
        For one variable, ``forcing_weights_[k]`` is ``a_k``; for ``d``
        variables, ``forcing_weights_[i][k]`` multiplies forcing function ``k``
        of equation ``i`` (R's ``awtlist``).  Empty when fitted without
        forcing.
    residuals_ : FData
        ``L x - Σ_k a_k u_k`` for the fitted curves, in their basis.
    operator_ : LDO or None
        The fitted operator ``L`` for one variable (so that ``fd(t,
        pda.operator_)`` evaluates ``L x``); ``None`` for a system.
    n_vars_ : int
        Number of variables of the fitted curves.
    domain_ : tuple of float
        Domain of the fitted curves.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.dynamics import PDA
    >>> t = np.linspace(0.0, 1.0, 101)
    >>> basis = fdt.BSpline(domain=(0.0, 1.0), n_basis=24, order=5)
    >>> from fdatools.smoothing import smooth
    >>> fd = smooth(np.exp(-4 * t), t, basis=basis, lam=0.0).fd
    >>> pda = PDA(order=1).fit(fd)  # Dx + βx = 0
    >>> round(float(pda.weights_[0].coefs[0, 0]), 6)
    4.0

    A forcing function: ``x = 0.5 (1 - exp(-4t))`` solves ``Dx = -4x + 2u``
    for the constant input ``u = 1``.

    >>> fd = smooth(0.5 * (1 - np.exp(-4 * t)), t, basis=basis, lam=0.0).fd
    >>> u = fdt.FData(np.array([1.0]), fdt.Constant(domain=(0.0, 1.0)))
    >>> pda = PDA(order=1).fit(fd, forcing=u)
    >>> (
    ...     round(float(pda.weights_[0].coefs[0, 0]), 6),
    ...     round(float(pda.forcing_weights_[0].coefs[0, 0]), 6),
    ... )
    (4.0, 2.0)
    """

    def __init__(
        self,
        order: int = 2,
        *,
        weight_basis: Basis | Sequence[Basis] | None = None,
        lam: float | Sequence[float] = 0.0,
        penalty: int | LDO = 2,
        n_grid: int | None = _R_GRID,
        forcing_basis: Basis | Sequence[Any] | None = None,
        forcing_lam: float | Sequence[Any] = 0.0,
    ) -> None:
        self.order = order
        self.weight_basis = weight_basis
        self.lam = lam
        self.penalty = penalty
        self.n_grid = n_grid
        self.forcing_basis = forcing_basis
        self.forcing_lam = forcing_lam

    def __sklearn_tags__(self) -> Any:
        """Declare an unsupervised transformer of functional data."""
        tags = super().__sklearn_tags__()
        tags.target_tags.required = False
        tags.input_tags.sparse = False
        return tags

    # ------------------------------------------------------------ validation

    def _checked_order(self) -> int:
        order = int(self.order)
        if order < 1:
            raise ValueError(f"order must be at least 1, got {self.order}")
        return order

    def _weight_bases(self, domain: tuple[float, float], order: int) -> list[Basis]:
        spec = self.weight_basis
        if spec is None:
            return [Constant(domain=domain)] * order
        bases = [spec] * order if isinstance(spec, Basis) else list(spec)
        if len(bases) != order:
            raise ValueError(f"weight_basis must hold {order} bases, got {len(bases)}")
        for basis in bases:
            if not _same_domain(basis.domain, domain):
                raise ValueError(
                    f"weight basis domain {basis.domain} differs from the curves' domain {domain}"
                )
        return bases

    def _lambdas(self, order: int) -> list[float]:
        spec = self.lam
        lams = [float(spec)] * order if isinstance(spec, (int, float)) else [float(v) for v in spec]
        if len(lams) != order:
            raise ValueError(f"lam must hold {order} values, got {len(lams)}")
        if any(value < 0.0 for value in lams):
            raise ValueError(f"lam must be non-negative, got {lams}")
        return lams

    def _forcing_bases(self, counts: list[int], domain: tuple[float, float]) -> list[list[Basis]]:
        groups = _per_forcing(self.forcing_basis, counts, _is_basis_spec, "forcing_basis")
        out: list[list[Basis]] = []
        for group in groups:
            row = []
            for basis in group:
                chosen = Constant(domain=domain) if basis is None else basis
                if not _same_domain(chosen.domain, domain):
                    raise ValueError(
                        f"forcing basis domain {chosen.domain} differs from the curves' "
                        f"domain {domain}"
                    )
                row.append(chosen)
            out.append(row)
        return out

    def _forcing_lambdas(self, counts: list[int]) -> list[list[float]]:
        groups = _per_forcing(self.forcing_lam, counts, _is_lambda_spec, "forcing_lam")
        lams = [[float(value) for value in group] for group in groups]
        if any(value < 0.0 for group in lams for value in group):
            raise ValueError(f"forcing_lam must be non-negative, got {lams}")
        return lams

    @staticmethod
    def _forcing_functions(
        forcing: Any, n_vars: int, domain: tuple[float, float], n_curves: int | None
    ) -> list[list[FData]]:
        """Normalise ``forcing`` to one list of forcing functions per equation.

        One equation takes an :class:`FData` or a sequence of them; a system
        takes one entry per equation, each ``None``, an :class:`FData` or a
        sequence of them.  Each forcing function has one variable, the curves'
        domain, and one curve or ``n_curves`` curves (``n_curves=None`` asks
        for exactly one curve).
        """
        if forcing is None:
            return [[] for _ in range(n_vars)]
        if isinstance(forcing, FData):
            if n_vars > 1:
                raise ValueError(
                    f"a system of {n_vars} equations takes one forcing entry per equation"
                )
            groups: list[list[Any]] = [[forcing]]
        elif n_vars == 1:
            groups = [list(forcing)]
        else:
            entries = list(forcing)
            if len(entries) != n_vars:
                raise ValueError(
                    f"forcing must hold {n_vars} entries (one per equation), got {len(entries)}"
                )
            groups = [
                [] if entry is None else [entry] if isinstance(entry, FData) else list(entry)
                for entry in entries
            ]
        for group in groups:
            for u in group:
                if not isinstance(u, FData):
                    raise TypeError(f"forcing functions must be FData, got {type(u).__name__}")
                if u.n_vars != 1:
                    raise ValueError(f"a forcing function has one variable, got {u.n_vars}")
                if not _same_domain(u.domain, domain):
                    raise ValueError(
                        f"forcing function domain {u.domain} differs from the curves' "
                        f"domain {domain}"
                    )
                allowed = (1,) if n_curves is None else (1, n_curves)
                if u.n_curves not in allowed:
                    wanted = "1" if n_curves is None else f"1 or {n_curves}"
                    raise ValueError(
                        f"a forcing function must have {wanted} curves, got {u.n_curves}"
                    )
        return groups

    def _matched_forcing(self, forcing: Any, n_curves: int | None) -> list[list[FData]]:
        """Return the forcing functions for a fitted model, checked against the fit."""
        groups = self._forcing_functions(forcing, self.n_vars_, self.domain_, n_curves)
        counts = [len(group) for group in groups]
        if counts != self._counts:
            raise ValueError(
                f"the model was fitted with {self._counts} forcing functions per equation, "
                f"got {counts}"
            )
        return groups

    def _rule(self, curves: FData, bases: list[Basis], xp: ModuleType) -> tuple[Array, Array]:
        if self.n_grid is None:
            nodes, weights = _quadrature(curves.basis, *bases)
            return asarray(nodes), asarray(weights)
        n_grid = int(self.n_grid)
        if n_grid < 2:
            raise ValueError(f"n_grid must be at least 2, got {self.n_grid}")
        if n_grid < curves.basis.n_basis:
            raise ValueError(
                f"n_grid={n_grid} cannot resolve a basis of {curves.basis.n_basis} functions"
            )
        return _trapezoid_rule(curves.domain, n_grid, xp)

    # --------------------------------------------------------------- fitting

    def fit(self, X: FData, y: Any = None, *, forcing: Any = None) -> PDA:  # noqa: N803
        """Estimate the weight functions from a sample of curves.

        Parameters
        ----------
        X : FData
            The curves; a multivariate :class:`~fdatools.core.FData` (``n_vars >
            1``) is fitted as a system of coupled equations.
        y : None
            Ignored; present for scikit-learn compatibility.
        forcing : FData or sequence, optional
            Forcing functions ``u_k`` (R's ``ufdlist``).  For one equation an
            :class:`~fdatools.core.FData` or a sequence of them; for a system one
            entry per equation, each ``None``, an :class:`~fdatools.core.FData` or
            a sequence of them.  Each has one variable and either one curve per
            curve of ``X`` or a single curve shared by all of them.

        Returns
        -------
        PDA
            The fitted estimator.

        Raises
        ------
        TypeError
            If ``X`` or a forcing function is not an :class:`~fdatools.core.FData`.
        ValueError
            If a parameter is invalid or the normal equations are singular.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> len(PDA(order=2).fit(fd).weights_)
        2
        """
        if not isinstance(X, FData):
            raise TypeError(f"PDA.fit expects an FData, got {type(X).__name__}")
        xp = default_namespace()
        order = self._checked_order()
        bases = self._weight_bases(X.domain, order)
        lams = self._lambdas(order)
        n_vars = X.n_vars
        forcing_fns = self._forcing_functions(forcing, n_vars, X.domain, X.n_curves)
        counts = [len(group) for group in forcing_fns]
        forcing_bases = self._forcing_bases(counts, X.domain)
        forcing_lams = self._forcing_lambdas(counts)
        penalties: list[Array | None] = [
            asarray(basis.penalty(self.penalty)) if lam > 0.0 else None
            for basis, lam in zip(bases, lams, strict=True)
        ]
        extra = [b for group in forcing_bases for b in group]
        extra += [u.basis for group in forcing_fns for u in group]
        nodes, quad = self._rule(X, bases + extra, xp)
        derivs = _derivatives(X, nodes, order)
        blocks = [asarray(basis(nodes)) for basis in bases]
        u_values = _forcing_values(forcing_fns, nodes, X.n_curves)
        weights: list[list[list[FData]]] = []
        alphas: list[list[FData]] = []
        for i in range(n_vars):
            terms = [
                _Term(derivs[j][:, :, k], blocks[j], penalties[j], lams[j])
                for k in range(n_vars)
                for j in range(order)
            ]
            for basis, lam, values in zip(
                forcing_bases[i], forcing_lams[i], u_values[i], strict=True
            ):
                pen = asarray(basis.penalty(self.penalty)) if lam > 0.0 else None
                terms.append(_Term(-values, asarray(basis(nodes)), pen, lam))
            coefs = _fit_equation(derivs[order][:, :, i], terms, quad, xp)
            weights.append(
                [
                    [FData(coefs[k * order + j], bases[j]) for j in range(order)]
                    for k in range(n_vars)
                ]
            )
            alphas.append(
                [
                    FData(coefs[n_vars * order + position], basis)
                    for position, basis in enumerate(forcing_bases[i])
                ]
            )
        self.n_vars_ = n_vars
        self.domain_ = X.domain
        if n_vars == 1:
            self.weights_: tuple[Any, ...] = tuple(weights[0][0])
            self.operator_: LDO | None = LDO(weights=list(self.weights_))
            self.forcing_weights_: tuple[Any, ...] = tuple(alphas[0])
        else:
            self.weights_ = tuple(tuple(tuple(row) for row in eq) for eq in weights)
            self.operator_ = None
            self.forcing_weights_ = tuple(tuple(row) for row in alphas) if any(counts) else ()
        self._nested = weights
        self._alphas = alphas
        self._counts = counts
        self._forcing = forcing_fns
        self.residuals_ = self._residuals(X, nodes, quad, derivs, u_values)
        return self

    def _residuals(
        self,
        curves: FData,
        nodes: Array,
        quad: Array,
        derivs: list[Array],
        u_values: list[list[Array]],
    ) -> FData:
        """Return ``L x - Σ a u`` for ``curves`` expressed in their own basis."""
        xp = default_namespace()
        order = len(derivs) - 1
        n_vars = self.n_vars_
        columns = []
        for i in range(n_vars):
            value = derivs[order][:, :, i]
            for k in range(n_vars):
                for j in range(order):
                    beta = asarray(to_numpy(self._nested[i][k][j](nodes)))
                    value = value + beta * derivs[j][:, :, k]
            for alpha, u in zip(self._alphas[i], u_values[i], strict=True):
                value = value - asarray(to_numpy(alpha(nodes))) * u
            columns.append(value)
        values = xp.stack(columns, axis=2)
        design = asarray(curves.basis(nodes))
        flat = xp.reshape(values, (values.shape[0], -1))
        if self.n_grid is None:
            gram = xp.matmul(xp.matrix_transpose(design), quad[:, None] * design)
            rhs = xp.matmul(xp.matrix_transpose(design), quad[:, None] * flat)
            coefs = asarray(_linalg.solve_spd(gram, rhs))
        else:
            coefs = xp.linalg.lstsq(design, flat, rcond=None)[0]
        shape = (curves.basis.n_basis, values.shape[1], n_vars)
        coefs = xp.reshape(coefs, shape)
        return FData(coefs if n_vars > 1 else coefs[:, :, 0], curves.basis)

    def transform(self, X: FData, *, forcing: Any = None) -> FData:  # noqa: N803
        """Apply the fitted equation to curves: return ``L x - Σ a u``.

        Parameters
        ----------
        X : FData
            Curves with the same number of variables and domain as the fitted
            ones.
        forcing : FData or sequence, optional
            Forcing functions for ``X``, structured as in :meth:`fit`; required
            exactly when the model was fitted with forcing.

        Returns
        -------
        FData
            The residual functions, in the basis of ``X``.

        Raises
        ------
        TypeError
            If ``X`` is not an :class:`~fdatools.core.FData`.
        ValueError
            If ``X`` or ``forcing`` does not match the fitted model.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> pda = PDA(order=2, n_grid=None).fit(fd)
        >>> new = fdt.FData(np.array([0.0, 2.0, -1.0]), basis)
        >>> bool(np.max(np.abs(pda.transform(new)(np.linspace(0, 6, 7)))) < 1e-10)
        True
        """
        check_is_fitted(self, "weights_")
        if not isinstance(X, FData):
            raise TypeError(f"PDA.transform expects an FData, got {type(X).__name__}")
        if X.n_vars != self.n_vars_:
            raise ValueError(f"X has {X.n_vars} variables, the model was fitted on {self.n_vars_}")
        if not _same_domain(X.domain, self.domain_):
            raise ValueError(f"X has domain {X.domain}, the model was fitted on {self.domain_}")
        forcing_fns = self._matched_forcing(forcing, X.n_curves)
        xp = default_namespace()
        order = self._checked_order()
        bases = [w.basis for w in self._nested[0][0]]
        bases += [a.basis for group in self._alphas for a in group]
        bases += [u.basis for group in forcing_fns for u in group]
        nodes, quad = self._rule(X, bases, xp)
        u_values = _forcing_values(forcing_fns, nodes, X.n_curves)
        return self._residuals(X, nodes, quad, _derivatives(X, nodes, order), u_values)

    def fit_transform(self, X: FData, y: Any = None, *, forcing: Any = None) -> FData:  # noqa: N803
        """Fit to ``X`` (with its forcing) and return the residual functions.

        Parameters
        ----------
        X : FData
            The curves.
        y : None
            Ignored; present for scikit-learn compatibility.
        forcing : FData or sequence, optional
            Forcing functions, as in :meth:`fit`.

        Returns
        -------
        FData
            ``L x - Σ a u`` for the fitted curves (equal to ``residuals_``).

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> PDA(order=2).fit_transform(fd).n_curves
        2
        """
        return self.fit(X, y, forcing=forcing).transform(X, forcing=forcing)

    # --------------------------------------------------------------- solving

    def _checked_points(self, t: Any) -> Array:
        """Validate evaluation points: one-dimensional, non-empty, in the domain."""
        xp = default_namespace()
        points = asarray(to_numpy(t))
        if len(points.shape) != 1 or points.shape[0] == 0:
            raise ValueError("t must be a one-dimensional array with at least one point")
        lower, upper = self.domain_
        slack = _DOMAIN_TOL * max(1.0, upper - lower)
        if float(xp.min(points)) < lower - slack or float(xp.max(points)) > upper + slack:
            raise ValueError(f"t must lie in the domain {self.domain_}")
        return points

    def solve(
        self,
        t: Any,
        initial: Any,
        *,
        forcing: Any = None,
        rtol: float = _ODE_RTOL,
        atol: float = _ODE_ATOL,
    ) -> Array:
        """Integrate the fitted equation from an initial state.

        The equation is rewritten as a first-order system in ``(x, Dx, ...,
        D^{m-1} x)`` and integrated with :func:`scipy.integrate.solve_ivp`
        (method ``DOP853``) from ``t[0]``, forwards or backwards.  Without
        ``forcing`` this is the homogeneous equation ``L x = 0`` (the free
        response, also for a model fitted with forcing); with it, ``L x =
        Σ a_k u_k``.

        Parameters
        ----------
        t : array_like
            Strictly monotone points inside the domain; the initial condition
            holds at ``t[0]``.
        initial : array_like
            ``(x, Dx, ..., D^{m-1} x)`` at ``t[0]``: shape ``(order,)`` for one
            variable, ``(n_vars, order)`` for a system.
        forcing : FData or sequence, optional
            The input functions, structured as in :meth:`fit`, one curve each.
        rtol, atol : float, optional
            Relative and absolute tolerances of the integrator.

        Returns
        -------
        numpy.ndarray
            ``x(t)``: shape ``(n_t,)`` for one variable, ``(n_t, n_vars)`` for a
            system.

        Raises
        ------
        ValueError
            If ``t`` is empty, not monotone or outside the domain, if
            ``initial`` has the wrong shape, if ``forcing`` does not match the
            fitted model, or if the integrator fails.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> pda = PDA(order=2, n_grid=None).fit(fd)
        >>> t = np.linspace(0.0, np.pi, 5)
        >>> bool(np.allclose(pda.solve(t, [0.0, 1.0]), np.sin(t), atol=1e-8))
        True
        """
        check_is_fitted(self, "weights_")
        xp = default_namespace()
        order = self._checked_order()
        n_vars = self.n_vars_
        points = self._checked_points(t)
        steps = points[1:] - points[:-1]
        if points.shape[0] > 1 and not (bool(xp.all(steps > 0.0)) or bool(xp.all(steps < 0.0))):
            raise ValueError("t must be strictly monotone")
        start = asarray(to_numpy(initial))
        expected = (order,) if n_vars == 1 else (n_vars, order)
        if tuple(start.shape) != expected:
            raise ValueError(f"initial must have shape {expected}, got {tuple(start.shape)}")
        inputs = (
            [[] for _ in range(n_vars)] if forcing is None else self._matched_forcing(forcing, None)
        )
        state0 = xp.reshape(start, (n_vars * order,))
        if points.shape[0] == 1:
            first = xp.reshape(state0, (n_vars, order))[:, 0]
            return first if n_vars == 1 else first[None, :]
        clip_lo, clip_hi = self.domain_
        nested = self._nested
        alphas = self._alphas

        def rhs(time: float, state: Array) -> Array:
            when = xp.asarray([min(max(time, clip_lo), clip_hi)], dtype=xp.float64)
            z = xp.reshape(asarray(state), (n_vars, order))
            dz = xp.zeros((n_vars, order), dtype=xp.float64)
            if order > 1:
                dz[:, :-1] = z[:, 1:]
            for i in range(n_vars):
                total = 0.0
                for k in range(n_vars):
                    for j in range(order):
                        total += float(nested[i][k][j](when)[0, 0]) * float(z[k, j])
                drive = 0.0
                for alpha, u in zip(alphas[i], inputs[i], strict=False):
                    drive += float(alpha(when)[0, 0]) * float(u(when)[0, 0])
                dz[i, -1] = drive - total
            return xp.reshape(dz, (n_vars * order,))

        span = (float(points[0]), float(points[-1]))
        result = solve_ivp(
            rhs,
            span,
            to_numpy(state0),
            method="DOP853",
            t_eval=to_numpy(points),
            rtol=rtol,
            atol=atol,
        )
        if not result.success:  # pragma: no cover - DOP853 on a linear ODE does not fail
            raise ValueError(f"the ODE integration failed: {result.message}")
        values = asarray(result.y)
        positions = xp.reshape(values, (n_vars, order, points.shape[0]))[:, 0, :]
        solution = xp.matrix_transpose(positions)
        return solution[:, 0] if n_vars == 1 else solution

    # ------------------------------------------------------------- stability

    def stability(
        self, t: Any = None, *, n_points: int = _R_GRID, forcing: Any = None
    ) -> PDAStability:
        """Eigenvalues and equilibria of the fitted equation over time.

        Replaces R's ``eigen.pda``.  At each time the equations are the
        first-order system ``Dz = A(t) z + f(t)`` in the state ``z = (x_1,
        Dx_1, …, D^{m-1} x_1, x_2, …)``: ``A(t)`` has ones on the
        superdiagonal of each variable's block and ``-β_{ikj}(t)`` in the last
        row of block ``i``, column ``(k, j)``; ``f(t)`` holds ``Σ_k a_k(t)
        u_k(t)`` in the last row of each block.  See :class:`PDAStability`.

        Parameters
        ----------
        t : array_like, optional
            Time points in the domain.  Defaults to ``n_points`` equally spaced
            points spanning it (R's ``argvals``).
        n_points : int, optional
            Number of default time points.  Default ``501``, as in R.
        forcing : FData or sequence, optional
            The input functions, structured as in :meth:`fit`, one curve each.
            Defaults to the forcing the model was fitted with when every one of
            those functions has a single curve.

        Returns
        -------
        PDAStability
            Eigenvalues and equilibrium states at every time point.

        Raises
        ------
        ValueError
            If ``t`` is invalid, or the forcing is missing or does not match
            the fitted model.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> from fdatools.smoothing import smooth
        >>> t = np.linspace(0.0, 1.0, 101)
        >>> basis = fdt.BSpline(domain=(0.0, 1.0), n_basis=24, order=5)
        >>> fd = smooth(0.5 * (1 - np.exp(-4 * t)), t, basis=basis, lam=0.0).fd
        >>> u = fdt.FData(np.array([1.0]), fdt.Constant(domain=(0.0, 1.0)))
        >>> result = PDA(order=1).fit(fd, forcing=u).stability(n_points=3)
        >>> np.round(result.eigenvalues.real, 6) + 0.0, np.round(result.limits, 6) + 0.0
        (array([[-4.],
               [-4.],
               [-4.]]), array([[0.5],
               [0.5],
               [0.5]]))
        """
        check_is_fitted(self, "weights_")
        xp = default_namespace()
        order = self._checked_order()
        n_vars = self.n_vars_
        if t is None:
            if int(n_points) < 1:
                raise ValueError(f"n_points must be at least 1, got {n_points}")
            points = xp.linspace(self.domain_[0], self.domain_[1], int(n_points), dtype=xp.float64)
        else:
            points = self._checked_points(t)
        forced = any(self._counts)
        if forcing is not None:
            inputs = self._matched_forcing(forcing, None)
        elif forced:
            if any(u.n_curves != 1 for group in self._forcing for u in group):
                raise ValueError(
                    "the model was fitted with forcing functions of several curves; pass "
                    "forcing= with one curve per forcing function"
                )
            inputs = self._forcing
        else:
            inputs = [[] for _ in range(n_vars)]
        size = n_vars * order
        n_t = points.shape[0]
        matrix = xp.zeros((n_t, size, size), dtype=xp.float64)
        drive = xp.zeros((n_t, size), dtype=xp.float64)
        for i in range(n_vars):
            for j in range(order - 1):
                matrix[:, i * order + j, i * order + j + 1] = 1.0
            last = i * order + order - 1
            for k in range(n_vars):
                for j in range(order):
                    beta = asarray(to_numpy(self._nested[i][k][j](points)))[:, 0]
                    matrix[:, last, k * order + j] = -beta
            for alpha, u in zip(self._alphas[i], inputs[i], strict=True):
                a_values = asarray(to_numpy(alpha(points)))[:, 0]
                drive[:, last] += a_values * asarray(to_numpy(u(points)))[:, 0]
        values = xp.linalg.eigvals(matrix)
        ranking = xp.argsort(-xp.abs(values), axis=1, stable=True)
        values = xp.take_along_axis(values, ranking, axis=1)
        if not any(inputs):
            limits = xp.zeros((n_t, size), dtype=xp.float64)
        else:
            regular = xp.linalg.det(matrix) != 0.0
            identity = xp.eye(size, dtype=xp.float64)
            safe = xp.where(regular[:, None, None], matrix, identity)
            limits = -xp.linalg.solve(safe, drive[:, :, None])[:, :, 0]
            limits = xp.where(regular[:, None], limits, xp.nan) + 0.0
        return PDAStability(t=points, eigenvalues=values, limits=limits)

    # --------------------------------------------------------------- plotting

    def plot_overlay(
        self,
        ax: Axes | None = None,
        *,
        n_points: int = _PLOT_POINTS,
        labels: Mapping[float, str] | None = None,
        **kwargs: Any,
    ) -> Axes:
        r"""Plot the weights of a second-order equation on its stability diagram.

        Replaces R's ``pda.overlay``.  The trajectory ``(β₁(t), β₀(t))`` is drawn
        over the parabola ``β₀ = β₁² / 4``: above it the characteristic roots of
        ``D²x + β₁ Dx + β₀ x = 0`` are complex (the system oscillates), below it
        they are real; ``β₁ > 0`` damps and ``β₁ < 0`` amplifies.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when omitted.
        n_points : int, optional
            Number of time points on the trajectory.
        labels : mapping of float to str, optional
            Text to write at the trajectory's position at given times.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.plot` for the trajectory.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on; ``lines[0]`` is the trajectory and ``lines[1]``
            the parabola.

        Raises
        ------
        ValueError
            If the model is not a single second-order equation.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.dynamics import PDA
        >>> basis = fdt.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fdt.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> len(PDA(order=2).fit(fd).plot_overlay().lines) >= 2
        True
        """
        check_is_fitted(self, "weights_")
        if self._checked_order() != 2 or self.n_vars_ != 1:
            raise ValueError("plot_overlay needs a single second-order equation")
        xp = default_namespace()
        grid = xp.linspace(self.domain_[0], self.domain_[1], n_points, dtype=xp.float64)
        beta0 = to_numpy(self.weights_[0](grid))[:, 0]
        beta1 = to_numpy(self.weights_[1](grid))[:, 0]
        axes = _axes(ax)
        axes.plot(beta1, beta0, **kwargs)
        reach = max(float(abs(beta1).max()), 1e-12) * 1.1
        edge = to_numpy(xp.linspace(-reach, reach, _PLOT_POINTS, dtype=xp.float64))
        axes.plot(edge, edge**2 / 4.0, linestyle="--", color="grey")
        axes.axhline(0.0, color="grey", linewidth=0.5)
        axes.axvline(0.0, color="grey", linewidth=0.5)
        if labels:

            def weight_at(j: int) -> Any:
                return lambda when: to_numpy(
                    self.weights_[j](xp.asarray([when], dtype=xp.float64))
                )[0, 0]

            _annotate(axes, labels, weight_at(1), weight_at(0))
        axes.set_xlabel("beta 1")
        axes.set_ylabel("beta 0")
        return axes


# --------------------------------------------------------------------------- #
# phase-plane plot
# --------------------------------------------------------------------------- #


def phase_plane(
    fd: FData,
    t: Any = None,
    *,
    deriv: tuple[int, int] = (1, 2),
    labels: Mapping[float, str] | None = None,
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes:
    """Plot one derivative of each curve against another.

    Replaces R's ``phaseplanePlot``.  The default plots acceleration ``D²x``
    against velocity ``Dx``: kinetic energy grows with the horizontal distance
    from the origin and potential energy with the vertical one, so an
    oscillation traces a loop around the origin.

    Parameters
    ----------
    fd : FData
        Single-variable curves.
    t : array_like, optional
        Evaluation points.  Defaults to 201 equally spaced points on the domain.
    deriv : tuple of int, optional
        Derivative orders on the horizontal and vertical axes.  Default
        ``(1, 2)``.
    labels : mapping of float to str, optional
        Text written on every curve at its position at the given times, e.g.
        month names on a year of data.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  A new figure is created when omitted.
    **kwargs
        Passed to :meth:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes drawn on, with one line per curve.

    Raises
    ------
    ValueError
        If ``fd`` has more than one variable.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.dynamics import phase_plane
    >>> basis = fdt.Fourier(domain=(0.0, 1.0), n_basis=5)
    >>> fd = fdt.FData(np.array([[0.0, 1.0, 0.0, 0.2, 0.0]]).T, basis)
    >>> len(phase_plane(fd).lines)
    1
    """
    if fd.n_vars != 1:
        raise ValueError(f"phase_plane plots single-variable curves, got {fd.n_vars} variables")
    xp = default_namespace()
    points = (
        xp.linspace(fd.domain[0], fd.domain[1], _PLOT_POINTS, dtype=xp.float64)
        if t is None
        else asarray(to_numpy(t))
    )
    horizontal = to_numpy(fd(points, deriv[0]))
    vertical = to_numpy(fd(points, deriv[1]))
    axes = _axes(ax)
    for n in range(fd.n_curves):
        axes.plot(horizontal[:, n], vertical[:, n], **kwargs)
        if labels:

            def at(order: int, curve: int = n) -> Any:
                return lambda when: to_numpy(fd(xp.asarray([when], dtype=xp.float64), order))[
                    0, curve
                ]

            _annotate(axes, labels, at(deriv[0]), at(deriv[1]))
    axes.set_xlabel(f"D{deriv[0]} x")
    axes.set_ylabel(f"D{deriv[1]} x")
    return axes
