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

Replaces R's ``pda.fd`` (:class:`PDA`), ``pda.overlay``
(:meth:`PDA.plot_overlay`) and ``phaseplanePlot`` (:func:`phase_plane`); the
fitted equation is integrated with :func:`scipy.integrate.solve_ivp`
(:meth:`PDA.solve`).

Examples
--------
>>> import numpy as np
>>> import fabel as fb
>>> from fabel.dynamics import PDA
>>> basis = fb.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
>>> curves = fb.FData(np.array([[0.0, 0.0], [1.0, 0.5], [0.0, 2.0]]), basis)
>>> pda = PDA(order=2, n_grid=None).fit(curves)  # a sin t + b cos t solve D²x + x = 0
>>> [round(float(w.coefs[0, 0]), 10) + 0.0 for w in pda.weights_]
[1.0, 0.0]
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from types import ModuleType
from typing import TYPE_CHECKING, Any

from scipy.integrate import solve_ivp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from fabel import _linalg
from fabel._backend import asarray, default_namespace, to_numpy
from fabel._operator import LDO
from fabel.basis import Basis, Constant, _same_domain
from fabel.core import FData, _quadrature

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes

__all__ = ["PDA", "phase_plane"]

Array = Any

#: Points on the equally spaced grid R's ``pda.fd`` integrates with the
#: trapezoidal rule; the default, so that Fabel reproduces ``pda.fd``.
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


def _fit_equation(
    target: Array,
    derivs: list[Array],
    blocks: list[Array],
    penalties: list[Array | None],
    lams: list[float],
    quad: Array,
    xp: ModuleType,
) -> list[list[Array]]:
    """Estimate the weight coefficients of one equation.

    ``target`` is ``D^m x_i`` of shape ``(n_nodes, N)``; ``derivs[j]`` is
    ``D^j x`` of shape ``(n_nodes, N, d)``; ``blocks[j]`` the weight basis
    ``j`` at the nodes.  Returns ``coefs[k][j]`` for variable ``k`` and
    derivative ``j``.
    """
    order = len(blocks)
    n_curves = target.shape[1]
    n_vars = derivs[0].shape[2]
    slots = [(k, j) for k in range(n_vars) for j in range(order)]
    sizes = [blocks[j].shape[1] for _, j in slots]
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1] + size)
    matrix = xp.zeros((offsets[-1], offsets[-1]), dtype=xp.float64)
    rhs = xp.zeros((offsets[-1],), dtype=xp.float64)
    for p, (k, j) in enumerate(slots):
        left = derivs[j][:, :, k]
        for q, (kk, jj) in enumerate(slots[: p + 1]):
            right = derivs[jj][:, :, kk]
            kernel = xp.sum(left * right, axis=1) * quad / n_curves
            block = xp.matmul(xp.matrix_transpose(blocks[j]), kernel[:, None] * blocks[jj])
            matrix[offsets[p] : offsets[p + 1], offsets[q] : offsets[q + 1]] = block
            matrix[offsets[q] : offsets[q + 1], offsets[p] : offsets[p + 1]] = xp.matrix_transpose(
                block
            )
        forcing = xp.sum(left * target, axis=1) * quad / n_curves
        rhs[offsets[p] : offsets[p + 1]] = -xp.matmul(xp.matrix_transpose(blocks[j]), forcing)
        penalty = penalties[j]
        if penalty is not None:
            matrix[offsets[p] : offsets[p + 1], offsets[p] : offsets[p + 1]] += lams[j] * penalty
    solution = _solve_normal_equations(matrix, rhs, xp)
    coefs: list[list[Array]] = [[None] * order for _ in range(n_vars)]
    for p, (k, j) in enumerate(slots):
        coefs[k][j] = solution[offsets[p] : offsets[p + 1]]
    return coefs


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
# the estimator
# --------------------------------------------------------------------------- #


class PDA(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    r"""Principal differential analysis: fit a linear ODE to a sample of curves.

    Replaces R's ``pda.fd``.  Estimates weight functions ``β_j`` so that
    ``D^m x + Σ_j β_j D^j x`` is as small as possible over the sample, in the
    mean integrated square, with an optional roughness penalty on each weight
    (see the module docstring for the criterion).  Curves with several
    variables (``n_vars = d > 1``) define a system of ``d`` coupled equations.

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
        Roughness operator applied to the weights.  Default ``2``.
    n_grid : int or None, optional
        How the integrals are computed.  An integer is the number of equally
        spaced points of a trapezoidal rule; the default ``501`` is the rule
        R's ``pda.fd`` uses, so Fabel reproduces it.  The residual functions
        are then least-squares fits of ``Lx`` on the same grid.  ``None``
        integrates exactly (composite Gauss-Legendre on the break points,
        exact for spline curves and weights) and projects the residuals in
        ``L²``; it differs from R's rule by about ``1e-5`` relative on
        realistic data.

    Attributes
    ----------
    weights_ : tuple
        For one variable, ``weights_[j]`` is ``β_j`` as a one-curve
        :class:`~fabel.core.FData`.  For ``d`` variables, ``weights_[i][k][j]``
        multiplies ``D^j x_k`` in equation ``i``.  This mirrors the nesting of
        R's ``bwtlist``.
    residuals_ : FData
        ``L x`` for the fitted curves, in their basis.
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
    >>> import fabel as fb
    >>> from fabel.dynamics import PDA
    >>> t = np.linspace(0.0, 1.0, 101)
    >>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=24, order=5)
    >>> from fabel.smoothing import smooth
    >>> fd = smooth(np.exp(-4 * t), t, basis=basis, lam=0.0).fd
    >>> pda = PDA(order=1).fit(fd)  # Dx + βx = 0
    >>> round(float(pda.weights_[0].coefs[0, 0]), 6)
    4.0
    """

    def __init__(
        self,
        order: int = 2,
        *,
        weight_basis: Basis | Sequence[Basis] | None = None,
        lam: float | Sequence[float] = 0.0,
        penalty: int | LDO = 2,
        n_grid: int | None = _R_GRID,
    ) -> None:
        self.order = order
        self.weight_basis = weight_basis
        self.lam = lam
        self.penalty = penalty
        self.n_grid = n_grid

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

    def fit(self, X: FData, y: Any = None) -> PDA:  # noqa: N803
        """Estimate the weight functions from a sample of curves.

        Parameters
        ----------
        X : FData
            The curves; a multivariate :class:`~fabel.core.FData` (``n_vars >
            1``) is fitted as a system of coupled equations.
        y : None
            Ignored; present for scikit-learn compatibility.

        Returns
        -------
        PDA
            The fitted estimator.

        Raises
        ------
        TypeError
            If ``X`` is not an :class:`~fabel.core.FData`.
        ValueError
            If a parameter is invalid or the normal equations are singular.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> from fabel.dynamics import PDA
        >>> basis = fb.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fb.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> len(PDA(order=2).fit(fd).weights_)
        2
        """
        if not isinstance(X, FData):
            raise TypeError(f"PDA.fit expects an FData, got {type(X).__name__}")
        xp = default_namespace()
        order = self._checked_order()
        bases = self._weight_bases(X.domain, order)
        lams = self._lambdas(order)
        penalties: list[Array | None] = [
            asarray(basis.penalty(self.penalty)) if lam > 0.0 else None
            for basis, lam in zip(bases, lams, strict=True)
        ]
        nodes, quad = self._rule(X, bases, xp)
        derivs = _derivatives(X, nodes, order)
        blocks = [asarray(basis(nodes)) for basis in bases]
        n_vars = derivs[0].shape[2]
        weights: list[list[list[FData]]] = []
        for i in range(n_vars):
            coefs = _fit_equation(
                derivs[order][:, :, i], derivs[:order], blocks, penalties, lams, quad, xp
            )
            weights.append(
                [[FData(coefs[k][j], bases[j]) for j in range(order)] for k in range(n_vars)]
            )
        self.n_vars_ = n_vars
        self.domain_ = X.domain
        if n_vars == 1:
            self.weights_: tuple[Any, ...] = tuple(weights[0][0])
            self.operator_: LDO | None = LDO(weights=list(self.weights_))
        else:
            self.weights_ = tuple(tuple(tuple(row) for row in eq) for eq in weights)
            self.operator_ = None
        self._nested = weights
        self.residuals_ = self._residuals(X, nodes, quad, derivs)
        return self

    def _residuals(self, curves: FData, nodes: Array, quad: Array, derivs: list[Array]) -> FData:
        """Return ``L x`` for ``curves`` expressed in their own basis."""
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

    def transform(self, X: FData) -> FData:  # noqa: N803
        """Apply the fitted operator to curves: return ``L x``.

        Parameters
        ----------
        X : FData
            Curves with the same number of variables and domain as the fitted
            ones.

        Returns
        -------
        FData
            The residual functions ``L x``, in the basis of ``X``.

        Raises
        ------
        TypeError
            If ``X`` is not an :class:`~fabel.core.FData`.
        ValueError
            If ``X`` does not match the fitted curves.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> from fabel.dynamics import PDA
        >>> basis = fb.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fb.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> pda = PDA(order=2, n_grid=None).fit(fd)
        >>> new = fb.FData(np.array([0.0, 2.0, -1.0]), basis)
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
        xp = default_namespace()
        order = self._checked_order()
        bases = [w.basis for w in self._nested[0][0]]
        nodes, quad = self._rule(X, bases, xp)
        return self._residuals(X, nodes, quad, _derivatives(X, nodes, order))

    # --------------------------------------------------------------- solving

    def solve(
        self,
        t: Any,
        initial: Any,
        *,
        rtol: float = _ODE_RTOL,
        atol: float = _ODE_ATOL,
    ) -> Array:
        """Integrate the fitted homogeneous equation ``L x = 0``.

        The equation is rewritten as a first-order system in ``(x, Dx, ...,
        D^{m-1} x)`` and integrated with :func:`scipy.integrate.solve_ivp`
        (method ``DOP853``) from ``t[0]``, forwards or backwards.

        Parameters
        ----------
        t : array_like
            Strictly monotone points inside the domain; the initial condition
            holds at ``t[0]``.
        initial : array_like
            ``(x, Dx, ..., D^{m-1} x)`` at ``t[0]``: shape ``(order,)`` for one
            variable, ``(n_vars, order)`` for a system.
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
            ``initial`` has the wrong shape, or if the integrator fails.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> from fabel.dynamics import PDA
        >>> basis = fb.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fb.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
        >>> pda = PDA(order=2, n_grid=None).fit(fd)
        >>> t = np.linspace(0.0, np.pi, 5)
        >>> bool(np.allclose(pda.solve(t, [0.0, 1.0]), np.sin(t), atol=1e-8))
        True
        """
        check_is_fitted(self, "weights_")
        xp = default_namespace()
        order = self._checked_order()
        n_vars = self.n_vars_
        points = asarray(to_numpy(t))
        if len(points.shape) != 1 or points.shape[0] == 0:
            raise ValueError("t must be a one-dimensional array with at least one point")
        steps = points[1:] - points[:-1]
        if points.shape[0] > 1 and not (bool(xp.all(steps > 0.0)) or bool(xp.all(steps < 0.0))):
            raise ValueError("t must be strictly monotone")
        lower, upper = self.domain_
        slack = _DOMAIN_TOL * max(1.0, upper - lower)
        if float(xp.min(points)) < lower - slack or float(xp.max(points)) > upper + slack:
            raise ValueError(f"t must lie in the domain {self.domain_}")
        start = asarray(to_numpy(initial))
        expected = (order,) if n_vars == 1 else (n_vars, order)
        if tuple(start.shape) != expected:
            raise ValueError(f"initial must have shape {expected}, got {tuple(start.shape)}")
        state0 = xp.reshape(start, (n_vars * order,))
        if points.shape[0] == 1:
            first = xp.reshape(state0, (n_vars, order))[:, 0]
            return first if n_vars == 1 else first[None, :]
        clip_lo, clip_hi = lower, upper
        nested = self._nested

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
                dz[i, -1] = -total
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
        >>> import fabel as fb
        >>> from fabel.dynamics import PDA
        >>> basis = fb.Fourier(domain=(0.0, 2 * np.pi), n_basis=3)
        >>> fd = fb.FData(np.array([[0.0, 0.0], [1.0, 0.3], [0.2, 1.0]]), basis)
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
    >>> import fabel as fb
    >>> from fabel.dynamics import phase_plane
    >>> basis = fb.Fourier(domain=(0.0, 1.0), n_basis=5)
    >>> fd = fb.FData(np.array([[0.0, 1.0, 0.0, 0.2, 0.0]]).T, basis)
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
