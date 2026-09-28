r"""Penalised density and intensity estimation.

Two estimators replace R ``fda``'s ``density.fd`` (dropped from fda 6.3.0) and
``intensity.fd``.  Both write a positive function as the exponential of a basis
expansion ``W(x) = Σ_k c_k φ_k(x)`` and maximise a penalised log-likelihood.

:func:`fit_density` estimates a probability density
``p(x) = C exp W(x)`` from a sample ``x_1 .. x_n`` by minimising

.. math::

    F(c) = -\sum_{i=1}^{n} W(x_i) + n \log \int e^{W(u)}\,du
           + \lambda \int (LW)^2 ,

with ``C = 1 / ∫ e^W`` the normalising constant (Ramsay, Hooker & Graves
2009, § 5.4; Silverman 1986 for penalised log-density estimation).

:func:`fit_intensity` estimates the intensity ``μ(t) = exp W(t)`` of a
non-homogeneous Poisson process from its event times by minimising the
penalised negative log-likelihood of the process,

.. math::

    F(c) = -\sum_{i=1}^{n} W(x_i) + \int e^{W(u)}\,du + \lambda \int (LW)^2 .

The penalty is ``λ cᵀRc`` with ``R`` the roughness matrix of the linear
differential operator ``L`` (:meth:`fdatools.Basis.penalty`), the same scaling as
R.  Both criteria are convex in ``c``.  They are minimised by Newton's method
with the exact Hessian and a backtracking (Armijo) line search; the integrals
are computed by composite Gauss-Legendre quadrature, 12 nodes on each of the
basis's refined quadrature panels, which is exact to rounding for the
exponential of a polynomial piece of the size met in practice.

Convergence: the iteration stops once a full Newton step, which is then
taken, moves no coefficient by more than ``tol * (1 + max|c|)`` (default
``tol = 1e-10``); Newton's quadratic convergence leaves the result about
``tol²`` from the optimum.  If rounding hides every further decrease of ``F``
first, the fit counts as converged when the last proposed step was within
``sqrt(tol) * (1 + max|c|)``.

When the basis reproduces the constant function and the penalty does not
charge it (every derivative operator ``Dᵐ`` with ``m >= 1``), the density
criterion does not change under ``W -> W + s``.  The fit then returns the
member of that family with ``∫ W = 0`` over the domain; the density itself,
``C exp W``, is unique.

Examples
--------
>>> import numpy as np
>>> import fdatools as fdt
>>> from fdatools.density import fit_density
>>> rng = np.random.default_rng(0)
>>> x = rng.normal(size=300)
>>> basis = fdt.BSpline(domain=(-5.0, 5.0), n_basis=11)
>>> result = fit_density(x, basis=basis, lam=1e-2)
>>> grid = np.linspace(-5.0, 5.0, 2001)
>>> mass = float(np.sum(result(grid)) * (grid[1] - grid[0]))
>>> round(mass, 3)
1.0
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from math import exp, isfinite, log
from types import ModuleType
from typing import Any

from fdatools import _linalg
from fdatools._backend import array_namespace, asarray, default_namespace, to_numpy
from fdatools._operator import LDO
from fdatools.basis import Basis, BSpline
from fdatools.core import FData

__all__ = ["DensityResult", "IntensityResult", "fit_density", "fit_intensity"]

Array = Any

#: Gauss-Legendre nodes per quadrature panel for ``∫ exp W``.
_QUAD_DEGREE = 12

#: Default relative step tolerance of the Newton iteration.
DEFAULT_TOL = 1e-10

#: Default iteration cap of the Newton iteration.
DEFAULT_MAX_ITER = 100

#: Line-search parameters: sufficient-decrease constant and halving budget.
_ARMIJO = 1e-4
_MAX_HALVINGS = 60

#: ``log`` of the largest finite double, rounded down: ``exp(W)`` above it
#: overflows, so a trial point reaching it is rejected by the line search.
_LOG_OVERFLOW = 700.0

#: Relative thresholds for "the basis reproduces constants" and "the penalty
#: does not charge them".
_CONSTANT_TOL = 1e-12
_NULL_PENALTY_TOL = 1e-10

#: Default basis: cubic B-splines, ``min(max(4, n // 10), 20)`` functions.
_AUTO_ORDER = 4
_AUTO_MIN = 4
_AUTO_MAX = 20
_AUTO_PER = 10


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


def _single_curve(values: Array) -> Array:
    """Drop the curve axis of a one-curve evaluation ``(n_t, 1) -> (n_t,)``."""
    return values[:, 0]


@dataclass(frozen=True)
class DensityResult:
    r"""Outcome of :func:`fit_density`.

    The estimated density is ``p(x) = normaliser * exp(W(x))``.

    Attributes
    ----------
    fd : FData
        The log-density function ``W`` up to the additive constant
        ``log(normaliser)``.  When ``W`` is only defined up to a constant (see
        the module notes) it is centred so that ``∫ W = 0``.
    normaliser : float
        ``C = 1 / ∫ exp W`` over the basis domain.
    lam : float
        The smoothing parameter.
    penalty_matrix : array
        The roughness matrix ``R``; the penalty is ``lam * cᵀ R c``.
    criterion : float
        The minimised penalised negative log-likelihood ``F(c)``.
    log_likelihood : float
        ``Σ log p(x_i)`` at the estimate.
    gradient_norm : float
        Euclidean norm of the gradient of ``F`` at the estimate.
    n_iter : int
        Number of Newton iterations taken.
    converged : bool
        Whether the step tolerance was met within ``max_iter`` iterations.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.density import fit_density
    >>> x = np.array([0.1, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.8])
    >>> res = fit_density(x, basis=fdt.BSpline(domain=(0.0, 1.0), n_basis=6), lam=1e-3)
    >>> res.converged
    True
    >>> bool(np.all(res(np.array([0.0, 0.5, 1.0])) > 0))
    True
    """

    fd: FData
    normaliser: float
    lam: float
    penalty_matrix: Array
    criterion: float
    log_likelihood: float
    gradient_norm: float
    n_iter: int
    converged: bool

    def __call__(self, t: Any) -> Array:
        """Evaluate the density ``p(t) = C exp W(t)``.

        Parameters
        ----------
        t : array_like
            One-dimensional evaluation points inside the basis domain.  A
            :class:`torch.Tensor` gives a tensor back.

        Returns
        -------
        array
            Density values, shape ``(n_t,)``.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.density import fit_density
        >>> x = np.array([0.2, 0.3, 0.5, 0.6, 0.65, 0.9])
        >>> res = fit_density(x, basis=fdt.BSpline(domain=(0.0, 1.0), n_basis=5), lam=1.0)
        >>> res(np.array([0.5])).shape
        (1,)
        """
        values = self.log_density(t)
        return array_namespace(values).exp(values)

    def log_density(self, t: Any) -> Array:
        """Evaluate the log-density ``log p(t) = W(t) + log C``.

        Parameters
        ----------
        t : array_like
            One-dimensional evaluation points inside the basis domain.

        Returns
        -------
        array
            Log-density values, shape ``(n_t,)``.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.density import fit_density
        >>> x = np.array([0.2, 0.3, 0.5, 0.6, 0.65, 0.9])
        >>> res = fit_density(x, basis=fdt.BSpline(domain=(0.0, 1.0), n_basis=5), lam=1.0)
        >>> t = np.array([0.25, 0.75])
        >>> bool(np.allclose(np.exp(res.log_density(t)), res(t)))
        True
        """
        return _single_curve(self.fd(t)) + log(self.normaliser)


@dataclass(frozen=True)
class IntensityResult:
    r"""Outcome of :func:`fit_intensity`.

    The estimated intensity is ``μ(t) = exp(W(t))``.

    Attributes
    ----------
    fd : FData
        The log-intensity function ``W``.
    lam : float
        The smoothing parameter.
    penalty_matrix : array
        The roughness matrix ``R``; the penalty is ``lam * cᵀ R c``.
    criterion : float
        The minimised penalised negative log-likelihood ``F(c)`` (R's
        ``Flist$f``).
    log_likelihood : float
        The Poisson-process log-likelihood ``Σ W(x_i) - ∫ exp W``.
    expected_count : float
        ``∫ exp W`` over the domain, the expected number of events.
    gradient_norm : float
        Euclidean norm of the gradient of ``F`` at the estimate.
    n_iter : int
        Number of Newton iterations taken.
    converged : bool
        Whether the step tolerance was met within ``max_iter`` iterations.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.density import fit_intensity
    >>> events = np.cumsum(np.full(40, 0.25))
    >>> res = fit_intensity(events, basis=fdt.BSpline(domain=(0.0, 10.0), n_basis=5), lam=10.0)
    >>> round(res.expected_count, 6)
    40.0
    """

    fd: FData
    lam: float
    penalty_matrix: Array
    criterion: float
    log_likelihood: float
    expected_count: float
    gradient_norm: float
    n_iter: int
    converged: bool

    def __call__(self, t: Any) -> Array:
        """Evaluate the intensity ``μ(t) = exp W(t)``.

        Parameters
        ----------
        t : array_like
            One-dimensional evaluation points inside the basis domain.  A
            :class:`torch.Tensor` gives a tensor back.

        Returns
        -------
        array
            Intensity values (events per unit time), shape ``(n_t,)``.

        Examples
        --------
        >>> import numpy as np
        >>> import fdatools as fdt
        >>> from fdatools.density import fit_intensity
        >>> events = np.linspace(0.5, 9.5, 19)
        >>> basis = fdt.BSpline(domain=(0.0, 10.0), n_basis=4)
        >>> res = fit_intensity(events, basis=basis, penalty=1, lam=1e6)
        >>> np.round(res(np.array([2.0, 8.0])), 3).tolist()
        [1.9, 1.9]
        """
        values = _single_curve(self.fd(t))
        return array_namespace(values).exp(values)


# --------------------------------------------------------------------------- #
# the optimisation problem
# --------------------------------------------------------------------------- #


class _Problem:
    """Penalised log-likelihood of a density or intensity, with derivatives.

    ``kind`` is ``"density"`` or ``"intensity"``.  All arrays are in the
    default (NumPy) namespace.
    """

    def __init__(self, kind: str, basis: Basis, x: Array, pen: Array, lam: float) -> None:
        xp = default_namespace()
        self.xp = xp
        self.kind = kind
        self.n = int(x.shape[0])
        self.data_sum = xp.sum(basis(x), axis=0)
        nodes, weights = _linalg.composite_gauss_legendre(
            to_numpy(basis._quadrature_panels()), _QUAD_DEGREE
        )
        self.design = asarray(basis(asarray(nodes, xp)), xp)
        self.weights = asarray(weights, xp)
        self.pen = pen
        self.lam = lam

    def latent(self, coefs: Array) -> Array:
        """``W`` at the quadrature nodes."""
        return self.xp.matmul(self.design, coefs)

    def value(self, coefs: Array) -> float:
        """Return the criterion ``F(c)``; ``inf`` where ``exp W`` would overflow."""
        xp = self.xp
        w = self.latent(coefs)
        top = float(xp.max(w))
        fit = -float(xp.sum(self.data_sum * coefs))
        rough = self.lam * float(xp.sum(coefs * xp.matmul(self.pen, coefs)))
        if self.kind == "density":
            total = float(xp.sum(self.weights * xp.exp(w - top)))
            return fit + self.n * (top + log(total)) + rough
        if top > _LOG_OVERFLOW:
            return float("inf")
        return fit + float(xp.sum(self.weights * xp.exp(w))) + rough

    def derivatives(self, coefs: Array) -> tuple[Array, Array]:
        """Gradient and Hessian of ``F`` at ``coefs``."""
        xp = self.xp
        w = self.latent(coefs)
        if self.kind == "density":
            mass = self.weights * xp.exp(w - float(xp.max(w)))
            mass = mass / xp.sum(mass)
            scale = float(self.n)
        else:
            mass = self.weights * xp.exp(w)
            scale = 1.0
        first = xp.matmul(xp.matrix_transpose(self.design), mass)
        second = xp.matmul(xp.matrix_transpose(self.design), mass[:, None] * self.design)
        if self.kind == "density":
            second = second - first[:, None] * first[None, :]
        rough = 2.0 * self.lam
        gradient = -self.data_sum + scale * first + rough * xp.matmul(self.pen, coefs)
        hessian = scale * second + rough * self.pen
        return gradient, hessian

    def integral(self, coefs: Array) -> float:
        """``∫ exp W`` over the domain."""
        xp = self.xp
        return float(xp.sum(self.weights * xp.exp(self.latent(coefs))))

    def log_integral(self, coefs: Array) -> float:
        """``log ∫ exp W`` over the domain, without overflow."""
        xp = self.xp
        w = self.latent(coefs)
        top = float(xp.max(w))
        return top + log(float(xp.sum(self.weights * xp.exp(w - top))))

    def mean_latent(self, coefs: Array) -> float:
        """``∫ W`` over the domain."""
        return float(self.xp.sum(self.weights * self.latent(coefs)))

    def constant_coefs(self) -> Array | None:
        """Coefficients of the constant ``1`` if the basis reproduces it, else ``None``."""
        xp = self.xp
        gram = xp.matmul(xp.matrix_transpose(self.design), self.weights[:, None] * self.design)
        rhs = xp.matmul(xp.matrix_transpose(self.design), self.weights)
        coefs = _linalg.lstsq(gram, rhs)
        resid = 1.0 - xp.matmul(self.design, coefs)
        length = float(xp.sum(self.weights))
        if float(xp.sum(self.weights * resid * resid)) > _CONSTANT_TOL * length:
            return None
        return coefs


def _null_direction(problem: _Problem) -> Array | None:
    """Return the unidentified direction of a density criterion, if there is one."""
    if problem.kind != "density":
        return None
    unit = problem.constant_coefs()
    if unit is None:
        return None
    xp = problem.xp
    if problem.lam > 0.0:
        charge = float(xp.sum(unit * xp.matmul(problem.pen, unit)))
        size = float(xp.max(xp.abs(problem.pen))) * float(xp.sum(unit * unit))
        if charge > _NULL_PENALTY_TOL * max(size, 1.0):
            return None
    return unit


def _newton(
    problem: _Problem, start: Array, null: Array | None, tol: float, max_iter: int
) -> tuple[Array, int, bool]:
    """Minimise ``problem`` by damped Newton from ``start``.

    Each step is halved until the Armijo sufficient-decrease condition holds.
    Close to the optimum the decrease of ``F`` along weakly determined
    directions (tail coefficients carrying almost no mass) drops below the
    rounding error of ``F`` itself; the full step is then accepted when it
    reduces the gradient norm instead.

    With a null direction ``u`` the Hessian is singular along ``u`` and the
    gradient is orthogonal to it; solving with ``H + s u uᵀ/|u|²`` then gives
    the Newton step orthogonal to ``u``.
    """
    xp = problem.xp
    coefs = start
    value = problem.value(coefs)
    if not isfinite(value):
        raise ValueError("the starting coefficients overflow exp(W); pass a smaller start")
    gradient, hessian = problem.derivatives(coefs)
    for iteration in range(1, max_iter + 1):
        if null is not None:
            size = float(xp.sum(xp.abs(xp.linalg.diagonal(hessian)))) / hessian.shape[0]
            direction = null / float(xp.sqrt(xp.sum(null * null)))
            hessian = hessian + max(size, 1.0) * direction[:, None] * direction[None, :]
        step = -_linalg.solve_spd(hessian, gradient)
        small = _max_abs(step, xp) <= tol * (1.0 + _max_abs(coefs, xp))
        slope = float(xp.sum(gradient * step))
        scale = 1.0
        accepted = False
        for _ in range(_MAX_HALVINGS):
            candidate = coefs + scale * step
            trial = problem.value(candidate)
            if trial < value and trial <= value + _ARMIJO * scale * slope:
                accepted = True
                break
            scale *= 0.5
        if accepted:
            new_gradient, new_hessian = problem.derivatives(candidate)
        else:
            scale = 1.0
            candidate = coefs + step
            trial = problem.value(candidate)
            progress = isfinite(trial)
            if progress:
                new_gradient, new_hessian = problem.derivatives(candidate)
                progress = _norm(new_gradient, xp) < _norm(gradient, xp)
            if not progress:
                # Neither test can see progress: the iterate is optimal to
                # rounding.  Converged when the proposed step was small.
                moved = _max_abs(step, xp)
                return coefs, iteration, moved <= _relaxed(tol) * (1.0 + _max_abs(coefs, xp))
        coefs, value, gradient, hessian = candidate, trial, new_gradient, new_hessian
        if small and scale == 1.0:
            return coefs, iteration, True
    return coefs, max_iter, False


def _norm(values: Array, xp: ModuleType) -> float:
    return float(xp.sqrt(xp.sum(values * values)))


def _relaxed(tol: float) -> float:
    """Relaxed tolerance ``sqrt(tol)`` used when the line search stalls."""
    return float(tol**0.5)


def _max_abs(values: Array, xp: ModuleType) -> float:
    return float(xp.max(xp.abs(values)))


# --------------------------------------------------------------------------- #
# argument handling
# --------------------------------------------------------------------------- #


def _sample(x: Any, name: str) -> Array:
    """Validate a one-dimensional finite sample and return it as float64."""
    xp = default_namespace()
    values = asarray(to_numpy(x), xp)
    if values.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {tuple(values.shape)}")
    if values.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one value")
    if not bool(xp.all(xp.isfinite(values))):
        raise ValueError(f"{name} must be finite")
    return values


def _default_basis(domain: tuple[float, float], n: int) -> BSpline:
    n_basis = min(max(_AUTO_MIN, n // _AUTO_PER), _AUTO_MAX)
    return BSpline(domain=domain, n_basis=n_basis, order=_AUTO_ORDER)


def _resolve_domain(
    domain: tuple[float, float] | None, default: tuple[float, float]
) -> tuple[float, float]:
    lower, upper = (float(v) for v in (default if domain is None else domain))
    if not upper > lower:
        raise ValueError(f"domain must have lower < upper, got ({lower}, {upper})")
    return lower, upper


def _check_inside(x: Array, basis: Basis, name: str) -> None:
    xp = default_namespace()
    lower, upper = basis.domain
    low, high = float(xp.min(x)), float(xp.max(x))
    if low < lower or high > upper:
        raise ValueError(f"{name} must lie in the basis domain {basis.domain}; got [{low}, {high}]")


def _check_settings(lam: float, tol: float, max_iter: int) -> float:
    value = float(lam)
    if not isfinite(value) or value < 0.0:
        raise ValueError(f"lam must be a finite non-negative number, got {lam}")
    if not tol > 0.0:
        raise ValueError(f"tol must be positive, got {tol}")
    if int(max_iter) < 1:
        raise ValueError(f"max_iter must be at least 1, got {max_iter}")
    return value


def _start(start: Any, basis: Basis) -> Array:
    xp = default_namespace()
    if start is None:
        return xp.zeros(basis.n_basis, dtype=xp.float64)
    coefs = start.coefs if isinstance(start, FData) else start
    values = xp.reshape(asarray(to_numpy(coefs), xp), (-1,))
    if values.shape[0] != basis.n_basis:
        raise ValueError(f"start must have {basis.n_basis} coefficients, got {values.shape[0]}")
    return values


def _fit(
    kind: str,
    values: Array,
    basis: Basis,
    lam: float,
    penalty: int | LDO,
    start: Any,
    tol: float,
    max_iter: int,
) -> tuple[_Problem, Array, Array, int, bool, float]:
    """Shared driver: build the problem, run Newton, report convergence."""
    xp = default_namespace()
    operator = penalty if isinstance(penalty, LDO) else LDO(int(penalty))
    pen = asarray(basis.penalty(operator), xp)
    problem = _Problem(kind, basis, values, pen, lam)
    null = _null_direction(problem)
    coefs, n_iter, converged = _newton(problem, _start(start, basis), null, tol, int(max_iter))
    if not converged:
        warnings.warn(
            f"fit_{kind} did not converge in {n_iter} Newton iterations "
            f"(tol={tol}); the result is the last iterate",
            RuntimeWarning,
            stacklevel=3,
        )
    if null is not None:
        length = float(xp.sum(problem.weights))
        coefs = coefs - (problem.mean_latent(coefs) / length) * null
    gradient, _ = problem.derivatives(coefs)
    grad_norm = _norm(gradient, xp)
    return problem, coefs, pen, n_iter, converged, grad_norm


# --------------------------------------------------------------------------- #
# public entry points
# --------------------------------------------------------------------------- #


def fit_density(
    x: Any,
    *,
    basis: Basis | None = None,
    domain: tuple[float, float] | None = None,
    lam: float = 0.0,
    penalty: int | LDO = 2,
    start: Any = None,
    tol: float = DEFAULT_TOL,
    max_iter: int = DEFAULT_MAX_ITER,
) -> DensityResult:
    r"""Estimate a probability density by penalised maximum likelihood.

    Replaces R ``fda``'s ``density.fd``.  The density is
    ``p(x) = C exp W(x)`` with ``W`` a basis expansion; ``W`` minimises
    ``-Σ W(x_i) + n log ∫ exp W + lam * ∫ (L W)²``.

    Parameters
    ----------
    x : array_like
        One-dimensional sample (a :class:`torch.Tensor` is accepted and read as
        data; the fit itself is not differentiable).
    basis : Basis, optional
        Basis for ``W``; its domain is the support of the density.  Default:
        cubic B-splines with ``min(max(4, n // 10), 20)`` functions on
        ``domain``.
    domain : tuple of float, optional
        Support used for the default basis; default ``(min x, max x)``.
        Ignored when ``basis`` is given.
    lam : float, optional
        Smoothing parameter ``λ >= 0`` (default ``0``, no penalty).
    penalty : int or LDO, optional
        Derivative order or linear differential operator ``L`` of the roughness
        penalty (default ``2``).  A heavy first-derivative penalty pulls the
        density towards the uniform one; a heavy ``D³`` penalty towards a
        (truncated) normal.
    start : FData or array_like, optional
        Starting coefficients of ``W``; default zero (the uniform density).
    tol : float, optional
        Newton step tolerance: stop once a full Newton step moves no
        coefficient by more than ``tol * (1 + max|c|)``.
    max_iter : int, optional
        Iteration cap.

    Returns
    -------
    DensityResult
        ``result(t)`` evaluates the density, ``result.fd`` is ``W`` and
        ``result.normaliser`` is ``C``.

    Raises
    ------
    ValueError
        If ``x`` is empty, not finite or not one-dimensional, lies outside the
        basis domain, or a setting is out of range.

    Warns
    -----
    RuntimeWarning
        If the Newton iteration has not converged after ``max_iter`` steps.

    Notes
    -----
    If the basis reproduces constants and ``L`` annihilates them, ``W`` is
    fixed only up to an additive constant; the returned ``W`` has ``∫ W = 0``.
    In that case the fit coincides with :func:`fit_intensity` on the same
    data, normalised: ``p = μ / n``.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.density import fit_density
    >>> rng = np.random.default_rng(1)
    >>> x = rng.exponential(size=500)
    >>> res = fit_density(x, domain=(0.0, 10.0), lam=1e-2)
    >>> bool(res(np.array([0.1]))[0] > res(np.array([3.0]))[0])
    True
    """
    values = _sample(x, "x")
    lam_value = _check_settings(lam, tol, max_iter)
    if basis is None:
        xp = default_namespace()
        support = _resolve_domain(domain, (float(xp.min(values)), float(xp.max(values))))
        basis = _default_basis(support, int(values.shape[0]))
    _check_inside(values, basis, "x")
    problem, coefs, pen, n_iter, converged, grad_norm = _fit(
        "density", values, basis, lam_value, penalty, start, tol, max_iter
    )
    xp = problem.xp
    log_total = problem.log_integral(coefs)
    normaliser = exp(-log_total)
    criterion = problem.value(coefs)
    fitted = float(xp.sum(problem.data_sum * coefs))
    return DensityResult(
        fd=FData(coefs, basis),
        normaliser=normaliser,
        lam=lam_value,
        penalty_matrix=pen,
        criterion=criterion,
        log_likelihood=fitted - problem.n * log_total,
        gradient_norm=grad_norm,
        n_iter=n_iter,
        converged=converged,
    )


def fit_intensity(
    x: Any,
    *,
    basis: Basis | None = None,
    domain: tuple[float, float] | None = None,
    lam: float = 0.0,
    penalty: int | LDO = 2,
    start: Any = None,
    tol: float = DEFAULT_TOL,
    max_iter: int = DEFAULT_MAX_ITER,
) -> IntensityResult:
    r"""Estimate the intensity of a Poisson process by penalised likelihood.

    Replaces R ``fda``'s ``intensity.fd``.  The intensity is
    ``μ(t) = exp W(t)``; ``W`` minimises
    ``-Σ W(x_i) + ∫ exp W + lam * ∫ (L W)²`` over the basis domain.

    Parameters
    ----------
    x : array_like
        One-dimensional event times (a :class:`torch.Tensor` is accepted and
        read as data).  Their order does not matter.
    basis : Basis, optional
        Basis for ``W``; its domain is the observation window.  Default: cubic
        B-splines with ``min(max(4, n // 10), 20)`` functions on ``domain``.
    domain : tuple of float, optional
        Observation window for the default basis; default ``(0, max x)``, R's
        convention that observation starts at time 0.  Ignored when ``basis``
        is given.
    lam : float, optional
        Smoothing parameter ``λ >= 0`` (default ``0``).
    penalty : int or LDO, optional
        Derivative order or operator ``L`` of the roughness penalty (default
        ``2``).  A heavy first-derivative penalty pulls the process towards a
        homogeneous one.
    start : FData or array_like, optional
        Starting coefficients of ``W``; default zero.
    tol : float, optional
        Newton step tolerance: stop once a full Newton step moves no
        coefficient by more than ``tol * (1 + max|c|)``.
    max_iter : int, optional
        Iteration cap.

    Returns
    -------
    IntensityResult
        ``result(t)`` evaluates ``μ(t)``; ``result.fd`` is ``W``.

    Raises
    ------
    ValueError
        If ``x`` is empty, not finite or not one-dimensional, lies outside the
        basis domain, or a setting is out of range.

    Warns
    -----
    RuntimeWarning
        If the Newton iteration has not converged after ``max_iter`` steps.

    Examples
    --------
    >>> import numpy as np
    >>> import fdatools as fdt
    >>> from fdatools.density import fit_intensity
    >>> rng = np.random.default_rng(2)
    >>> events = np.cumsum(rng.exponential(scale=0.5, size=200))
    >>> res = fit_intensity(events, penalty=1, lam=10.0)
    >>> round(res.expected_count, 6)
    200.0
    """
    values = _sample(x, "x")
    lam_value = _check_settings(lam, tol, max_iter)
    if basis is None:
        xp = default_namespace()
        window = _resolve_domain(domain, (min(0.0, float(xp.min(values))), float(xp.max(values))))
        basis = _default_basis(window, int(values.shape[0]))
    _check_inside(values, basis, "x")
    problem, coefs, pen, n_iter, converged, grad_norm = _fit(
        "intensity", values, basis, lam_value, penalty, start, tol, max_iter
    )
    xp = problem.xp
    count = problem.integral(coefs)
    return IntensityResult(
        fd=FData(coefs, basis),
        lam=lam_value,
        penalty_matrix=pen,
        criterion=problem.value(coefs),
        log_likelihood=float(xp.sum(problem.data_sum * coefs)) - count,
        expected_count=count,
        gradient_norm=grad_norm,
        n_iter=n_iter,
        converged=converged,
    )
