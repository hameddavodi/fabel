"""Functional bases: B-spline, Fourier, monomial, exponential, power, polygonal.

Every basis is an immutable dataclass, is callable (``b(t, deriv=0)`` returns the
``(len(t), n_basis)`` basis matrix in the caller's array namespace) and knows how
to build its own penalty and Gram matrices.

The mathematics is implemented from first principles: the Cox-de Boor recursion
for B-splines, the analytic Fourier/monomial/exponential/power definitions, and
Gauss-Legendre quadrature of sufficient order for the inner products, following
Ramsay & Silverman, *Functional Data Analysis* (2nd ed.), chapter 3.

Examples
--------
>>> import numpy as np
>>> from fabel import BSpline
>>> b = BSpline(domain=(0.0, 1.0), n_basis=6)
>>> b(np.array([0.0, 0.5, 1.0])).shape
(3, 6)
>>> float(np.max(np.abs(b(np.linspace(0, 1, 11)).sum(axis=1) - 1.0))) < 1e-14
True
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from math import isclose, pi, sqrt
from types import ModuleType
from typing import Any, ClassVar

from fabel import _linalg
from fabel._backend import array_namespace, asarray, default_namespace, to_numpy
from fabel._operator import LDO

__all__ = [
    "BSpline",
    "Basis",
    "Constant",
    "Exponential",
    "Fourier",
    "Monomial",
    "Polygonal",
    "Power",
]

Array = Any

#: Relative tolerance used when deciding whether an evaluation point lies inside
#: the domain, and when comparing two domains for compatibility.
_TOL = 1e-10

#: Order of the B-spline used when a basis product has no exact closed form.
_FALLBACK_PRODUCT_ORDER = 8

#: Gauss-Legendre nodes per panel for penalties without a closed form.
_NUMERIC_QUAD_DEGREE = 12


# --------------------------------------------------------------------------- #
# B-spline kernels (shared by BSpline and Polygonal)
# --------------------------------------------------------------------------- #


def _knot_vector(breaks: tuple[float, ...], order: int) -> tuple[float, ...]:
    """Return the open knot vector: endpoints repeated ``order`` times in total."""
    pad = order - 1
    return (breaks[0],) * pad + breaks + (breaks[-1],) * pad


def _bspline_matrix(
    t: Array, knots: Array, order: int, xp: ModuleType, *, closed_right: bool = True
) -> Array:
    """Cox-de Boor recursion for the order-``order`` basis on ``knots``.

    Returns an ``(n_t, len(knots) - order)`` matrix.  Every knot interval is
    half-open on the right, so the basis is right-continuous at the knots.  When
    ``closed_right`` is true the rightmost knot is additionally counted as part
    of the last non-degenerate interval, so the final basis function equals one
    at the right end of the domain instead of dropping to zero.
    """
    m = knots.shape[0]
    t2 = t[:, None]
    lower = knots[0 : m - 1]
    upper = knots[1:m]
    right = knots[m - 1]
    mat = xp.astype((t2 >= lower) & (t2 < upper), t.dtype)
    if closed_right:
        at_right_end = (t2 >= right) & (lower < upper) & (upper >= right)
        mat = mat + xp.astype(at_right_end, t.dtype)
    for degree in range(2, order + 1):
        ncol = m - degree
        left_knot = knots[0:ncol]
        right_knot = knots[degree : degree + ncol]
        span1 = knots[degree - 1 : degree - 1 + ncol] - left_knot
        span2 = right_knot - knots[1 : 1 + ncol]
        num1 = t2 - left_knot
        num2 = right_knot - t2
        ok1 = span1 > 0.0
        ok2 = span2 > 0.0
        w1 = xp.where(ok1, num1 / xp.where(ok1, span1, xp.ones_like(span1)), num1 * 0.0)
        w2 = xp.where(ok2, num2 / xp.where(ok2, span2, xp.ones_like(span2)), num2 * 0.0)
        mat = w1 * mat[:, 0:ncol] + w2 * mat[:, 1 : ncol + 1]
    return mat


def _bspline_derivative_operator(knots: Array, order: int, deriv: int, xp: ModuleType) -> Array:
    r"""Matrix ``A`` with ``N_order^{(deriv)} = N_{order-deriv} @ A``.

    Uses ``D B_{i,d} = g_i B_{i,d-1} - g_{i+1} B_{i+1,d-1}`` with
    ``g_j = (d - 1) / (knots[j + d - 1] - knots[j])``.
    """
    m = knots.shape[0]
    acc: Array | None = None
    for degree in range(order - deriv + 1, order + 1):
        ncol = m - degree
        nrow = ncol + 1
        gap = knots[degree - 1 : degree - 1 + nrow] - knots[0:nrow]
        ok = gap > 0.0
        g = xp.where(ok, (degree - 1.0) / xp.where(ok, gap, xp.ones_like(gap)), gap * 0.0)
        step = xp.eye(nrow, ncol, k=0, dtype=knots.dtype) * g[None, 0:ncol]
        step = step - xp.eye(nrow, ncol, k=-1, dtype=knots.dtype) * g[None, 1 : ncol + 1]
        acc = step if acc is None else xp.matmul(acc, step)
    return acc


def _bspline_eval(t: Array, knots: Array, order: int, deriv: int, xp: ModuleType) -> Array:
    """Evaluate the ``deriv``-th derivative of the B-spline basis at ``t``.

    Derivative order ``order - 1`` turns the spline into a step function, which
    has no value at a knot other than its right-hand limit.  At the right end of
    the domain there is no interval to the right, so that derivative is reported
    as zero -- the convention R's ``fda`` uses.  Lower derivatives are
    continuous, so the right endpoint is included in the usual way.
    """
    n_basis = knots.shape[0] - order
    if deriv >= order:
        return xp.zeros((t.shape[0], n_basis), dtype=t.dtype)
    closed_right = deriv == 0 or order - deriv > 1
    base = _bspline_matrix(t, knots, order - deriv, xp, closed_right=closed_right)
    if deriv == 0:
        return base
    return xp.matmul(base, _bspline_derivative_operator(knots, order, deriv, xp))


# --------------------------------------------------------------------------- #
# Basis ABC
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, init=False)
class Basis(ABC):
    """Abstract base class for every Fabel basis system.

    A basis is a frozen dataclass, so it is hashable, comparable and safe to use
    as a cache key.  Subclasses provide the evaluation rule; this class provides
    the shared machinery for derivatives, penalties, Gram matrices and products.

    Attributes
    ----------
    domain : tuple of float
        The interval ``(a, b)`` the basis is defined on.
    """

    domain: tuple[float, float]

    #: Whether evaluation outside ``domain`` is an error (piecewise bases).
    _bounded: ClassVar[bool] = False

    # ---------------------------------------------------------------- abstract

    @property
    @abstractmethod
    def n_basis(self) -> int:
        """Number of basis functions."""

    @property
    @abstractmethod
    def names(self) -> tuple[str, ...]:
        """Human-readable name of each basis function."""

    @abstractmethod
    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        """Evaluate the ``deriv``-th derivative at the validated points ``t``."""

    # ------------------------------------------------------------------ public

    def __call__(self, t: Array, deriv: int | LDO = 0) -> Array:
        """Evaluate the basis (or a derivative of it) at ``t``.

        Parameters
        ----------
        t : array_like
            Evaluation points, shape ``(n_t,)``.
        deriv : int or LDO, optional
            Derivative order, or a linear differential operator to apply.

        Returns
        -------
        array
            Basis matrix of shape ``(n_t, n_basis)`` in the namespace of ``t``.

        Raises
        ------
        ValueError
            If ``deriv`` is negative, or if ``t`` leaves the domain of a
            piecewise-defined basis (B-spline, polygonal).

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import Monomial
        >>> Monomial(domain=(0.0, 1.0), n_basis=3)(np.array([2.0])).tolist()
        [[1.0, 2.0, 4.0]]
        """
        xp = array_namespace(t)
        points = asarray(t, xp=xp)
        if len(points.shape) != 1:
            points = xp.reshape(points, (-1,))
        points = self._validate_points(points, xp)
        if isinstance(deriv, LDO):
            return deriv.apply(lambda j: self._evaluate(points, j, xp), points, xp)
        order = int(deriv)
        if order < 0:
            raise ValueError(f"deriv must be non-negative, got {order}")
        return self._evaluate(points, order, xp)

    def penalty(self, op: int | LDO = 2, *, xp: ModuleType | None = None) -> Array:
        r"""Roughness penalty matrix ``R[i, j] = ∫ (L φ_i)(L φ_j) dt``.

        Parameters
        ----------
        op : int or LDO, optional
            Derivative order (default ``2``, the curvature penalty) or a full
            linear differential operator.
        xp : module, optional
            Array namespace for the result.  Defaults to NumPy; the matrix is a
            constant of the basis, so it carries no gradient either way.

        Returns
        -------
        array
            Symmetric positive-semidefinite matrix of shape
            ``(n_basis, n_basis)``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import Fourier
        >>> pen = Fourier(domain=(0.0, 1.0), n_basis=3).penalty(1)
        >>> np.round(np.diag(pen), 6).tolist()
        [0.0, 39.478418, 39.478418]
        """
        target = default_namespace() if xp is None else xp
        operator = op if isinstance(op, LDO) else LDO(int(op))
        if operator.is_derivative:
            key: Any = (type(self).__name__, self)
            slot = operator.order
        else:
            key = (type(self).__name__, self, operator)
            slot = -1
        values = _linalg.cached_gram(key, slot, lambda: self._penalty_matrix(operator))
        return _linalg.as_backend(values, target)

    def gram(self, deriv: int = 0, *, xp: ModuleType | None = None) -> Array:
        r"""Return the cached Gram matrix ``G[i, j] = ∫ φ_i^{(deriv)} φ_j^{(deriv)} dt``.

        Parameters
        ----------
        deriv : int, optional
            Derivative order applied to both factors.
        xp : module, optional
            Array namespace for the result.  Defaults to NumPy.

        Returns
        -------
        array
            Symmetric matrix of shape ``(n_basis, n_basis)``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel import Fourier
        >>> np.round(Fourier(domain=(0.0, 1.0), n_basis=3).gram(), 12).tolist()
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        """
        return self.penalty(int(deriv), xp=xp)

    def __mul__(self, other: Basis) -> Basis:
        """Return a basis whose span contains every pointwise product ``φ_i ψ_j``.

        Parameters
        ----------
        other : Basis
            The second factor; must share this basis's domain.

        Returns
        -------
        Basis
            The product basis.  Exact for spline x spline (merged breaks, order
            ``o1 + o2 - 1``), Fourier x Fourier of equal period, and for the
            monomial / power / exponential families (summed exponents or rates).
            Any other combination falls back to a rich B-spline space.
        """
        if not isinstance(other, Basis):
            return NotImplemented
        if not _same_domain(self.domain, other.domain):
            raise ValueError(f"domains differ: {self.domain} vs {other.domain}")
        return _product_basis(self, other)

    # --------------------------------------------------------------- internals

    def _validate_points(self, t: Array, xp: ModuleType) -> Array:
        """Check ``t`` against the domain and clip rounding errors at the ends."""
        if not self._bounded:
            return t
        lower, upper = self.domain
        slack = _TOL * max(1.0, abs(upper - lower))
        if bool(xp.any((t < lower - slack) | (t > upper + slack))):
            span = to_numpy(t)
            raise ValueError(
                f"evaluation points must lie in {self.domain}; got [{span.min()}, {span.max()}]"
            )
        return xp.clip(t, lower, upper)

    def _natural_breaks(self) -> tuple[float, ...]:
        """Panel boundaries on which the basis is smooth; the domain by default."""
        return self.domain

    def _quadrature_panels(self) -> Array:
        """Refined panel boundaries (NumPy) for penalties without a closed form."""
        xp = default_namespace()
        natural = self._natural_breaks()
        target = max(4 * self.n_basis, 32)
        splits = max(1, -(-target // (len(natural) - 1)))
        pieces = [
            xp.linspace(natural[i], natural[i + 1], splits + 1, dtype=xp.float64)[:-1]
            for i in range(len(natural) - 1)
        ]
        pieces.append(xp.asarray([natural[-1]], dtype=xp.float64))
        return xp.concat(pieces)

    def _penalty_matrix(self, op: LDO) -> Array:
        """Build the penalty matrix as a NumPy array."""
        if op.is_derivative:
            return self._derivative_penalty(op.order)
        return self._quadrature_penalty(op, self._quadrature_panels(), _NUMERIC_QUAD_DEGREE)

    def _derivative_penalty(self, deriv: int) -> Array:
        """Penalty for a plain derivative; overridden with closed forms."""
        return self._quadrature_penalty(LDO(deriv), self._quadrature_panels(), _NUMERIC_QUAD_DEGREE)

    def _quadrature_penalty(self, op: LDO, panels: Array, degree: int) -> Array:
        """Composite Gauss-Legendre approximation of the penalty integral."""
        xp = default_namespace()
        nodes, weights = _linalg.composite_gauss_legendre(panels, degree)
        mat = self(nodes, op)
        raw = xp.matmul(xp.matrix_transpose(mat), weights[:, None] * mat)
        return 0.5 * (raw + xp.matrix_transpose(raw))


# --------------------------------------------------------------------------- #
# B-spline
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, init=False)
class BSpline(Basis):
    """B-spline basis of a given order on a set of break points.

    Replaces R's ``create.bspline.basis``.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)``.  Defaults to ``(0.0, 1.0)``.
    n_basis : int, optional
        Number of basis functions.  Ignored when ``breaks`` is supplied, where
        it becomes ``len(breaks) + order - 2``.  Defaults to ``order``, i.e. a
        single polynomial piece.
    order : int, optional
        Spline order = polynomial degree + 1.  Defaults to ``4`` (cubic).
    breaks : sequence of float, optional
        Strictly increasing break points including both endpoints.  Defaults to
        ``n_basis - order + 2`` equally spaced points.

    Attributes
    ----------
    order : int
        Spline order.
    breaks : tuple of float
        The break points actually used.

    Raises
    ------
    ValueError
        If the arguments are inconsistent, or ``n_basis < order``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> b = BSpline(domain=(0.0, 1.0), n_basis=5, order=3)
    >>> b.n_basis
    5
    >>> b.breaks
    (0.0, 0.3333333333333333, 0.6666666666666666, 1.0)
    >>> float(b(np.array([1.0]))[0, -1])
    1.0
    """

    order: int
    breaks: tuple[float, ...]

    _bounded: ClassVar[bool] = True

    def __init__(
        self,
        domain: tuple[float, float] = (0.0, 1.0),
        n_basis: int | None = None,
        order: int = 4,
        breaks: Sequence[float] | None = None,
    ) -> None:
        span = _clean_domain(domain)
        k = int(order)
        if k < 1:
            raise ValueError(f"order must be at least 1, got {k}")
        if breaks is not None:
            knots = _clean_breaks(breaks, span, k)
            derived = len(knots) + k - 2
            if n_basis is not None and int(n_basis) != derived:
                raise ValueError(f"n_basis={n_basis} conflicts with len(breaks)+order-2={derived}")
        else:
            wanted = k if n_basis is None else int(n_basis)
            if wanted < k:
                raise ValueError(f"n_basis must be at least order={k}, got {wanted}")
            xp = default_namespace()
            grid = xp.linspace(span[0], span[1], wanted - k + 2, dtype=xp.float64)
            knots = tuple(float(v) for v in grid)
        object.__setattr__(self, "domain", span)
        object.__setattr__(self, "order", k)
        object.__setattr__(self, "breaks", knots)

    @property
    def n_basis(self) -> int:
        """Number of B-spline basis functions, ``len(breaks) + order - 2``.

        Returns
        -------
        int
            The basis dimension.

        Examples
        --------
        >>> from fabel import BSpline
        >>> BSpline(breaks=[0.0, 0.5, 1.0], order=4).n_basis
        5
        """
        return len(self.breaks) + self.order - 2

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``bspl<order>.<i>`` following R's ``fda`` convention.

        Returns
        -------
        tuple of str
            One name per basis function.

        Examples
        --------
        >>> from fabel import BSpline
        >>> BSpline(n_basis=4).names
        ('bspl4.1', 'bspl4.2', 'bspl4.3', 'bspl4.4')
        """
        return tuple(f"bspl{self.order}.{i + 1}" for i in range(self.n_basis))

    @property
    def knots(self) -> tuple[float, ...]:
        """The open knot vector, endpoints repeated ``order`` times.

        Returns
        -------
        tuple of float
            Knot sequence of length ``n_basis + order``.

        Examples
        --------
        >>> from fabel import BSpline
        >>> BSpline(breaks=[0.0, 1.0], order=2).knots
        (0.0, 0.0, 1.0, 1.0)
        """
        return _knot_vector(self.breaks, self.order)

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        knots = asarray(self.knots, xp=xp)
        return _bspline_eval(t, knots, self.order, deriv, xp)

    def _natural_breaks(self) -> tuple[float, ...]:
        return self.breaks

    def _derivative_penalty(self, deriv: int) -> Array:
        # The integrand is a polynomial of degree 2 * (order - 1 - deriv) on each
        # knot interval, so an `order`-point Gauss-Legendre rule per interval is
        # exact (it integrates degree 2 * order - 1 exactly).
        xp = default_namespace()
        panels = xp.asarray(self.breaks, dtype=xp.float64)
        raw = self._quadrature_penalty(LDO(deriv), panels, self.order)
        return _band_mask(raw, self.order - 1, xp)


# --------------------------------------------------------------------------- #
# Fourier
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, init=False)
class Fourier(Basis):
    r"""Fourier basis, orthonormal over one period.

    The functions are ``1/sqrt(P)``, then ``sqrt(2/P) sin(h ω t)`` and
    ``sqrt(2/P) cos(h ω t)`` for ``h = 1 .. H`` with ``ω = 2π / P``.  Replaces
    R's ``create.fourier.basis``.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)``.  Defaults to ``(0.0, 1.0)``.
    n_basis : int, optional
        Number of basis functions.  An even value is raised by one so that the
        sine/cosine pairs are complete, matching R; read back the effective size
        from :attr:`n_basis`.
    period : float, optional
        Period ``P``.  Defaults to the width of the domain.

    Attributes
    ----------
    n_harmonics : int
        Number of sine/cosine pairs ``H``.
    period : float
        The period.

    Raises
    ------
    ValueError
        If ``n_basis < 1`` or ``period <= 0``.

    Examples
    --------
    >>> from fabel import Fourier
    >>> Fourier(domain=(0.0, 365.0), n_basis=4).n_basis
    5
    >>> Fourier(domain=(0.0, 365.0), n_basis=5).period
    365.0
    """

    n_harmonics: int
    period: float

    def __init__(
        self,
        domain: tuple[float, float] = (0.0, 1.0),
        n_basis: int = 3,
        period: float | None = None,
    ) -> None:
        span = _clean_domain(domain)
        size = int(n_basis)
        if size < 1:
            raise ValueError(f"n_basis must be at least 1, got {size}")
        if size % 2 == 0:
            size += 1
        length = (span[1] - span[0]) if period is None else float(period)
        if length <= 0.0:
            raise ValueError(f"period must be positive, got {length}")
        object.__setattr__(self, "domain", span)
        object.__setattr__(self, "n_harmonics", (size - 1) // 2)
        object.__setattr__(self, "period", length)

    @property
    def n_basis(self) -> int:
        """Number of basis functions, always odd: ``2 * n_harmonics + 1``.

        Returns
        -------
        int
            The basis dimension.

        Examples
        --------
        >>> from fabel import Fourier
        >>> Fourier(n_basis=6).n_basis
        7
        """
        return 2 * self.n_harmonics + 1

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``const``, ``sin1``, ``cos1``, ... following R's convention.

        Returns
        -------
        tuple of str
            One name per basis function.

        Examples
        --------
        >>> from fabel import Fourier
        >>> Fourier(n_basis=3).names
        ('const', 'sin1', 'cos1')
        """
        out = ["const"]
        for h in range(1, self.n_harmonics + 1):
            out += [f"sin{h}", f"cos{h}"]
        return tuple(out)

    def _terms(self, deriv: int) -> tuple[list[float], list[float], list[bool]]:
        r"""Write the ``deriv``-th derivative as ``amp * {sin, cos}(freq * t)``.

        The constant function is represented as ``cos(0 * t)``, which keeps the
        product integrals below uniform across all basis functions.
        """
        omega = 2.0 * pi / self.period
        amps = [1.0 / sqrt(self.period) if deriv == 0 else 0.0]
        freqs = [0.0]
        sines = [False]
        scale = sqrt(2.0 / self.period)
        for h in range(1, self.n_harmonics + 1):
            rate = h * omega
            size = scale * rate**deriv
            # d^m/dt^m sin(a t) = a^m sin(a t + m π / 2)
            sin_sign, sin_is_sin = ((1.0, True), (1.0, False), (-1.0, True), (-1.0, False))[
                deriv % 4
            ]
            # d^m/dt^m cos(a t) = a^m cos(a t + m π / 2)
            cos_sign, cos_is_sin = ((1.0, False), (-1.0, True), (-1.0, False), (1.0, True))[
                deriv % 4
            ]
            amps += [sin_sign * size, cos_sign * size]
            freqs += [rate, rate]
            sines += [sin_is_sin, cos_is_sin]
        return amps, freqs, sines

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        amps, freqs, sines = self._terms(deriv)
        amp = asarray(amps, xp=xp)
        freq = asarray(freqs, xp=xp)
        is_sin = xp.asarray(sines)
        arg = t[:, None] * freq[None, :]
        return amp[None, :] * xp.where(is_sin[None, :], xp.sin(arg), xp.cos(arg))

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        amps, freqs, sines = self._terms(deriv)
        amp = xp.asarray(amps, dtype=xp.float64)
        freq = xp.asarray(freqs, dtype=xp.float64)
        is_sin = xp.asarray(sines)
        lower, upper = self.domain
        harmonic = xp.asarray([0, *[h for h in range(1, self.n_harmonics + 1) for _ in (0, 1)]])
        # Every frequency is an integer multiple of 2*pi/period, so over a whole
        # number of periods the harmonics are exactly orthogonal.  Saying so
        # explicitly avoids cancellation noise of order eps * n_harmonics**deriv.
        whole = _whole_periods(lower, upper, self.period)
        p = freq[:, None]
        q = freq[None, :]
        si = is_sin[:, None]
        sj = is_sin[None, :]
        diff_k = harmonic[:, None] - harmonic[None, :]
        sum_k = harmonic[:, None] + harmonic[None, :]
        cos_diff = _cos_integral(xp, p - q, lower, upper, diff_k if whole else None)
        cos_sum = _cos_integral(xp, p + q, lower, upper, sum_k if whole else None)
        sin_diff = _sin_integral(xp, p - q, lower, upper, diff_k if whole else None)
        sin_sum = _sin_integral(xp, p + q, lower, upper, sum_k if whole else None)
        sin_sin = 0.5 * (cos_diff - cos_sum)
        cos_cos = 0.5 * (cos_diff + cos_sum)
        sin_cos = 0.5 * (sin_sum + sin_diff)
        cos_sin = 0.5 * (sin_sum - sin_diff)
        integral = xp.where(
            si & sj, sin_sin, xp.where(si & ~sj, sin_cos, xp.where(~si & sj, cos_sin, cos_cos))
        )
        raw = amp[:, None] * amp[None, :] * integral
        return 0.5 * (raw + xp.matrix_transpose(raw))


def _whole_periods(lower: float, upper: float, period: float) -> bool:
    """Return whether ``[lower, upper]`` spans an exact whole number of periods."""
    count = (upper - lower) / period
    return bool(abs(count - round(count)) <= _TOL * max(1.0, abs(count)))


def _cos_integral(
    xp: ModuleType, rate: Array, lower: float, upper: float, harmonics: Array | None = None
) -> Array:
    """Return ``∫_lower^upper cos(rate t) dt`` elementwise, safe at ``rate = 0``.

    When ``harmonics`` is given the integration range is known to span a whole
    number of periods and ``rate`` is ``harmonics`` times the fundamental, so the
    integral is exactly zero away from the zeroth harmonic.
    """
    zero_rate = xp.abs(rate) <= _TOL if harmonics is None else harmonics == 0
    if harmonics is not None:
        return xp.where(zero_rate, xp.full_like(rate, upper - lower), xp.zeros_like(rate))
    safe = xp.where(zero_rate, xp.ones_like(rate), rate)
    return xp.where(
        zero_rate,
        xp.full_like(rate, upper - lower),
        (xp.sin(safe * upper) - xp.sin(safe * lower)) / safe,
    )


def _sin_integral(
    xp: ModuleType, rate: Array, lower: float, upper: float, harmonics: Array | None = None
) -> Array:
    """Return ``∫_lower^upper sin(rate t) dt`` elementwise, safe at ``rate = 0``.

    ``harmonics`` carries the same meaning as in :func:`_cos_integral`; over a
    whole number of periods every sine integrates to exactly zero.
    """
    if harmonics is not None:
        return xp.zeros_like(rate)
    zero_rate = xp.abs(rate) <= _TOL
    safe = xp.where(zero_rate, xp.ones_like(rate), rate)
    return xp.where(
        zero_rate,
        xp.zeros_like(rate),
        (xp.cos(safe * lower) - xp.cos(safe * upper)) / safe,
    )


# --------------------------------------------------------------------------- #
# Monomial / Power
# --------------------------------------------------------------------------- #


def _falling_factorial(exponent: float, deriv: int) -> float:
    """Return ``e (e-1) ... (e-deriv+1)``, the coefficient of ``D^deriv t^e``."""
    out = 1.0
    for r in range(deriv):
        out *= exponent - r
    return out


@dataclass(frozen=True, init=False)
class Monomial(Basis):
    """Monomial basis ``t^e`` for non-negative integer exponents.

    Replaces R's ``create.monomial.basis``.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)``.  Defaults to ``(0.0, 1.0)``.
    n_basis : int, optional
        Number of basis functions; ignored when ``exponents`` is given.
        Defaults to ``2`` (constant and linear).
    exponents : sequence of int, optional
        Exponents to use.  Defaults to ``0 .. n_basis - 1``.

    Attributes
    ----------
    exponents : tuple of int
        The exponents.

    Raises
    ------
    ValueError
        If exponents are negative or repeated, or conflict with ``n_basis``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import Monomial
    >>> Monomial(n_basis=3)(np.array([3.0]), deriv=1).tolist()
    [[0.0, 1.0, 6.0]]
    """

    exponents: tuple[int, ...]

    def __init__(
        self,
        domain: tuple[float, float] = (0.0, 1.0),
        n_basis: int = 2,
        exponents: Sequence[int] | None = None,
    ) -> None:
        span = _clean_domain(domain)
        if exponents is None:
            size = int(n_basis)
            if size < 1:
                raise ValueError(f"n_basis must be at least 1, got {size}")
            powers = tuple(range(size))
        else:
            powers = tuple(int(e) for e in exponents)
            if any(e < 0 for e in powers):
                raise ValueError(f"monomial exponents must be non-negative, got {powers}")
            if len(set(powers)) != len(powers):
                raise ValueError(f"monomial exponents must be distinct, got {powers}")
        object.__setattr__(self, "domain", span)
        object.__setattr__(self, "exponents", powers)

    @property
    def n_basis(self) -> int:
        """Number of monomials.

        Returns
        -------
        int
            ``len(self.exponents)``.

        Examples
        --------
        >>> from fabel import Monomial
        >>> Monomial(exponents=[0, 2, 4]).n_basis
        3
        """
        return len(self.exponents)

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``monomial<e>``.

        Returns
        -------
        tuple of str
            One name per basis function.

        Examples
        --------
        >>> from fabel import Monomial
        >>> Monomial(n_basis=2).names
        ('monomial0', 'monomial1')
        """
        return tuple(f"monomial{e}" for e in self.exponents)

    def _terms(self, deriv: int) -> tuple[list[float], list[float]]:
        coefs: list[float] = []
        powers: list[float] = []
        for e in self.exponents:
            if deriv > e:
                coefs.append(0.0)
                powers.append(0.0)
            else:
                coefs.append(_falling_factorial(float(e), deriv))
                powers.append(float(e - deriv))
        return coefs, powers

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        coefs, powers = self._terms(deriv)
        coef = asarray(coefs, xp=xp)
        power = asarray(powers, xp=xp)
        return coef[None, :] * xp.pow(t[:, None], power[None, :])

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        coefs, powers = self._terms(deriv)
        coef = xp.asarray(coefs, dtype=xp.float64)
        power = xp.asarray(powers, dtype=xp.float64)
        lower, upper = self.domain
        total = power[:, None] + power[None, :] + 1.0
        integral = (upper**total - lower**total) / total
        raw = coef[:, None] * coef[None, :] * integral
        return 0.5 * (raw + xp.matrix_transpose(raw))


@dataclass(frozen=True, init=False)
class Power(Basis):
    """Power basis ``t^e`` for arbitrary real exponents.

    Replaces R's ``create.power.basis``.  The domain must not contain negative
    values, and must exclude zero whenever a negative power appears.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)`` with ``a >= 0``.  Defaults to ``(1.0, 2.0)``.
    exponents : sequence of float, optional
        The exponents.  Defaults to ``(0.0, 1.0)``.

    Attributes
    ----------
    exponents : tuple of float
        The exponents.

    Raises
    ------
    ValueError
        If the domain is negative, the exponents repeat, or a negative exponent
        is combined with a domain touching zero.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import Power
    >>> Power(domain=(1.0, 4.0), exponents=[0.5])(np.array([4.0])).tolist()
    [[2.0]]
    """

    exponents: tuple[float, ...]

    def __init__(
        self,
        domain: tuple[float, float] = (1.0, 2.0),
        exponents: Sequence[float] = (0.0, 1.0),
    ) -> None:
        span = _clean_domain(domain)
        if span[0] < 0.0:
            raise ValueError(f"a power basis needs a non-negative domain, got {span}")
        powers = tuple(float(e) for e in exponents)
        if not powers:
            raise ValueError("a power basis needs at least one exponent")
        if len(set(powers)) != len(powers):
            raise ValueError(f"power exponents must be distinct, got {powers}")
        if span[0] == 0.0 and any(e < 0.0 for e in powers):
            raise ValueError(f"negative exponents {powers} need a domain excluding zero")
        object.__setattr__(self, "domain", span)
        object.__setattr__(self, "exponents", powers)

    @property
    def n_basis(self) -> int:
        """Number of power functions.

        Returns
        -------
        int
            ``len(self.exponents)``.

        Examples
        --------
        >>> from fabel import Power
        >>> Power(exponents=[0.0, 0.5, 1.0]).n_basis
        3
        """
        return len(self.exponents)

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``power<e>``.

        Returns
        -------
        tuple of str
            One name per basis function.

        Examples
        --------
        >>> from fabel import Power
        >>> Power(exponents=[0.0, 0.5]).names
        ('power0.0', 'power0.5')
        """
        return tuple(f"power{e}" for e in self.exponents)

    def _terms(self, deriv: int) -> tuple[list[float], list[float]]:
        coefs = [_falling_factorial(e, deriv) for e in self.exponents]
        powers = [e - deriv for e in self.exponents]
        return coefs, powers

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        coefs, powers = self._terms(deriv)
        coef = asarray(coefs, xp=xp)
        power = asarray(powers, xp=xp)
        return coef[None, :] * xp.pow(t[:, None], power[None, :])

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        coefs, powers = self._terms(deriv)
        coef = xp.asarray(coefs, dtype=xp.float64)
        power = xp.asarray(powers, dtype=xp.float64)
        lower, upper = self.domain
        total = power[:, None] + power[None, :] + 1.0
        singular = xp.abs(total) <= _TOL
        safe = xp.where(singular, xp.ones_like(total), total)
        if lower <= 0.0 and bool(xp.any(singular | (total < 0.0))):
            raise ValueError(
                f"penalty of order {deriv} diverges on domain {self.domain}; "
                "use a domain bounded away from zero"
            )
        log_case = xp.full_like(total, _safe_log_ratio(lower, upper))
        integral = xp.where(singular, log_case, (upper**safe - lower**safe) / safe)
        raw = coef[:, None] * coef[None, :] * integral
        return 0.5 * (raw + xp.matrix_transpose(raw))


def _safe_log_ratio(lower: float, upper: float) -> float:
    """Return ``log(upper / lower)``, or ``0`` when the log branch is unused."""
    if lower <= 0.0:
        return 0.0
    from math import log

    return log(upper) - log(lower)


# --------------------------------------------------------------------------- #
# Exponential
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, init=False)
class Exponential(Basis):
    """Exponential basis ``exp(r t)`` for a set of rates.

    Replaces R's ``create.exponential.basis``.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)``.  Defaults to ``(0.0, 1.0)``.
    rates : sequence of float, optional
        The rates ``r``.  Defaults to ``(0.0, 1.0)``.

    Attributes
    ----------
    rates : tuple of float
        The rates.

    Raises
    ------
    ValueError
        If the rates are empty or repeated.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import Exponential
    >>> Exponential(rates=[0.0, 1.0])(np.array([0.0])).tolist()
    [[1.0, 1.0]]
    """

    rates: tuple[float, ...]

    def __init__(
        self,
        domain: tuple[float, float] = (0.0, 1.0),
        rates: Sequence[float] = (0.0, 1.0),
    ) -> None:
        span = _clean_domain(domain)
        values = tuple(float(r) for r in rates)
        if not values:
            raise ValueError("an exponential basis needs at least one rate")
        if len(set(values)) != len(values):
            raise ValueError(f"exponential rates must be distinct, got {values}")
        object.__setattr__(self, "domain", span)
        object.__setattr__(self, "rates", values)

    @property
    def n_basis(self) -> int:
        """Number of exponentials.

        Returns
        -------
        int
            ``len(self.rates)``.

        Examples
        --------
        >>> from fabel import Exponential
        >>> Exponential(rates=[0.0, 1.0, 2.0]).n_basis
        3
        """
        return len(self.rates)

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``exp<rate>``.

        Returns
        -------
        tuple of str
            One name per basis function.

        Examples
        --------
        >>> from fabel import Exponential
        >>> Exponential(rates=[0.0, 1.0]).names
        ('exp0.0', 'exp1.0')
        """
        return tuple(f"exp{r}" for r in self.rates)

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        rate = asarray(self.rates, xp=xp)
        scale = asarray([r**deriv for r in self.rates], xp=xp)
        return scale[None, :] * xp.exp(t[:, None] * rate[None, :])

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        rate = xp.asarray(self.rates, dtype=xp.float64)
        scale = xp.asarray([r**deriv for r in self.rates], dtype=xp.float64)
        lower, upper = self.domain
        total = rate[:, None] + rate[None, :]
        singular = xp.abs(total) <= _TOL
        safe = xp.where(singular, xp.ones_like(total), total)
        integral = xp.where(
            singular,
            xp.full_like(total, upper - lower),
            (xp.exp(safe * upper) - xp.exp(safe * lower)) / safe,
        )
        raw = scale[:, None] * scale[None, :] * integral
        return 0.5 * (raw + xp.matrix_transpose(raw))


# --------------------------------------------------------------------------- #
# Constant / Polygonal
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, init=False)
class Constant(Basis):
    """The one-function basis ``φ(t) = 1``.

    Replaces R's ``create.constant.basis``.

    Parameters
    ----------
    domain : tuple of float, optional
        Interval ``(a, b)``.  Defaults to ``(0.0, 1.0)``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import Constant
    >>> Constant(domain=(0.0, 2.0)).gram().tolist()
    [[2.0]]
    """

    def __init__(self, domain: tuple[float, float] = (0.0, 1.0)) -> None:
        object.__setattr__(self, "domain", _clean_domain(domain))

    @property
    def n_basis(self) -> int:
        """Always one.

        Returns
        -------
        int
            ``1``.

        Examples
        --------
        >>> from fabel import Constant
        >>> Constant().n_basis
        1
        """
        return 1

    @property
    def names(self) -> tuple[str, ...]:
        """The single name ``const``.

        Returns
        -------
        tuple of str
            ``("const",)``.

        Examples
        --------
        >>> from fabel import Constant
        >>> Constant().names
        ('const',)
        """
        return ("const",)

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        if deriv == 0:
            return xp.ones((t.shape[0], 1), dtype=t.dtype)
        return xp.zeros((t.shape[0], 1), dtype=t.dtype)

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        width = 0.0 if deriv else self.domain[1] - self.domain[0]
        return xp.asarray([[width]], dtype=xp.float64)


@dataclass(frozen=True, init=False)
class Polygonal(Basis):
    """Piecewise-linear ("hat function") basis on a set of arguments.

    Replaces R's ``create.polygonal.basis``.  Mathematically this is the
    order-2 B-spline basis with breaks at ``argvals``, and it is evaluated with
    the same Cox-de Boor kernel.

    Parameters
    ----------
    argvals : sequence of float, optional
        Strictly increasing vertex locations; the domain is
        ``(argvals[0], argvals[-1])``.  Defaults to ``(0.0, 0.5, 1.0)``.

    Attributes
    ----------
    argvals : tuple of float
        The vertices.

    Raises
    ------
    ValueError
        If fewer than two arguments are given, or they are not increasing.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import Polygonal
    >>> Polygonal([0.0, 1.0, 2.0])(np.array([0.5])).tolist()
    [[0.5, 0.5, 0.0]]
    """

    argvals: tuple[float, ...]

    _bounded: ClassVar[bool] = True

    def __init__(self, argvals: Sequence[float] = (0.0, 0.5, 1.0)) -> None:
        points = tuple(float(v) for v in argvals)
        if len(points) < 2:
            raise ValueError("a polygonal basis needs at least two arguments")
        if any(b <= a for a, b in pairwise(points)):
            raise ValueError(f"argvals must be strictly increasing, got {points}")
        object.__setattr__(self, "domain", (points[0], points[-1]))
        object.__setattr__(self, "argvals", points)

    @property
    def n_basis(self) -> int:
        """One basis function per vertex.

        Returns
        -------
        int
            ``len(self.argvals)``.

        Examples
        --------
        >>> from fabel import Polygonal
        >>> Polygonal([0.0, 1.0, 2.0]).n_basis
        3
        """
        return len(self.argvals)

    @property
    def names(self) -> tuple[str, ...]:
        """Names ``polyg<i>``.

        Returns
        -------
        tuple of str
            One name per vertex.

        Examples
        --------
        >>> from fabel import Polygonal
        >>> Polygonal([0.0, 1.0]).names
        ('polyg1', 'polyg2')
        """
        return tuple(f"polyg{i + 1}" for i in range(self.n_basis))

    @property
    def order(self) -> int:
        """The equivalent B-spline order, always ``2``.

        Returns
        -------
        int
            ``2``.

        Examples
        --------
        >>> from fabel import Polygonal
        >>> Polygonal([0.0, 1.0]).order
        2
        """
        return 2

    def _evaluate(self, t: Array, deriv: int, xp: ModuleType) -> Array:
        knots = asarray(_knot_vector(self.argvals, 2), xp=xp)
        return _bspline_eval(t, knots, 2, deriv, xp)

    def _natural_breaks(self) -> tuple[float, ...]:
        return self.argvals

    def _derivative_penalty(self, deriv: int) -> Array:
        xp = default_namespace()
        panels = xp.asarray(self.argvals, dtype=xp.float64)
        raw = self._quadrature_penalty(LDO(deriv), panels, 2)
        return _band_mask(raw, 1, xp)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _band_mask(a: Array, bandwidth: int, xp: ModuleType) -> Array:
    """Zero entries further than ``bandwidth`` off the diagonal.

    Two B-spline basis functions whose indices differ by at least the spline
    order have disjoint supports, so their inner product is exactly zero.  The
    quadrature returns rounding noise there instead; this removes it.
    """
    n = a.shape[0]
    index = xp.arange(n)
    keep = xp.abs(index[:, None] - index[None, :]) <= bandwidth
    return xp.where(keep, a, xp.zeros_like(a))


def _clean_domain(domain: Sequence[float]) -> tuple[float, float]:
    """Validate and normalise a domain to a ``(float, float)`` tuple."""
    values = tuple(float(v) for v in domain)
    if len(values) != 2:
        raise ValueError(f"domain must have two entries, got {domain!r}")
    if values[1] <= values[0]:
        raise ValueError(f"domain must be increasing, got {values}")
    return (values[0], values[1])


def _clean_breaks(
    breaks: Sequence[float], domain: tuple[float, float], order: int
) -> tuple[float, ...]:
    """Validate break points against ``domain`` and return them as a tuple.

    Interior breaks may repeat: a break of multiplicity ``m`` drops the spline's
    smoothness there to ``C^(order - 1 - m)``, which is how the exact product of
    two spline spaces is expressed.  The two endpoints must each appear once.
    """
    values = tuple(float(v) for v in breaks)
    if len(values) < 2:
        raise ValueError("breaks must contain at least the two endpoints")
    if any(b < a for a, b in pairwise(values)):
        raise ValueError(f"breaks must be non-decreasing, got {values}")
    if values[0] == values[1] or values[-2] == values[-1]:
        raise ValueError(f"the end breaks must not repeat, got {values}")
    for value in set(values[1:-1]):
        if values.count(value) > order - 1:
            raise ValueError(
                f"break {value} repeats {values.count(value)} times, "
                f"which exceeds the limit of order - 1 = {order - 1}"
            )
    if not isclose(values[0], domain[0], rel_tol=_TOL, abs_tol=_TOL) or not isclose(
        values[-1], domain[1], rel_tol=_TOL, abs_tol=_TOL
    ):
        raise ValueError(f"breaks {values[0], values[-1]} must span the domain {domain}")
    return (domain[0], *values[1:-1], domain[1])


def _same_domain(left: tuple[float, float], right: tuple[float, float]) -> bool:
    """Return whether two domains agree to within the module tolerance."""
    scale = max(1.0, abs(left[1] - left[0]))
    return all(abs(a - b) <= _TOL * scale for a, b in zip(left, right, strict=True))


def _spline_spec(basis: Basis) -> tuple[int, tuple[float, ...]] | None:
    """Return ``(order, breaks)`` if ``basis`` is a spline space, else ``None``."""
    if isinstance(basis, BSpline):
        return basis.order, basis.breaks
    if isinstance(basis, Polygonal):
        return 2, basis.argvals
    return None


def _smoothness(spec: tuple[int, tuple[float, ...]], break_point: float) -> int:
    """Return the continuity class of a spline space at ``break_point``.

    A break of multiplicity ``m`` in an order-``k`` space leaves the spline
    ``k - 1 - m`` times continuously differentiable there; away from every break
    it is a polynomial, reported here as a large finite number.
    """
    order, breaks = spec
    multiplicity = breaks[1:-1].count(break_point)
    if multiplicity == 0:
        return order  # smooth enough to impose no constraint
    return order - 1 - multiplicity


def _product_breaks(
    left: tuple[int, tuple[float, ...]],
    right: tuple[int, tuple[float, ...]],
    order: int,
) -> tuple[float, ...]:
    """Return the break sequence whose spline space holds the product exactly.

    The product of two piecewise polynomials is ``C^min(s1, s2)`` at a shared
    break, so the product space needs multiplicity ``order - 1 - min(s1, s2)``
    there.
    """
    interior = sorted(set(left[1][1:-1]) | set(right[1][1:-1]))
    out: list[float] = [left[1][0]]
    for point in interior:
        smoothness = min(_smoothness(left, point), _smoothness(right, point))
        out += [point] * max(1, min(order - 1, order - 1 - smoothness))
    out.append(left[1][-1])
    return tuple(out)


def _product_basis(left: Basis, right: Basis) -> Basis:
    """Return a basis spanning every product of a ``left`` and a ``right`` function."""
    if isinstance(left, Constant):
        return right
    if isinstance(right, Constant):
        return left
    left_spline = _spline_spec(left)
    right_spline = _spline_spec(right)
    if left_spline is not None and right_spline is not None:
        order = left_spline[0] + right_spline[0] - 1
        return BSpline(
            domain=left.domain,
            order=order,
            breaks=_product_breaks(left_spline, right_spline, order),
        )
    if (
        isinstance(left, Fourier)
        and isinstance(right, Fourier)
        and isclose(left.period, right.period, rel_tol=_TOL)
    ):
        return Fourier(
            domain=left.domain,
            n_basis=2 * max(left.n_basis, right.n_basis) - 1,
            period=left.period,
        )
    if isinstance(left, Monomial) and isinstance(right, Monomial):
        sums = sorted({a + b for a in left.exponents for b in right.exponents})
        return Monomial(domain=left.domain, exponents=sums)
    if isinstance(left, Power) and isinstance(right, Power):
        powers = sorted({a + b for a in left.exponents for b in right.exponents})
        return Power(domain=left.domain, exponents=powers)
    if isinstance(left, Exponential) and isinstance(right, Exponential):
        rates = sorted({a + b for a in left.rates for b in right.rates})
        return Exponential(domain=left.domain, rates=rates)
    size = max(left.n_basis + right.n_basis, _FALLBACK_PRODUCT_ORDER)
    return BSpline(domain=left.domain, n_basis=size, order=_FALLBACK_PRODUCT_ORDER)
