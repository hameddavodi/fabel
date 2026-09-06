"""Functional data objects: :class:`FData`, :class:`BiFData` and :func:`inprod`.

An :class:`FData` pairs a coefficient array with a :class:`~fabel.basis.Basis`;
everything else -- evaluation, derivatives, statistics, arithmetic and inner
products -- follows from that pair.  All operations stay inside the array
namespace of their inputs, so NumPy arrays give NumPy results and PyTorch
tensors give differentiable PyTorch results.
"""

from __future__ import annotations

import operator
from collections.abc import Iterator
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

from fabel import _linalg
from fabel._backend import array_namespace, asarray, default_namespace, result_namespace, to_numpy
from fabel._operator import LDO
from fabel._plot import PlotMixin
from fabel.basis import Basis, Constant, _same_domain

if TYPE_CHECKING:  # pragma: no cover - typing only
    import numpy as np

    Array = Any
    NDArray = np.ndarray[Any, np.dtype[Any]]
else:
    Array = Any
    NDArray = Any

__all__ = ["LDO", "BiFData", "FData", "inprod"]

#: Gauss-Legendre nodes per quadrature panel.  A degree-``d`` rule is exact for
#: polynomials of degree ``2 d - 1``, so 12 covers every spline product that
#: arises from bases of order up to 12.
_QUAD_DEGREE = 12

#: Panels used when a projection has to resolve a smooth (non-piecewise) basis.
_MIN_PANELS = 16

#: Break refinement applied when a fractional power has to be projected.
_POWER_REFINEMENT = 8

#: Order of the spline used to hold a fractional power.
_POWER_ORDER = 6

#: Largest projection residual accepted when a scalar is expanded in a basis.
_CONSTANT_TOL = 1e-8


def _quadrature(*bases: Basis) -> tuple[NDArray, NDArray]:
    """Return composite Gauss-Legendre nodes and weights for ``bases``.

    The panels are the union of every basis's own break points, refined until
    there are enough of them to resolve the richest basis involved.  For
    piecewise polynomials this makes the rule exact.
    """
    xp = default_namespace()
    natural = sorted({float(b) for basis in bases for b in basis._natural_breaks()})
    target = max(_MIN_PANELS, 2 * sum(basis.n_basis for basis in bases))
    per_interval = max(1, -(-target // (len(natural) - 1)))
    pieces = [
        xp.linspace(natural[i], natural[i + 1], per_interval + 1, dtype=xp.float64)[:-1]
        for i in range(len(natural) - 1)
    ]
    pieces.append(xp.asarray([natural[-1]], dtype=xp.float64))
    return _linalg.composite_gauss_legendre(xp.concat(pieces), _QUAD_DEGREE)


def _cross_gram(left: Basis, right: Basis, op1: LDO, op2: LDO) -> NDArray:
    """Return ``M[i, j] = integral of (op1 phi_i) (op2 psi_j)`` as a NumPy array."""
    xp = default_namespace()
    nodes, weights = _quadrature(left, right)
    a = left(nodes, op1)
    b = right(nodes, op2)
    return cast("NDArray", xp.matmul(xp.matrix_transpose(a), weights[:, None] * b))


def _project(basis: Basis, nodes: NDArray, weights: NDArray, values: Array) -> Array:
    """Return the coefficients of the L2 projection of ``values`` onto ``basis``.

    The design matrix, quadrature weights and Gram matrix depend only on the
    basis, so they are built once in NumPy and converted into the namespace of
    ``values``.  Every step that touches ``values`` then stays in that
    namespace, which keeps a torch autograd graph intact.
    """
    xp = array_namespace(values)
    mat = asarray(basis(nodes), xp=xp)
    quad = asarray(weights, xp=xp)
    flat = xp.reshape(values, (values.shape[0], -1))
    rhs = xp.matmul(xp.matrix_transpose(mat), quad[:, None] * flat)
    coefs = _linalg.solve_spd(asarray(basis.gram(), xp=xp), rhs)
    return xp.reshape(coefs, (basis.n_basis, *values.shape[1:]))


def _curve_positions(index: Any, n_curves: int) -> list[int] | None:
    """Resolve a curve selector to a list of positions, or ``None`` for a scalar.

    Accepts Python and NumPy integers, sequences and arrays of integers, and
    boolean masks.  Negative positions count from the end.

    Raises
    ------
    IndexError
        If a position is out of range, or a boolean mask has the wrong length.
    """
    try:
        operator.index(index)
    except TypeError:
        pass
    else:
        return None
    entries = [entry.item() if hasattr(entry, "item") else entry for entry in index]
    if entries and all(isinstance(entry, bool) for entry in entries):
        if len(entries) != n_curves:
            raise IndexError(f"mask of length {len(entries)} for {n_curves} curves")
        return [i for i, keep in enumerate(entries) if keep]
    positions = [operator.index(entry) for entry in entries]
    resolved = [place + n_curves if place < 0 else place for place in positions]
    if any(not 0 <= place < n_curves for place in resolved):
        raise IndexError(f"curves {positions} out of range for {n_curves} curves")
    return resolved


def _promote(coefs: Array, value: Any) -> tuple[Array, Any]:
    """Return ``coefs`` and ``value`` in one namespace, ready to combine.

    An array scalar keeps its identity -- and so its autograd graph -- and the
    coefficients move to the promoted namespace to meet it; anything else
    becomes a plain float and leaves the coefficients untouched.
    """
    if not _is_array(value):
        return coefs, float(value)
    xp = result_namespace(coefs, value)
    return asarray(coefs, xp=xp), asarray(value, xp=xp)


def _is_array(value: Any) -> bool:
    """Return whether ``value`` is an array API object."""
    return hasattr(value, "__array_namespace__") or hasattr(value, "__array__")


def _as_operator(op: int | LDO) -> LDO:
    """Coerce an integer derivative order to an :class:`LDO`."""
    return op if isinstance(op, LDO) else LDO(int(op))


def _refined_spline(basis: Basis) -> Basis:
    """Return a spline basis fine enough to hold a nonlinear function of ``basis``."""
    from fabel.basis import BSpline

    xp = default_namespace()
    natural = sorted(set(basis._natural_breaks()))
    pieces = [
        xp.linspace(natural[i], natural[i + 1], _POWER_REFINEMENT + 1, dtype=xp.float64)[:-1]
        for i in range(len(natural) - 1)
    ]
    pieces.append(xp.asarray([natural[-1]], dtype=xp.float64))
    breaks = [float(v) for v in xp.concat(pieces)]
    return BSpline(domain=basis.domain, order=_POWER_ORDER, breaks=breaks)


@dataclass(frozen=True, eq=False, init=False)
class FData(PlotMixin):
    """A set of curves expressed as ``x_j(t) = sum_i coefs[i, j] phi_i(t)``.

    Replaces R's ``fd`` object and the functions that operate on it.

    Parameters
    ----------
    coefs : array_like
        Coefficients of shape ``(n_basis,)``, ``(n_basis, n_curves)`` or
        ``(n_basis, n_curves, n_vars)``.  A one-dimensional array is read as a
        single curve.
    basis : Basis
        The expansion basis.

    Attributes
    ----------
    coefs : array
        The coefficient array, always at least two-dimensional.
    basis : Basis
        The expansion basis.

    Raises
    ------
    ValueError
        If ``coefs`` has more than three axes, or its first axis does not match
        ``basis.n_basis``.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> fd = fb.FData(np.eye(5), fb.BSpline(domain=(0.0, 1.0), n_basis=5))
    >>> fd.n_curves
    5
    >>> fd(np.array([0.0, 1.0])).shape
    (2, 5)
    """

    coefs: Array
    basis: Basis

    def __init__(self, coefs: Any, basis: Basis) -> None:
        values = asarray(coefs)
        if len(values.shape) == 1:
            values = values[:, None]
        if len(values.shape) > 3:
            raise ValueError(f"coefs may have at most three axes, got {len(values.shape)}")
        if values.shape[0] != basis.n_basis:
            raise ValueError(
                f"coefs must have {basis.n_basis} rows to match the basis, got {values.shape[0]}"
            )
        object.__setattr__(self, "coefs", values)
        object.__setattr__(self, "basis", basis)

    # ----------------------------------------------------------------- shape

    @property
    def n_curves(self) -> int:
        """Number of curves held by the object.

        Returns
        -------
        int
            ``coefs.shape[1]``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.ones((4, 3)), fb.BSpline(n_basis=4)).n_curves
        3
        """
        return int(self.coefs.shape[1])

    @property
    def n_vars(self) -> int:
        """Number of variables per curve.

        Returns
        -------
        int
            ``1`` for two-dimensional coefficients, else ``coefs.shape[2]``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.ones((4, 3, 2)), fb.BSpline(n_basis=4)).n_vars
        2
        """
        return int(self.coefs.shape[2]) if len(self.coefs.shape) == 3 else 1

    @property
    def domain(self) -> tuple[float, float]:
        """Interval on which the curves are defined.

        Returns
        -------
        tuple of float
            ``self.basis.domain``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.ones((4, 1)), fb.BSpline(domain=(0.0, 2.0), n_basis=4)).domain
        (0.0, 2.0)
        """
        return self.basis.domain

    def __len__(self) -> int:
        """Return the number of curves."""
        return self.n_curves

    def __iter__(self) -> Iterator[FData]:
        """Iterate over the curves one at a time."""
        return (self[i] for i in range(self.n_curves))

    def __repr__(self) -> str:
        """Return a short description of the object."""
        return (
            f"FData(n_curves={self.n_curves}, n_vars={self.n_vars}, "
            f"basis={type(self.basis).__name__}(n_basis={self.basis.n_basis}))"
        )

    def __getitem__(self, index: Any) -> FData:
        """Select curves by position.

        Parameters
        ----------
        index : int, slice, sequence of int, array or boolean mask
            Curve selector applied to the second axis of ``coefs``.  Negative
            positions count from the end.  A boolean mask must have one entry
            per curve.

        Returns
        -------
        FData
            The selected curves, always with a curve axis.

        Raises
        ------
        IndexError
            If a position is out of range, or a boolean mask has the wrong
            length.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.eye(4), fb.BSpline(n_basis=4))
        >>> len(fd[1:3])
        2
        >>> len(fd[-1])
        1
        """
        xp = array_namespace(self.coefs)
        if isinstance(index, slice):
            return FData(self.coefs[:, index, ...], self.basis)
        positions = _curve_positions(index, self.n_curves)
        if positions is None:
            # A single curve: slice rather than take, so the curve axis is kept
            # without a copy.
            place = operator.index(index)
            place += self.n_curves if place < 0 else 0
            if not 0 <= place < self.n_curves:
                raise IndexError(f"curve {index} is out of range for {self.n_curves} curves")
            return FData(self.coefs[:, place : place + 1, ...], self.basis)
        chosen = xp.asarray(positions, dtype=xp.int64)
        return FData(xp.take(self.coefs, chosen, axis=1), self.basis)

    # ------------------------------------------------------------ evaluation

    def __call__(self, t: Any, deriv: int | LDO = 0) -> Array:
        """Evaluate the curves at ``t``.

        Parameters
        ----------
        t : array_like
            Evaluation points, shape ``(n_points,)``.
        deriv : int or LDO, optional
            Derivative order, or a linear differential operator.  Defaults to
            ``0``.

        Returns
        -------
        array
            Values of shape ``(n_points, n_curves)``, or
            ``(n_points, n_curves, n_vars)`` for multivariate coefficients, in
            the array namespace of ``t``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((4, 1)), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> float(fd(np.array([0.5]))[0, 0])
        1.0
        """
        xp = result_namespace(t, self.coefs)
        mat = self.basis(asarray(t, xp=xp), deriv)
        coefs = asarray(self.coefs, xp=xp)
        flat = xp.reshape(coefs, (coefs.shape[0], -1))
        values = xp.matmul(mat, flat)
        return xp.reshape(values, (mat.shape[0], *coefs.shape[1:]))

    def to_numpy(self, t: Any) -> NDArray:
        """Evaluate at ``t`` and return a NumPy array.

        Parameters
        ----------
        t : array_like
            Evaluation points.

        Returns
        -------
        numpy.ndarray
            The evaluated curves.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((4, 1)), fb.BSpline(n_basis=4))
        >>> fd.to_numpy(np.array([0.5])).shape
        (1, 1)
        """
        return to_numpy(self(to_numpy(t)))

    def to_torch(self, t: Any) -> Array:
        """Evaluate at ``t`` and return a PyTorch tensor.

        Parameters
        ----------
        t : array_like
            Evaluation points.

        Returns
        -------
        torch.Tensor
            The evaluated curves, differentiable with respect to ``t`` and the
            coefficients.

        Raises
        ------
        ImportError
            If PyTorch is not installed.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((4, 1)), fb.BSpline(n_basis=4))
        >>> fd.to_torch(np.array([0.5])).shape  # doctest: +SKIP
        torch.Size([1, 1])
        """
        import torch

        points = torch.as_tensor(to_numpy(t), dtype=torch.float64)
        coefs = torch.as_tensor(to_numpy(self.coefs), dtype=torch.float64)
        return FData(coefs, self.basis)(points)

    # ------------------------------------------------------------ derivative

    def derivative(self, n: int = 1) -> FData:
        """Return the ``n``-th derivative as a new :class:`FData`.

        The derivative is exact: each basis reports the space its derivative
        lives in together with the matching coefficient map, so no numerical
        projection is involved.

        Parameters
        ----------
        n : int, optional
            Derivative order.  Defaults to ``1``.

        Returns
        -------
        FData
            The derivative, generally in a different basis.

        Raises
        ------
        ValueError
            If ``n`` is negative.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.arange(6.0), fb.BSpline(domain=(0.0, 1.0), n_basis=6))
        >>> fd.derivative().basis.order
        3
        """
        if n < 0:
            raise ValueError(f"derivative order must be non-negative, got {n}")
        basis, matrix = self.basis._derivative_map(int(n))
        xp = array_namespace(self.coefs)
        derivative_map = asarray(matrix, xp=xp)
        flat = xp.reshape(self.coefs, (self.coefs.shape[0], -1))
        coefs = xp.matmul(derivative_map, flat)
        return FData(xp.reshape(coefs, (basis.n_basis, *self.coefs.shape[1:])), basis)

    # ------------------------------------------------------------ statistics

    def mean(self) -> FData:
        """Return the pointwise mean curve.

        Returns
        -------
        FData
            A single curve, the average of the coefficients.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.eye(4), fb.BSpline(n_basis=4)).mean().n_curves
        1
        """
        xp = array_namespace(self.coefs)
        return FData(xp.mean(self.coefs, axis=1, keepdims=True), self.basis)

    def center(self) -> FData:
        """Return the curves with the mean curve removed.

        Returns
        -------
        FData
            Centred curves in the same basis.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.eye(4), fb.BSpline(n_basis=4)).center()
        >>> bool(abs(float(fd.coefs.sum())) < 1e-12)
        True
        """
        return FData(self.coefs - self.mean().coefs, self.basis)

    def std(self) -> FData:
        """Return the pointwise sample standard deviation curve.

        The pointwise standard deviation is not a member of the basis span --
        it involves a square root -- so it is evaluated on a fine grid and
        projected back by least squares, as R's ``sd.fd`` does.

        Returns
        -------
        FData
            A single curve in the same basis, keeping this object's variables:
            coefficients of shape ``(n_basis, 1)``, or ``(n_basis, 1, n_vars)``
            when several variables are present.  The deviation is taken over
            the replications only, never across variables.

        Raises
        ------
        ValueError
            If fewer than two curves are present.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.eye(4), fb.BSpline(n_basis=4)).std().n_curves
        1
        """
        if self.n_curves < 2:
            raise ValueError("a standard deviation needs at least two curves")
        size = max(201, 10 * self.basis.n_basis)
        lower, upper = self.domain
        grid = default_namespace().linspace(lower, upper, size, dtype=default_namespace().float64)
        xp = array_namespace(self.coefs)
        values = self(asarray(grid, xp=xp))
        deviation = xp.std(values, axis=1, correction=1)
        target = deviation[:, None] if len(deviation.shape) == 1 else deviation
        fitted = _linalg.lstsq(asarray(self.basis(grid), xp=xp), target)
        if len(values.shape) == 3:
            fitted = xp.reshape(fitted, (self.basis.n_basis, 1, self.n_vars))
        return FData(fitted, self.basis)

    def cov(self) -> BiFData:
        """Return the sample covariance surface.

        Returns
        -------
        BiFData
            The bivariate function ``c(s, t)`` estimating
            ``Cov(x(s), x(t))`` with denominator ``n_curves - 1``.  With several
            variables each gets its own surface and they are stacked on a
            trailing axis, so the coefficients have shape
            ``(n_basis, n_basis, n_vars)`` and evaluation returns
            ``(len(s), len(t), n_vars)``.  Variables are never pooled as extra
            replications, and cross-variable covariances are not formed.

        Raises
        ------
        ValueError
            If fewer than two curves are present.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fb.FData(np.eye(4), fb.BSpline(n_basis=4)).cov().coefs.shape
        (4, 4)
        """
        if self.n_curves < 2:
            raise ValueError("a covariance needs at least two curves")
        xp = array_namespace(self.coefs)
        centred = self.center().coefs
        if len(centred.shape) == 2:
            coefs = xp.matmul(centred, xp.matrix_transpose(centred)) / (self.n_curves - 1)
            return BiFData(coefs, self.basis, self.basis)
        per_variable = xp.permute_dims(centred, (2, 0, 1))
        gram = xp.matmul(per_variable, xp.matrix_transpose(per_variable))
        coefs = xp.permute_dims(gram, (1, 2, 0)) / (self.n_curves - 1)
        return BiFData(coefs, self.basis, self.basis)

    # ------------------------------------------------------------ arithmetic

    def _combine(self, other: FData, sign: float) -> FData:
        """Add ``sign * other`` to ``self``, requiring a shared basis."""
        if self.basis != other.basis:
            raise ValueError("adding two FData objects requires the same basis")
        if self.n_curves != other.n_curves and 1 not in (self.n_curves, other.n_curves):
            raise ValueError(f"cannot combine {self.n_curves} curves with {other.n_curves} curves")
        xp = result_namespace(self.coefs, other.coefs)
        left = asarray(self.coefs, xp=xp)
        right = asarray(other.coefs, xp=xp)
        return FData(left + sign * right, self.basis)

    def _constant(self, value: Any) -> FData:
        """Return the constant function ``value`` expanded in this basis.

        Adding a scalar is only meaningful if the basis spans the constants.
        Most do -- splines and polygonals through the partition of unity,
        Fourier through its first harmonic, monomials and powers through the
        zero exponent, exponentials through the zero rate -- but a basis such as
        ``Monomial(exponents=[1, 2])`` does not, and the L2 projection would
        silently return a different function.

        Raises
        ------
        ValueError
            If the basis cannot reproduce the constant.
        """
        nodes, weights = _quadrature(self.basis)
        coefs, scalar = _promote(self.coefs, value)
        xp = array_namespace(coefs)
        target = xp.ones((nodes.shape[0], 1), dtype=coefs.dtype) * scalar
        fitted = _project(self.basis, nodes, weights, target)
        nxp = default_namespace()
        wanted = to_numpy(target)
        residual = nxp.matmul(self.basis(nodes), to_numpy(fitted)) - wanted
        size = float(nxp.max(nxp.abs(wanted)))
        if float(nxp.max(nxp.abs(residual))) > _CONSTANT_TOL * max(1.0, size):
            raise ValueError(
                f"{type(self.basis).__name__} cannot represent the constant {size:g}; "
                "add an FData on the same basis instead, or use a basis that spans "
                "the constants (a spline, a Fourier basis, or a zero exponent or rate)"
            )
        return FData(fitted, self.basis)

    def __add__(self, other: FData | float) -> FData:
        """Add another function on the same basis, or a scalar."""
        if isinstance(other, FData):
            return self._combine(other, 1.0)
        return self._combine(self._constant(other), 1.0)

    __radd__ = __add__

    def __neg__(self) -> FData:
        """Return the pointwise negation."""
        return FData(-self.coefs, self.basis)

    def __sub__(self, other: FData | float) -> FData:
        """Subtract another function on the same basis, or a scalar."""
        if isinstance(other, FData):
            return self._combine(other, -1.0)
        return self._combine(self._constant(other), -1.0)

    def __rsub__(self, other: float) -> FData:
        """Subtract this function from a scalar."""
        return (-self) + other

    def __mul__(self, other: FData | float) -> FData:
        """Multiply by a scalar, or form the pointwise product of two functions.

        The product of two functions is expanded in ``basis1 * basis2``, which
        for splines carries the knot multiplicities that make the product exact.
        """
        if not isinstance(other, FData):
            coefs, scalar = _promote(self.coefs, other)
            return FData(coefs * scalar, self.basis)
        if not _same_domain(self.basis.domain, other.basis.domain):
            raise ValueError("multiplying two FData objects requires the same domain")
        if self.n_curves != other.n_curves and 1 not in (self.n_curves, other.n_curves):
            raise ValueError(f"cannot combine {self.n_curves} curves with {other.n_curves} curves")
        product = self.basis * other.basis
        nodes, weights = _quadrature(self.basis, other.basis, product)
        xp = result_namespace(self.coefs, other.coefs)
        values = self(asarray(nodes, xp=xp)) * other(asarray(nodes, xp=xp))
        return FData(_project(product, nodes, weights, values), product)

    __rmul__ = __mul__

    def __truediv__(self, other: float | FData) -> FData:
        """Divide by a scalar.

        Raises
        ------
        TypeError
            If ``other`` is an :class:`FData`: the ratio of two basis
            expansions is not closed in any of the supported bases.
        """
        if isinstance(other, FData):
            raise TypeError("division of two FData objects is not closed in a basis expansion")
        coefs, scalar = _promote(self.coefs, other)
        return FData(coefs / scalar, self.basis)

    def __pow__(self, power: float) -> FData:
        """Raise the curves to a power.

        Non-negative integer powers are exact repeated products.  Any other
        exponent is projected onto a refined spline basis, since the result is
        not a basis expansion of the original space.

        Parameters
        ----------
        power : int or float
            The exponent.

        Returns
        -------
        FData
            The powered curves.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.ones((4, 1)), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> round(float((fd**2)(np.array([0.5]))[0, 0]), 10)
        1.0
        """
        if isinstance(power, int) or float(power).is_integer():
            exponent = int(power)
            if exponent == 0:
                xp = array_namespace(self.coefs)
                shape = (1, *self.coefs.shape[1:])
                return FData(xp.ones(shape, dtype=self.coefs.dtype), Constant(self.domain))
            if exponent > 0:
                out = self
                for _ in range(exponent - 1):
                    out = out * self
                return out
        target = _refined_spline(self.basis)
        nodes, weights = _quadrature(target)
        xp = array_namespace(self.coefs)
        values = self(asarray(nodes, xp=xp)) ** float(power)
        return FData(_project(target, nodes, weights, values), target)

    def __matmul__(self, other: FData) -> Array:
        """Return the matrix of inner products with ``other``.

        Parameters
        ----------
        other : FData
            Right-hand functions.

        Returns
        -------
        array
            Matrix of shape ``(self.n_curves, other.n_curves)``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.eye(4), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> (fd @ fd).shape
        (4, 4)
        """
        return inprod(self, other)

    def to_pandas(self, t: Any) -> Any:
        """Evaluate on ``t`` and return a long-format :class:`pandas.DataFrame`.

        Thin delegation to :func:`fabel.io.to_pandas`; requires the ``pandas``
        extra.

        Parameters
        ----------
        t : array_like
            Argument values to evaluate at.

        Returns
        -------
        pandas.DataFrame
            One row per (t, curve[, var]) combination.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.eye(4), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> list(fd.to_pandas(np.array([0.0, 1.0])).columns)
        ['t', 'curve', 'value']
        """
        from fabel.io import to_pandas as _to_pandas

        return _to_pandas(self, t)

    def to_xarray(self, t: Any) -> Any:
        """Evaluate on ``t`` and return an :class:`xarray.DataArray`.

        Thin delegation to :func:`fabel.io.to_xarray`; requires the ``pandas``
        extra (which pulls in ``xarray``).

        Parameters
        ----------
        t : array_like
            Argument values to evaluate at.

        Returns
        -------
        xarray.DataArray
            Dims ``("t", "curve")``, or ``("t", "curve", "var")`` when
            ``self.n_vars > 1``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> fd = fb.FData(np.eye(4), fb.BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> fd.to_xarray(np.array([0.0, 1.0])).dims
        ('t', 'curve')
        """
        from fabel.io import to_xarray as _to_xarray

        return _to_xarray(self, t)


def inprod(
    first: FData | Basis,
    second: FData | Basis,
    lfd1: int | LDO = 0,
    lfd2: int | LDO = 0,
) -> Array:
    """Return the matrix of inner products ``integral (lfd1 x_i)(lfd2 y_j) dt``.

    Replaces R's ``inprod`` and ``inprod.bspline``.  A :class:`~fabel.basis.Basis`
    may be passed in place of an :class:`FData`, in which case its basis
    functions play the role of the curves.

    Parameters
    ----------
    first, second : FData or Basis
        Left- and right-hand functions.  They must share a domain.
    lfd1, lfd2 : int or LDO, optional
        Operators applied to the left and right arguments.  Default ``0``.

    Returns
    -------
    array
        Matrix of shape ``(n_left, n_right)``.

    Raises
    ------
    ValueError
        If the two arguments live on different domains, or either carries a
        variable axis.

    Examples
    --------
    >>> import fabel as fb
    >>> b = fb.BSpline(domain=(0.0, 1.0), n_basis=5)
    >>> fb.inprod(b, b).shape
    (5, 5)
    """
    left = first.basis if isinstance(first, FData) else first
    right = second.basis if isinstance(second, FData) else second
    if not _same_domain(left.domain, right.domain):
        raise ValueError(f"inner products need one domain, got {left.domain} and {right.domain}")
    for side in (first, second):
        if isinstance(side, FData) and len(side.coefs.shape) != 2:
            raise ValueError("inner products are defined for coefficients without a variable axis")
    # The cross-Gram matrix depends only on the two bases, so it is built in
    # NumPy once and then converted; the contractions with the coefficients run
    # in the coefficients' namespace so that gradients survive.
    sides = [side.coefs for side in (first, second) if isinstance(side, FData)]
    xp = result_namespace(*sides) if sides else default_namespace()
    matrix = asarray(_cross_gram(left, right, _as_operator(lfd1), _as_operator(lfd2)), xp=xp)
    if isinstance(first, FData):
        matrix = xp.matmul(xp.matrix_transpose(asarray(first.coefs, xp=xp)), matrix)
    if isinstance(second, FData):
        matrix = xp.matmul(matrix, asarray(second.coefs, xp=xp))
    return matrix


@dataclass(frozen=True, eq=False, init=False)
class BiFData:
    """A bivariate function ``c(s, t) = sum_ij coefs[i, j] phi_i(s) psi_j(t)``.

    Replaces R's ``bifd`` object.

    Parameters
    ----------
    coefs : array_like
        Coefficients whose first two axes match ``sbasis`` and ``tbasis``.
        Trailing axes (replications, variables) are carried through untouched.
    sbasis : Basis
        Basis for the first argument.
    tbasis : Basis
        Basis for the second argument.

    Raises
    ------
    ValueError
        If the first two axes of ``coefs`` do not match the two bases.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> b = fb.BSpline(domain=(0.0, 1.0), n_basis=4)
    >>> bifd = fb.BiFData(np.eye(4), b, b)
    >>> bifd(np.array([0.0, 1.0]), np.array([0.0, 0.5, 1.0])).shape
    (2, 3)
    """

    coefs: Array
    sbasis: Basis
    tbasis: Basis

    def __init__(self, coefs: Any, sbasis: Basis, tbasis: Basis) -> None:
        values = asarray(coefs)
        expected = (sbasis.n_basis, tbasis.n_basis)
        if len(values.shape) < 2 or tuple(values.shape[:2]) != expected:
            raise ValueError(
                f"coefs must have shape {expected} + trailing axes, got {values.shape}"
            )
        object.__setattr__(self, "coefs", values)
        object.__setattr__(self, "sbasis", sbasis)
        object.__setattr__(self, "tbasis", tbasis)

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:
        """The two argument domains.

        Returns
        -------
        tuple
            ``(sbasis.domain, tbasis.domain)``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> b = fb.BSpline(domain=(0.0, 1.0), n_basis=4)
        >>> fb.BiFData(np.eye(4), b, b).domain
        ((0.0, 1.0), (0.0, 1.0))
        """
        return self.sbasis.domain, self.tbasis.domain

    def __repr__(self) -> str:
        """Return a short description of the object."""
        return (
            f"BiFData(coefs={tuple(int(v) for v in self.coefs.shape)}, "
            f"sbasis={type(self.sbasis).__name__}, tbasis={type(self.tbasis).__name__})"
        )

    def __call__(self, s: Any, t: Any, deriv: tuple[int, int] = (0, 0)) -> Array:
        """Evaluate the surface on the grid ``s`` x ``t``.

        Parameters
        ----------
        s, t : array_like
            Evaluation points for the first and second argument.
        deriv : tuple of int, optional
            Derivative orders applied to each argument.  Defaults to ``(0, 0)``.

        Returns
        -------
        array
            Values of shape ``(len(s), len(t))`` followed by any trailing axes
            of ``coefs``.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> b = fb.BSpline(domain=(0.0, 1.0), n_basis=4)
        >>> fb.BiFData(np.eye(4), b, b)(np.array([0.5]), np.array([0.5])).shape
        (1, 1)
        """
        xp = result_namespace(s, t, self.coefs)
        smat = self.sbasis(asarray(s, xp=xp), deriv[0])
        tmat = self.tbasis(asarray(t, xp=xp), deriv[1])
        coefs = asarray(self.coefs, xp=xp)
        trailing = tuple(coefs.shape[2:])
        flat = xp.reshape(coefs, (coefs.shape[0], coefs.shape[1], -1))
        # sum_ij phi_i(s) c_ijk psi_j(t) -- contract the s axis, then the t axis.
        step = xp.matmul(smat, xp.reshape(flat, (flat.shape[0], -1)))
        step = xp.reshape(step, (smat.shape[0], flat.shape[1], flat.shape[2]))
        step = xp.permute_dims(step, (0, 2, 1))  # (n_s, n_trailing, n_tbasis)
        out = xp.matmul(step, xp.matrix_transpose(tmat))  # (n_s, n_trailing, n_t)
        out = xp.permute_dims(out, (0, 2, 1))
        return xp.reshape(out, (smat.shape[0], tmat.shape[0], *trailing))

    def transpose(self) -> BiFData:
        """Return the surface with its two arguments swapped.

        Returns
        -------
        BiFData
            The transposed surface.

        Examples
        --------
        >>> import numpy as np
        >>> import fabel as fb
        >>> b = fb.BSpline(domain=(0.0, 1.0), n_basis=4)
        >>> fb.BiFData(np.eye(4), b, b).transpose().coefs.shape
        (4, 4)
        """
        xp = array_namespace(self.coefs)
        axes = (1, 0, *range(2, len(self.coefs.shape)))
        return BiFData(xp.permute_dims(self.coefs, axes), self.tbasis, self.sbasis)
