"""Linear algebra helpers that have no array-API equivalent.

Together with :mod:`fabel._backend` this is the only module allowed to import
NumPy and SciPy directly.  Every routine dispatches on the concrete backend:
NumPy inputs use LAPACK through SciPy (banded Cholesky where the matrix is
banded), PyTorch inputs use :mod:`torch.linalg` so that gradients keep flowing.

Examples
--------
>>> import numpy as np
>>> from fabel._linalg import solve_spd
>>> a = np.array([[4.0, 1.0], [1.0, 3.0]])
>>> solve_spd(a, np.array([1.0, 2.0])).round(6).tolist()
[0.090909, 0.636364]
"""

from __future__ import annotations

from collections.abc import Callable, Hashable
from itertools import pairwise
from typing import Any

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from numpy.polynomial.legendre import leggauss

from fabel._backend import array_namespace, asarray, is_torch, to_numpy

__all__ = [
    "bandwidth_of",
    "cached_gram",
    "clear_gram_cache",
    "composite_gauss_legendre",
    "gauss_legendre",
    "gauss_legendre_reference",
    "lstsq",
    "pencil_eigh",
    "solve_spd",
    "sparse_penalty",
    "to_banded",
]

NDArray = np.ndarray[Any, np.dtype[np.float64]]

_BAND_TOL = 1e-14
_MAX_BANDED_BANDWIDTH_RATIO = 0.5


# --------------------------------------------------------------------------- #
# banded storage
# --------------------------------------------------------------------------- #


def bandwidth_of(a: NDArray, tol: float = _BAND_TOL) -> int:
    """Return the number of non-zero super-diagonals of ``a``.

    Parameters
    ----------
    a : numpy.ndarray
        Square matrix.
    tol : float, optional
        Entries with magnitude at or below ``tol`` count as zero.

    Returns
    -------
    int
        Largest ``k`` such that some ``|a[i, i + k]| > tol``; ``0`` for a
        diagonal matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import bandwidth_of
    >>> bandwidth_of(np.eye(4))
    0
    """
    n = a.shape[0]
    scale = float(np.max(np.abs(a))) if a.size else 0.0
    cutoff = max(tol, tol * scale)
    for k in range(n - 1, 0, -1):
        if float(np.max(np.abs(np.diag(a, k)))) > cutoff:
            return k
    return 0


def to_banded(a: NDArray, bandwidth: int) -> NDArray:
    """Pack the upper triangle of a symmetric banded matrix into LAPACK form.

    The result ``ab`` satisfies ``ab[bandwidth + i - j, j] == a[i, j]`` for
    ``i <= j``, which is the layout expected by
    :func:`scipy.linalg.cholesky_banded` with ``lower=False``.

    Parameters
    ----------
    a : numpy.ndarray
        Symmetric ``(n, n)`` matrix.
    bandwidth : int
        Number of super-diagonals to keep.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(bandwidth + 1, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import to_banded
    >>> to_banded(np.eye(3), 0).tolist()
    [[1.0, 1.0, 1.0]]
    """
    n = a.shape[0]
    ab = np.zeros((bandwidth + 1, n), dtype=np.float64)
    for k in range(bandwidth + 1):
        ab[bandwidth - k, k:] = np.diag(a, k)
    return ab


# --------------------------------------------------------------------------- #
# solvers
# --------------------------------------------------------------------------- #


def solve_spd(a: Any, b: Any, bandwidth: int | None = None) -> Any:
    """Solve ``a @ x = b`` for a symmetric, ideally positive-definite ``a``.

    For NumPy inputs a banded Cholesky factorisation is used when ``a`` is
    narrowly banded (the B-spline Gram and penalty matrices are), a dense
    Cholesky otherwise, and a symmetric LDL solve if ``a`` turns out not to be
    positive definite.  For PyTorch inputs :func:`torch.linalg.cholesky` is used
    with a :func:`torch.linalg.solve` fallback, so gradients flow.

    Parameters
    ----------
    a : array
        Symmetric ``(n, n)`` matrix.
    b : array
        Right-hand side of shape ``(n,)`` or ``(n, k)``.
    bandwidth : int, optional
        Known number of super-diagonals.  Detected automatically when omitted.

    Returns
    -------
    array
        Solution ``x`` with the shape of ``b``, in the namespace of the inputs.

    Raises
    ------
    ValueError
        If the shapes of ``a`` and ``b`` are incompatible.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import solve_spd
    >>> solve_spd(np.array([[2.0, 0.0], [0.0, 4.0]]), np.array([2.0, 4.0])).round(9).tolist()
    [1.0, 1.0]
    """
    xp = array_namespace(a, b)
    if a.shape[0] != a.shape[1] or a.shape[0] != b.shape[0]:
        raise ValueError(f"incompatible shape: a is {tuple(a.shape)}, b is {tuple(b.shape)}")
    if is_torch(a) or is_torch(b):
        return _solve_spd_torch(xp, a, b)
    return _solve_spd_numpy(
        np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64), bandwidth
    )


def _solve_spd_torch(xp: Any, a: Any, b: Any) -> Any:
    vector_rhs = b.ndim == 1
    rhs = b[:, None] if vector_rhs else b
    try:
        factor = xp.linalg.cholesky(a)
        out = xp.linalg.solve_triangular(
            xp.matrix_transpose(factor),
            xp.linalg.solve_triangular(factor, rhs, upper=False),
            upper=True,
        )
    except RuntimeError:
        out = xp.linalg.solve(a, rhs)
    return out[:, 0] if vector_rhs else out


def _solve_spd_numpy(a: NDArray, b: NDArray, bandwidth: int | None) -> NDArray:
    vector_rhs = b.ndim == 1
    rhs = b[:, None] if vector_rhs else b
    n = a.shape[0]
    band = bandwidth_of(a) if bandwidth is None else bandwidth
    out: NDArray
    if band < _MAX_BANDED_BANDWIDTH_RATIO * n:
        try:
            factor = sla.cholesky_banded(to_banded(a, band), lower=False)
            out = np.asarray(sla.cho_solve_banded((factor, False), rhs))
            return out[:, 0] if vector_rhs else out
        except sla.LinAlgError:
            pass
    try:
        out = np.asarray(sla.cho_solve(sla.cho_factor(a), rhs))
    except sla.LinAlgError:
        out = np.asarray(sla.solve(a, rhs, assume_a="sym"))
    return out[:, 0] if vector_rhs else out


def lstsq(a: Any, b: Any, ridge: float = 0.0) -> Any:
    """Least-squares solution of ``a @ x ~= b`` via the normal equations.

    The normal equations are used rather than a QR/SVD driver because they are
    differentiable on every backend and because the design matrices in Fabel are
    tall and very well conditioned (basis matrices on a fine grid).

    Parameters
    ----------
    a : array
        Design matrix of shape ``(m, n)`` with ``m >= n``.
    b : array
        Right-hand side of shape ``(m,)`` or ``(m, k)`` (further trailing axes
        are allowed and are flattened internally).
    ridge : float, optional
        Tikhonov weight added to the diagonal of ``aᵀa``.  Defaults to ``0``.

    Returns
    -------
    array
        Coefficients of shape ``(n,)`` or ``(n, k)``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import lstsq
    >>> a = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    >>> lstsq(a, np.array([1.0, 1.0, 2.0])).round(10).tolist()
    [1.0, 1.0]
    """
    xp = array_namespace(a, b)
    trailing = b.shape[1:]
    rhs = xp.reshape(b, (b.shape[0], -1)) if len(trailing) != 1 else b
    gram = xp.matmul(xp.matrix_transpose(a), a)
    if ridge:
        gram = gram + ridge * xp.eye(gram.shape[0], dtype=gram.dtype)
    coef = solve_spd(gram, xp.matmul(xp.matrix_transpose(a), rhs))
    if len(trailing) != 1:
        return xp.reshape(coef, (a.shape[1], *trailing))
    return coef


def pencil_eigh(a: NDArray, b: NDArray) -> tuple[NDArray, NDArray]:
    r"""Diagonalise the symmetric pencil ``a v = mu * b v`` with ``vᵀ b v = I``.

    A smoothing problem evaluates ``(S + λ R)⁻¹`` for many values of ``λ``.
    Taking ``b = S + R`` once gives ``vᵀ S v = diag(mu)`` and
    ``vᵀ R v = I - diag(mu)``, hence

    .. math:: (S + \lambda R)^{-1} = V\,\mathrm{diag}\!\big((\mu + \lambda(1-\mu))^{-1}\big)\,V^{T}

    so every further ``λ`` costs one diagonal scaling instead of a fresh
    factorisation.

    Parameters
    ----------
    a : numpy.ndarray
        Symmetric ``(n, n)`` matrix.
    b : numpy.ndarray
        Symmetric ``(n, n)`` matrix, positive definite.  When ``b`` is
        numerically singular a relative ridge of ``1e-12`` times its mean
        diagonal is added; the affected directions are the ones constrained
        neither by the data nor by the penalty, so the fit is unchanged.

    Returns
    -------
    mu : numpy.ndarray
        Eigenvalues, ascending, of shape ``(n,)``.
    v : numpy.ndarray
        Eigenvectors as columns, ``b``-orthonormal, of shape ``(n, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import pencil_eigh
    >>> s = np.array([[2.0, 0.0], [0.0, 1.0]])
    >>> mu, v = pencil_eigh(s, s + np.eye(2))
    >>> np.round(mu, 6).tolist()
    [0.5, 0.666667]
    """
    left = np.asarray(a, dtype=np.float64)
    right = np.asarray(b, dtype=np.float64)
    try:
        mu, vec = sla.eigh(left, right)
    except sla.LinAlgError:
        ridge = 1e-12 * float(np.mean(np.diag(right)))
        mu, vec = sla.eigh(left, right + ridge * np.eye(right.shape[0]))
    return np.asarray(mu, dtype=np.float64), np.asarray(vec, dtype=np.float64)


# --------------------------------------------------------------------------- #
# sparse assembly
# --------------------------------------------------------------------------- #


def sparse_penalty(a: NDArray, bandwidth: int | None = None) -> sp.dia_array:
    """Assemble a symmetric banded penalty matrix in SciPy DIA format.

    Parameters
    ----------
    a : numpy.ndarray
        Dense symmetric ``(n, n)`` penalty matrix.
    bandwidth : int, optional
        Number of super-diagonals to retain; detected automatically when
        omitted.  Entries outside the band are dropped.

    Returns
    -------
    scipy.sparse.dia_array
        Sparse representation with offsets ``-bandwidth .. bandwidth``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import sparse_penalty
    >>> sparse_penalty(np.eye(3)).toarray().tolist()
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    """
    dense = np.asarray(a, dtype=np.float64)
    n = dense.shape[0]
    band = bandwidth_of(dense) if bandwidth is None else bandwidth
    offsets = np.arange(-band, band + 1)
    data = np.zeros((offsets.size, n), dtype=np.float64)
    for row, k in enumerate(offsets):
        diag = np.diag(dense, int(k))
        if k >= 0:
            data[row, k:] = diag
        else:
            data[row, : n + int(k)] = diag
    return sp.dia_array((data, offsets), shape=(n, n))


# --------------------------------------------------------------------------- #
# quadrature
# --------------------------------------------------------------------------- #


def gauss_legendre_reference(deg: int) -> tuple[NDArray, NDArray]:
    """Return cached Gauss-Legendre nodes and weights on ``[-1, 1]``.

    Parameters
    ----------
    deg : int
        Number of nodes.  The rule is exact for polynomials of degree
        ``2 * deg - 1``.

    Returns
    -------
    nodes, weights : numpy.ndarray
        Read-only arrays of length ``deg``.

    Examples
    --------
    >>> from fabel._linalg import gauss_legendre_reference
    >>> nodes, weights = gauss_legendre_reference(2)
    >>> weights.tolist()
    [1.0, 1.0]
    """
    return _gauss_legendre_reference_cached(int(deg))


def _gauss_legendre_reference_cached(deg: int) -> tuple[NDArray, NDArray]:
    cached = _GL_CACHE.get(deg)
    if cached is None:
        nodes, weights = leggauss(deg)
        nodes = np.ascontiguousarray(nodes, dtype=np.float64)
        weights = np.ascontiguousarray(weights, dtype=np.float64)
        nodes.flags.writeable = False
        weights.flags.writeable = False
        cached = (nodes, weights)
        _GL_CACHE[deg] = cached
    return cached


_GL_CACHE: dict[int, tuple[NDArray, NDArray]] = {}


def gauss_legendre(deg: int, lower: float, upper: float) -> tuple[NDArray, NDArray]:
    """Gauss-Legendre nodes and weights mapped onto ``[lower, upper]``.

    Parameters
    ----------
    deg : int
        Number of nodes; exact for polynomials of degree ``2 * deg - 1``.
    lower, upper : float
        Interval endpoints.

    Returns
    -------
    nodes, weights : numpy.ndarray
        Arrays of length ``deg``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import gauss_legendre
    >>> nodes, weights = gauss_legendre(3, 0.0, 1.0)
    >>> float(np.sum(weights))
    1.0
    """
    ref_nodes, ref_weights = gauss_legendre_reference(deg)
    half = 0.5 * (upper - lower)
    mid = 0.5 * (upper + lower)
    return half * ref_nodes + mid, half * ref_weights


def composite_gauss_legendre(breaks: NDArray, deg: int) -> tuple[NDArray, NDArray]:
    """Composite Gauss-Legendre rule over the panels defined by ``breaks``.

    Zero-width panels (repeated knots) are skipped.

    Parameters
    ----------
    breaks : numpy.ndarray
        Non-decreasing panel boundaries, length at least 2.
    deg : int
        Number of nodes per panel.

    Returns
    -------
    nodes, weights : numpy.ndarray
        Concatenated nodes and weights over all non-degenerate panels.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import composite_gauss_legendre
    >>> nodes, weights = composite_gauss_legendre(np.array([0.0, 1.0, 2.0]), 2)
    >>> float(np.sum(weights))
    2.0
    """
    edges = np.asarray(breaks, dtype=np.float64)
    node_blocks: list[NDArray] = []
    weight_blocks: list[NDArray] = []
    for lower, upper in pairwise(edges):
        if upper <= lower:
            continue
        nodes, weights = gauss_legendre(deg, float(lower), float(upper))
        node_blocks.append(nodes)
        weight_blocks.append(weights)
    if not node_blocks:
        return np.zeros(0), np.zeros(0)
    return np.concatenate(node_blocks), np.concatenate(weight_blocks)


# --------------------------------------------------------------------------- #
# Gram-matrix cache
# --------------------------------------------------------------------------- #

_GRAM_CACHE: dict[tuple[Hashable, int], NDArray] = {}
_GRAM_CACHE_MAXSIZE = 256


def cached_gram(key: Hashable, deriv: int, compute: Callable[[], NDArray]) -> NDArray:
    """Memoise a basis Gram matrix on ``(basis key, derivative order)``.

    The cached array is made read-only so that callers cannot corrupt the cache.

    Parameters
    ----------
    key : hashable
        Identity of the basis (its frozen dataclass fields).
    deriv : int
        Derivative order the Gram matrix was computed for.
    compute : callable
        Zero-argument function producing the matrix on a cache miss.

    Returns
    -------
    numpy.ndarray
        The cached, read-only Gram matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import cached_gram, clear_gram_cache
    >>> clear_gram_cache()
    >>> cached_gram(("demo",), 0, lambda: np.eye(2)).tolist()
    [[1.0, 0.0], [0.0, 1.0]]
    """
    cache_key = (key, int(deriv))
    hit = _GRAM_CACHE.get(cache_key)
    if hit is not None:
        return hit
    value = np.array(compute(), dtype=np.float64, copy=True)
    value.flags.writeable = False
    if len(_GRAM_CACHE) >= _GRAM_CACHE_MAXSIZE:
        _GRAM_CACHE.pop(next(iter(_GRAM_CACHE)))
    _GRAM_CACHE[cache_key] = value
    return value


def clear_gram_cache() -> None:
    """Empty the Gram-matrix cache.

    Examples
    --------
    >>> from fabel._linalg import clear_gram_cache
    >>> clear_gram_cache()
    """
    _GRAM_CACHE.clear()


def as_backend(values: NDArray, like: Any) -> Any:
    """Move a NumPy constant into the namespace and dtype of ``like``.

    Parameters
    ----------
    values : numpy.ndarray
        Backend-independent constant, e.g. a Gram matrix.
    like : array or module
        Array whose namespace and dtype should be matched, or a namespace.

    Returns
    -------
    array
        ``values`` as an array of the target namespace.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel._linalg import as_backend
    >>> as_backend(np.eye(2), np.zeros(1)).shape
    (2, 2)
    """
    xp = like if hasattr(like, "asarray") else array_namespace(like)
    return asarray(to_numpy(values), xp=xp)
