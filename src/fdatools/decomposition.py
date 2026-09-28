r"""Functional principal components and canonical correlation.

Two estimators replace R ``fda``'s ``pca.fd``, ``varmx.pca.fd`` and ``cca.fd``.
Both work on the basis coefficients of an :class:`~fdatools.core.FData` and both
follow the scikit-learn estimator API, so they drop straight into a
:class:`~sklearn.pipeline.Pipeline` behind :class:`~fdatools.smoothing.Smoother`.

Functional PCA with a roughness penalty (Ramsay & Silverman §9.4) maximises the
sample variance of the scores :math:`\int \xi(t) x_i(t)\,dt` while charging the
harmonic :math:`\xi` for its own roughness.  Writing :math:`W = \int
\varphi\varphi^{T}` for the Gram matrix of the basis, :math:`R = \int
(L\varphi)(L\varphi)^{T}` for the roughness penalty and :math:`V = C C^{T}/N`
for the sample covariance of the centred coefficients, the harmonic
coefficients solve the symmetric generalised eigenproblem

.. math::

    W V W\,b = \mu\,(W + \lambda R)\,b ,
    \qquad b^{T}(W + \lambda R)\,b = 1 .

With :math:`\lambda = 0` this collapses to the unpenalised problem
:math:`V W b = \mu b`, i.e. the eigen-decomposition of the covariance operator.

Functional CCA (chapter 11) is the same construction on the cross-covariance of
two samples: with :math:`S_{xy} = W_1 C_x C_y^{T} W_2 / N` the canonical weights
solve

.. math::

    \begin{pmatrix} 0 & S_{xy} \\ S_{yx} & 0 \end{pmatrix} v
    = \rho \begin{pmatrix} S_{xx} + \lambda_1 R_1 & 0 \\
                           0 & S_{yy} + \lambda_2 R_2 \end{pmatrix} v .

Every Gram matrix is exact (Gauss-Legendre on the basis breaks), where R's
``pca.fd`` uses a Romberg approximation good to four or five digits; see
``tests/parity/test_decomposition.py`` for the measured consequences.

Examples
--------
>>> import numpy as np
>>> from fdatools import BSpline, FData
>>> from fdatools.decomposition import FPCA
>>> rng = np.random.default_rng(0)
>>> basis = BSpline(domain=(0.0, 1.0), n_basis=8)
>>> fd = FData(rng.standard_normal((8, 30)), basis)
>>> pca = FPCA(n=3).fit(fd)
>>> pca.scores.shape
(30, 3)
>>> bool(np.all(np.diff(pca.values) <= 1e-12))
True
"""

from __future__ import annotations

import operator
from math import atan2, cos, sin, sqrt
from types import ModuleType
from typing import Any, SupportsIndex

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted, validate_data

from fdatools import _linalg
from fdatools._backend import asarray, default_namespace, to_numpy
from fdatools._operator import LDO
from fdatools.basis import Basis, BSpline
from fdatools.core import FData

__all__ = ["FCCA", "FPCA"]

Array = Any

#: Number of grid points R's ``varmx.pca.fd`` rotates on, and the value fdatools
#: uses so that a rotation is comparable with R's.
_ROTATION_GRID = 501

#: Angle below which the varimax sweep stops rotating a pair of columns.
_ROTATION_TOL = 1e-12

#: Maximum number of varimax sweeps.
_ROTATION_MAX_SWEEPS = 1000

#: Log-spaced grid searched when ``lam="gcv"``: ``10 ** -8, ..., 10 ** 8``.
_LAMBDA_GRID = tuple(10.0**exponent for exponent in range(-8, 9))

#: Default B-spline order used when an estimator has to invent a basis.
_DEFAULT_ORDER = 4


def _xp() -> ModuleType:
    """Return the NumPy array namespace every estimator here computes in."""
    return default_namespace()


def _coefficients(fd: FData) -> Array:
    """Return the coefficients of ``fd`` as a float64 array with a curve axis."""
    xp = _xp()
    coefs = asarray(to_numpy(fd.coefs), xp)
    return coefs[:, None] if coefs.ndim == 1 else coefs


def _positive_int(value: object) -> int:
    """Return ``value`` as an ``int`` if it is a positive integer, else raise.

    Any object with ``__index__`` counts (``int``, ``numpy.int64``, ...), so
    the value is normalised with :func:`operator.index`.  Booleans (Python and
    NumPy, whose ``__index__`` raises) and floats -- even integral ones such
    as ``2.0`` -- are rejected.
    """
    # NumPy bools (type name ``bool`` / ``bool_``) are caught by name: NumPy < 2.3
    # warns instead of raising when one is passed to ``operator.index``.
    if isinstance(value, bool) or type(value).__name__ in {"bool", "bool_"}:
        raise ValueError(f"n must be a positive integer, got {value!r}")
    try:
        index = operator.index(value)  # type: ignore[arg-type]
    except TypeError:
        raise ValueError(f"n must be a positive integer, got {value!r}") from None
    if index < 1:
        raise ValueError(f"n must be a positive integer, got {value!r}")
    return index


def _default_basis(n_basis: int) -> Basis:
    """Return the fallback basis for a coefficient matrix with ``n_basis`` columns."""
    return BSpline(domain=(0.0, 1.0), n_basis=n_basis, order=min(_DEFAULT_ORDER, n_basis))


def _block_diagonal(matrix: Array, repeats: int) -> Array:
    """Return ``matrix`` repeated ``repeats`` times down a block diagonal."""
    if repeats == 1:
        return matrix
    xp = _xp()
    zero = xp.zeros_like(matrix)
    rows = [
        xp.concat([matrix if col == row else zero for col in range(repeats)], axis=1)
        for row in range(repeats)
    ]
    return xp.concat(rows, axis=0)


def _n_vars(coefs: Array) -> int:
    """Return the number of variables of a ``(n_basis, n_curves[, n_vars])`` array."""
    return 1 if coefs.ndim == 2 else int(coefs.shape[2])


def _stack_variables(coefs: Array) -> Array:
    """Flatten a ``(n_basis, n_curves[, n_vars])`` array to ``(n_basis * n_vars, n_curves)``."""
    if coefs.ndim == 2:
        return coefs
    return _xp().concat([coefs[:, :, k] for k in range(coefs.shape[2])], axis=0)


def _unstack_variables(stacked: Array, n_basis: int, n_vars: int) -> Array:
    """Undo :func:`_stack_variables` for a matrix of column vectors."""
    if n_vars == 1:
        return stacked
    parts = [stacked[k * n_basis : (k + 1) * n_basis, :] for k in range(n_vars)]
    return _xp().stack(parts, axis=2)


def _positive_sum_signs(vectors: Array) -> Array:
    """Return the ``±1`` per column that makes each column's coefficient sum positive.

    Eigenvector signs are arbitrary.  R's ``pca.fd`` reports every harmonic with
    a positive coefficient sum (measured on all 18 harmonics of the golden
    cases), so fdatools adopts that rule; a zero sum keeps the sign it has.
    """
    xp = _xp()
    totals = xp.sum(vectors, axis=0)
    return xp.where(totals < 0.0, -1.0, 1.0)


def _varimax_rotation(values: Array) -> Array:
    """Return the varimax rotation of the columns of ``values``.

    Uses the classical pairwise-Jacobi sweep on the *raw* (un-normalised)
    loadings, maximising ``sum_j [ sum_i a_ij^4 - (sum_i a_ij^2)^2 / n ]``.
    Each pair is rotated by the angle that maximises the criterion exactly, so
    the sweep stops at a stationary point of the criterion on the rotation
    group and the returned matrix is orthogonal to rounding.

    Parameters
    ----------
    values : array
        Loadings of shape ``(n_points, n_components)``.

    Returns
    -------
    array
        Orthogonal ``(n_components, n_components)`` rotation ``T`` such that
        ``values @ T`` maximises the varimax criterion.
    """
    xp = _xp()
    n_points, n_comp = int(values.shape[0]), int(values.shape[1])
    rotation = xp.eye(n_comp, dtype=xp.float64)
    loadings = xp.asarray(values, dtype=xp.float64, copy=True)
    if n_comp < 2:
        return rotation
    for _ in range(_ROTATION_MAX_SWEEPS):
        largest = 0.0
        for j in range(n_comp - 1):
            for k in range(j + 1, n_comp):
                x, y = loadings[:, j], loadings[:, k]
                u, v = x * x - y * y, 2.0 * x * y
                sum_u, sum_v = float(xp.sum(u)), float(xp.sum(v))
                num = 2.0 * (float(xp.sum(u * v)) - sum_u * sum_v / n_points)
                den = float(xp.sum(u * u)) - float(xp.sum(v * v)) - (sum_u**2 - sum_v**2) / n_points
                angle = atan2(num, den) / 4.0
                if abs(angle) <= _ROTATION_TOL:
                    continue
                c, s = cos(angle), sin(angle)
                loadings[:, j], loadings[:, k] = c * x + s * y, c * y - s * x
                col_j, col_k = rotation[:, j], rotation[:, k]
                rotation[:, j], rotation[:, k] = c * col_j + s * col_k, c * col_k - s * col_j
                largest = max(largest, abs(angle))
        if largest <= _ROTATION_TOL:
            break
    return rotation


def _grid(basis: Basis, n_points: int = _ROTATION_GRID) -> Array:
    """Return an evenly spaced evaluation grid over the domain of ``basis``."""
    lower, upper = basis.domain
    return _xp().linspace(lower, upper, n_points, dtype=_xp().float64)


class FPCA(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    r"""Functional principal component analysis, optionally with a roughness penalty.

    Replaces R's ``pca.fd`` and ``varmx.pca.fd``.

    Parameters
    ----------
    n : int or SupportsIndex, optional
        Number of harmonics to keep.  Default ``2``.  Any integer-like value
        (``numpy.int64``, ...) works; booleans and floats are rejected.
    lam : float or str, optional
        Roughness penalty on the harmonics.  ``0.0`` (default) is the
        unpenalised problem; ``"gcv"`` picks the penalty by leave-one-curve-out
        cross-validation of the ``n``-harmonic reconstruction error over the
        log grid ``10 ** -8, 10 ** -7, ..., 10 ** 8``.
    penalty : int or LDO, optional
        Roughness operator ``L``.  Default ``2`` (the second derivative).
    center : bool, optional
        Subtract the sample mean function before decomposing.  Default ``True``.
    basis : Basis, optional
        Basis to attach to a plain coefficient matrix.  Ignored when ``fit``
        receives an :class:`~fdatools.core.FData`.  ``None`` builds a cubic
        B-spline on ``(0, 1)`` with one basis function per column, which lets
        the estimator run inside a generic pipeline; pass the smoother's basis
        to decompose in the right metric.

    Attributes
    ----------
    harmonics : FData
        The ``n`` principal component functions.
    values : array
        The **full** eigenvalue spectrum, descending -- one entry per basis
        function (times the number of variables), as R reports it.  After
        :meth:`rotate` it holds the variance of each rotated component instead.
    scores : array
        Score matrix of shape ``(n_curves, n)``.  For multivariate curves the
        per-variable contributions are summed; see :attr:`scores_by_var`.
    varprop : array
        Proportion of total variance explained by each retained harmonic.
    mean_fd : FData
        The sample mean function (a zero function when ``center=False``).
    rotation : array or None
        The rotation applied by :meth:`rotate`, or ``None`` for an unrotated fit.

    Notes
    -----
    Eigenvector signs are mathematically arbitrary.  fdatools makes the
    coefficient sum of every harmonic positive -- the rule R's ``pca.fd``
    follows -- so repeated fits give identical output and match R's signs.

    Examples
    --------
    >>> import numpy as np
    >>> from fdatools import BSpline, FData
    >>> from fdatools.decomposition import FPCA
    >>> rng = np.random.default_rng(1)
    >>> fd = FData(rng.standard_normal((7, 40)), BSpline(domain=(0.0, 1.0), n_basis=7))
    >>> pca = FPCA(n=2).fit(fd)
    >>> pca.harmonics.n_curves
    2
    >>> bool(pca.varprop.sum() <= 1.0 + 1e-12)
    True
    """

    def __init__(
        self,
        n: SupportsIndex = 2,
        *,
        lam: float | str = 0.0,
        penalty: int | LDO = 2,
        center: bool = True,
        basis: Basis | None = None,
    ) -> None:
        self.n = n
        self.lam = lam
        self.penalty = penalty
        self.center = center
        self.basis = basis

    def __sklearn_tags__(self) -> Any:
        """Declare a dense, unsupervised transformer."""
        tags = super().__sklearn_tags__()
        tags.target_tags.required = False
        tags.input_tags.sparse = False
        return tags

    # ---------------------------------------------------------------- input

    def _as_fdata(self, X: Any, *, reset: bool) -> FData:  # noqa: N803
        """Interpret ``X`` as functional data, validating a plain array."""
        if isinstance(X, FData):
            if reset:
                self.n_features_in_ = X.basis.n_basis
            return X
        data = validate_data(self, X, reset=reset, ensure_min_samples=1)
        basis = self.basis if self.basis is not None else _default_basis(data.shape[1])
        if basis.n_basis != data.shape[1]:
            raise ValueError(
                f"X has {data.shape[1]} columns but the basis has {basis.n_basis} functions"
            )
        xp = _xp()
        return FData(xp.matrix_transpose(asarray(data, xp)), basis)

    def _centre(self, coefs: Array) -> tuple[Array, Array]:
        """Return ``(mean, coefs - mean)``; the mean is zero when not centring."""
        xp = _xp()
        mean = xp.mean(coefs, axis=1, keepdims=True) if self.center else xp.zeros_like(coefs[:, :1])
        return mean, coefs - mean

    def _gram(self, basis: Basis, n_vars: int) -> Array:
        """Return the (block-diagonal) Gram matrix for ``n_vars`` variables."""
        return _block_diagonal(asarray(to_numpy(basis.gram()), _xp()), n_vars)

    # ------------------------------------------------------------- fitting

    def _decompose(self, centred: Array, basis: Basis, lam: float) -> tuple[Array, Array]:
        """Return the eigenvalues and the ``(W + λR)``-orthonormal harmonic coefficients."""
        xp = _xp()
        n_curves = centred.shape[1]
        n_vars = _n_vars(centred)
        stacked = _stack_variables(centred)
        gram = self._gram(basis, n_vars)
        rough = _block_diagonal(asarray(to_numpy(basis.penalty(self.penalty)), xp), n_vars)
        cov = xp.matmul(stacked, xp.matrix_transpose(stacked)) / n_curves
        left = xp.matmul(gram, xp.matmul(cov, gram))
        left = 0.5 * (left + xp.matrix_transpose(left))
        right = gram + lam * rough
        right = 0.5 * (right + xp.matrix_transpose(right))
        mu, vec = _linalg.pencil_eigh(left, right)
        return xp.flip(asarray(mu, xp), axis=0), xp.flip(asarray(vec, xp), axis=1)

    def _reconstruction_sse(self, coefs: Array, basis: Basis, lam: float) -> float:
        """Return the leave-one-curve-out reconstruction error at penalty ``lam``.

        Each curve is projected, in the ``L²`` metric, onto the span of the ``n``
        harmonics estimated from the other curves; the error is the summed
        squared ``L²`` norm of the residuals.
        """
        xp = _xp()
        n_curves = coefs.shape[1]
        gram = self._gram(basis, _n_vars(coefs))
        stacked = _stack_variables(coefs)
        keep = min(_positive_int(self.n), stacked.shape[0])
        total = 0.0
        for i in range(n_curves):
            rest = xp.concat([coefs[:, :i], coefs[:, i + 1 :]], axis=1)
            mean, centred = self._centre(rest)
            _, vec = self._decompose(centred, basis, lam)
            harm = vec[:, :keep]
            residual = stacked[:, i : i + 1] - _stack_variables(mean)
            metric = xp.matmul(xp.matrix_transpose(harm), gram)
            weights = _linalg.solve_spd(xp.matmul(metric, harm), xp.matmul(metric, residual))
            err = residual - xp.matmul(harm, weights)
            total += float(xp.sum(err * xp.matmul(gram, err)))
        return total

    def _resolve_lambda(self, coefs: Array, basis: Basis) -> float:
        """Return the harmonic penalty requested by ``lam``."""
        lam = self.lam
        if isinstance(lam, str):
            if lam != "gcv":
                raise ValueError(f"lam must be a number or 'gcv', got {lam!r}")
            if coefs.shape[1] < 3:
                raise ValueError("lam='gcv' needs at least three curves")
            errors = [self._reconstruction_sse(coefs, basis, value) for value in _LAMBDA_GRID]
            return _LAMBDA_GRID[min(range(len(errors)), key=errors.__getitem__)]
        value = float(lam)
        if value < 0.0:
            raise ValueError(f"lam must be non-negative, got {value}")
        return value

    def fit(self, X: Any, y: Any = None) -> FPCA:  # noqa: N803
        """Decompose ``X`` into functional principal components.

        Parameters
        ----------
        X : FData or array of shape (n_samples, n_basis)
            Curves, or their basis coefficients one curve per row.
        y : ignored
            Present for API compatibility.

        Returns
        -------
        FPCA
            The fitted estimator.

        Raises
        ------
        ValueError
            If ``n`` is not a positive integer or ``lam`` is invalid.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import Fourier, FData
        >>> from fdatools.decomposition import FPCA
        >>> rng = np.random.default_rng(3)
        >>> fd = FData(rng.standard_normal((5, 20)), Fourier(domain=(0.0, 1.0), n_basis=5))
        >>> FPCA(n=2).fit(fd).harmonics.n_curves
        2
        """
        n_keep = _positive_int(self.n)
        xp = _xp()
        fd = self._as_fdata(X, reset=True)
        basis = fd.basis
        coefs = _coefficients(fd)
        n_curves = coefs.shape[1]
        n_vars = _n_vars(coefs)
        lam = self._resolve_lambda(coefs, basis)
        mean, centred = self._centre(coefs)

        values, vectors = self._decompose(centred, basis, lam)
        keep = min(n_keep, vectors.shape[1])
        harmonics = vectors[:, :keep]
        harmonics = harmonics * _positive_sum_signs(harmonics)

        gram = asarray(to_numpy(basis.gram()), xp)
        stacked = _stack_variables(centred)
        size = basis.n_basis
        per_var = xp.stack(
            [
                xp.matmul(
                    xp.matrix_transpose(stacked[k * size : (k + 1) * size, :]),
                    xp.matmul(gram, harmonics[k * size : (k + 1) * size, :]),
                )
                for k in range(n_vars)
            ],
            axis=2,
        )

        self.lam_ = lam
        self.n_components_ = keep
        self.n_curves_ = n_curves
        self.n_vars_ = n_vars
        self.values_ = values
        self.harmonics_ = FData(_unstack_variables(harmonics, size, n_vars), basis)
        self.scores_by_var_ = per_var
        self.scores_ = xp.sum(per_var, axis=2)
        total = float(xp.sum(values))
        self.total_variance_ = total
        self.varprop_ = values[:keep] / total if total != 0.0 else xp.zeros(keep)
        self.mean_fd_ = FData(mean, basis)
        self.rotation_: Array | None = None
        self.gram_ = self._gram(basis, n_vars)
        return self

    # ----------------------------------------------------------- accessors

    @property
    def harmonics(self) -> FData:
        """The principal component functions."""
        check_is_fitted(self)
        return self.harmonics_

    @property
    def values(self) -> Array:
        """The full eigenvalue spectrum, descending."""
        check_is_fitted(self)
        return self.values_

    @property
    def scores(self) -> Array:
        """Score matrix of shape ``(n_curves, n)``."""
        check_is_fitted(self)
        return self.scores_

    @property
    def scores_by_var(self) -> Array:
        """Score contributions of shape ``(n_curves, n, n_vars)``."""
        check_is_fitted(self)
        return self.scores_by_var_

    @property
    def varprop(self) -> Array:
        """Proportion of total variance carried by each retained harmonic."""
        check_is_fitted(self)
        return self.varprop_

    @property
    def mean_fd(self) -> FData:
        """The sample mean function."""
        check_is_fitted(self)
        return self.mean_fd_

    @property
    def rotation(self) -> Array | None:
        """The rotation applied by :meth:`rotate`, or ``None``."""
        check_is_fitted(self)
        return self.rotation_

    # ---------------------------------------------------------- transforms

    def transform(self, X: Any) -> Array:  # noqa: N803
        """Return the principal component scores of ``X``.

        Parameters
        ----------
        X : FData or array of shape (n_samples, n_basis)
            Curves to project onto the fitted harmonics.

        Returns
        -------
        array
            Scores of shape ``(n_samples, n)``, summed over variables.

        Raises
        ------
        ValueError
            If ``X`` does not carry the number of variables seen during ``fit``.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import BSpline, FData
        >>> from fdatools.decomposition import FPCA
        >>> rng = np.random.default_rng(4)
        >>> fd = FData(rng.standard_normal((6, 25)), BSpline(domain=(0.0, 1.0), n_basis=6))
        >>> pca = FPCA(n=2).fit(fd)
        >>> bool(np.allclose(pca.transform(fd), pca.scores))
        True
        """
        check_is_fitted(self)
        xp = _xp()
        coefs = _coefficients(self._as_fdata(X, reset=False))
        mean = _coefficients(self.mean_fd_)
        if _n_vars(coefs) != self.n_vars_ or coefs.shape[0] != mean.shape[0]:
            raise ValueError("X does not have the basis size and variables seen during fit")
        centred = _stack_variables(coefs - mean)
        harm = _stack_variables(_coefficients(self.harmonics_))
        return xp.matmul(xp.matrix_transpose(centred), xp.matmul(self.gram_, harm))

    def inverse_transform(self, X: Any) -> FData:  # noqa: N803
        """Rebuild curves from scores.

        Parameters
        ----------
        X : array of shape (n_samples, n)
            Principal component scores.

        Returns
        -------
        FData
            ``mean_fd`` plus the score-weighted harmonics.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import BSpline, FData
        >>> from fdatools.decomposition import FPCA
        >>> rng = np.random.default_rng(5)
        >>> fd = FData(rng.standard_normal((4, 30)), BSpline(domain=(0.0, 1.0), n_basis=4))
        >>> pca = FPCA(n=4).fit(fd)
        >>> bool(np.allclose(pca.inverse_transform(pca.scores).coefs, fd.coefs))
        True
        """
        check_is_fitted(self)
        xp = _xp()
        scores = asarray(X, xp)
        if scores.ndim == 1:
            scores = scores[None, :]
        harm = _coefficients(self.harmonics_)
        mean = _coefficients(self.mean_fd_)
        # Harmonics are (W + λR)-orthonormal, so the scores are converted back
        # to expansion weights through the Gram matrix of the harmonics.
        stacked = _stack_variables(harm)
        metric = xp.matmul(xp.matrix_transpose(stacked), xp.matmul(self.gram_, stacked))
        weights = _linalg.solve_spd(metric, xp.matrix_transpose(scores))
        if harm.ndim == 2:
            coefs = xp.matmul(harm, weights) + mean
        else:
            coefs = xp.stack(
                [xp.matmul(harm[:, :, k], weights) + mean[:, :, k] for k in range(harm.shape[2])],
                axis=2,
            )
        return FData(coefs, self.harmonics_.basis)

    # ------------------------------------------------------------ rotation

    def rotate(self, method: str = "varimax") -> FPCA:
        """Return a copy of this fit with rotated harmonics.

        Parameters
        ----------
        method : str, optional
            Only ``"varimax"`` is available; it maximises the variance of the
            squared harmonic values on a 501-point grid, exactly the quantity
            R's ``varmx.pca.fd`` rotates.

        Returns
        -------
        FPCA
            A new fitted estimator whose ``harmonics``, ``scores``, ``values``
            and ``varprop`` describe the rotated solution: ``values`` is the
            variance of each rotated component and ``varprop`` its share of the
            total variance of the unrotated spectrum.  ``mean_fd`` is unchanged.
            Rotated harmonics keep the positive-coefficient-sum sign rule.

        Raises
        ------
        ValueError
            If ``method`` is not ``"varimax"``.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import BSpline, FData
        >>> from fdatools.decomposition import FPCA
        >>> rng = np.random.default_rng(6)
        >>> fd = FData(rng.standard_normal((8, 30)), BSpline(domain=(0.0, 1.0), n_basis=8))
        >>> rotated = FPCA(n=3).fit(fd).rotate("varimax")
        >>> bool(np.allclose(rotated.rotation.T @ rotated.rotation, np.eye(3)))
        True
        """
        check_is_fitted(self)
        if method != "varimax":
            raise ValueError(f"unknown rotation {method!r}; only 'varimax' is available")
        xp = _xp()
        basis = self.harmonics_.basis
        harm = _coefficients(self.harmonics_)
        design = asarray(to_numpy(basis(_grid(basis))), xp)
        if harm.ndim == 2:
            loadings = xp.matmul(design, harm)
        else:
            loadings = xp.concat(
                [xp.matmul(design, harm[:, :, k]) for k in range(harm.shape[2])], axis=0
            )
        rot = _varimax_rotation(loadings)
        rot = rot * _positive_sum_signs(xp.matmul(_stack_variables(harm), rot))

        rotated = self.__class__(**self.get_params())
        rotated.__dict__.update({k: v for k, v in self.__dict__.items() if k.endswith("_")})
        if harm.ndim == 2:
            rotated.harmonics_ = FData(xp.matmul(harm, rot), basis)
        else:
            rotated.harmonics_ = FData(
                xp.stack([xp.matmul(harm[:, :, k], rot) for k in range(harm.shape[2])], axis=2),
                basis,
            )
        by_var = self.scores_by_var_
        rotated.scores_by_var_ = xp.stack(
            [xp.matmul(by_var[:, :, k], rot) for k in range(by_var.shape[2])], axis=2
        )
        rotated.scores_ = xp.sum(rotated.scores_by_var_, axis=2)
        rotated.values_ = xp.sum(rotated.scores_**2, axis=0) / self.n_curves_
        total = self.total_variance_
        rotated.varprop_ = (
            rotated.values_ / total if total != 0.0 else xp.zeros_like(rotated.values_)
        )
        rotated.rotation_ = rot
        return rotated

    # ------------------------------------------------------------ plotting

    def plot(self, ax: Any = None, *, n_points: int = 201) -> list[Any]:
        """Plot each harmonic as a perturbation of the mean.

        Every panel shows the mean function and the mean plus and minus the
        harmonic scaled by the standard deviation of its scores -- the display
        R's ``plot.pca.fd`` produces.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or sequence of Axes, optional
            Axes to draw on, one per harmonic.  ``None`` creates one row of
            panels.
        n_points : int, optional
            Number of evaluation points.  Default ``201``.

        Returns
        -------
        list of matplotlib.axes.Axes
            The axes drawn on.

        Raises
        ------
        ValueError
            If fewer axes than harmonics are supplied.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")
        >>> import numpy as np
        >>> from fdatools import BSpline, FData
        >>> from fdatools.decomposition import FPCA
        >>> rng = np.random.default_rng(7)
        >>> fd = FData(rng.standard_normal((6, 20)), BSpline(domain=(0.0, 1.0), n_basis=6))
        >>> len(FPCA(n=2).fit(fd).plot())
        2
        """
        check_is_fitted(self)
        import matplotlib.pyplot as plt

        xp = _xp()
        keep = self.n_components_
        if ax is None:
            _, grid_axes = plt.subplots(1, keep, figsize=(4.0 * keep, 3.0), squeeze=False)
            axes = list(grid_axes.ravel())
        elif hasattr(ax, "ravel"):
            axes = list(ax.ravel())
        elif isinstance(ax, (list, tuple)):
            axes = list(ax)
        else:
            axes = [ax]
        if len(axes) < keep:
            raise ValueError(f"{keep} harmonics need {keep} axes, got {len(axes)}")
        grid = _grid(self.harmonics_.basis, n_points)
        mean = xp.reshape(asarray(to_numpy(self.mean_fd_(grid)), xp), (n_points, -1))
        harm = xp.reshape(asarray(to_numpy(self.harmonics_(grid)), xp), (n_points, keep, -1))
        spread = xp.sum(self.scores_**2, axis=0) / self.n_curves_
        for j in range(keep):
            effect = sqrt(max(float(spread[j]), 0.0)) * harm[:, j, :]
            axes[j].plot(grid, mean, color="black", label="mean")
            axes[j].plot(grid, mean + effect, color="tab:blue", linestyle="--", label="+")
            axes[j].plot(grid, mean - effect, color="tab:red", linestyle=":", label="-")
            axes[j].set_title(f"PC {j + 1} ({100 * float(self.varprop_[j]):.1f}%)")
        return axes


class FCCA(BaseEstimator):  # type: ignore[misc]
    r"""Functional canonical correlation analysis.

    Replaces R's ``cca.fd``.

    Parameters
    ----------
    n : int or SupportsIndex, optional
        Number of canonical variate pairs to keep.  Default ``2``.  Any
        integer-like value (``numpy.int64``, ...) works; booleans and floats
        are rejected.
    lam1, lam2 : float, optional
        Roughness penalties on the first and second set of canonical weights.
        Without a penalty the problem is degenerate whenever there are fewer
        curves than basis functions (every correlation is one).
    penalty : int or LDO, optional
        Roughness operator ``L``.  Default ``2``.
    center : bool, optional
        Subtract each sample mean before decomposing.  Default ``True``.

    Attributes
    ----------
    weights1, weights2 : FData
        Canonical weight functions, normalised to unit :math:`L^2` norm.  The
        sign of each pair is fixed so the first weight function has a positive
        coefficient sum; both members of a pair share that sign, so their
        scores stay positively correlated.
    correlations : array
        The **full** canonical correlation spectrum, descending, as R reports it.
    scores1, scores2 : array
        Canonical variate scores of shape ``(n_curves, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> from fdatools import BSpline, FData
    >>> from fdatools.decomposition import FCCA
    >>> rng = np.random.default_rng(2)
    >>> basis = BSpline(domain=(0.0, 1.0), n_basis=6)
    >>> a = rng.standard_normal((6, 25))
    >>> fd1 = FData(a, basis)
    >>> fd2 = FData(a + 0.1 * rng.standard_normal((6, 25)), basis)
    >>> cca = FCCA(n=2, lam1=1e-6, lam2=1e-6).fit(fd1, fd2)
    >>> bool(cca.correlations[0] > 0.9)
    True
    """

    def __init__(
        self,
        n: SupportsIndex = 2,
        *,
        lam1: float = 0.0,
        lam2: float = 0.0,
        penalty: int | LDO = 2,
        center: bool = True,
    ) -> None:
        self.n = n
        self.lam1 = lam1
        self.lam2 = lam2
        self.penalty = penalty
        self.center = center

    def __sklearn_tags__(self) -> Any:
        """Declare a dense estimator whose second view arrives as ``y``."""
        tags = super().__sklearn_tags__()
        tags.target_tags.required = True
        tags.input_tags.sparse = False
        return tags

    @staticmethod
    def _univariate(fd: FData) -> Array:
        """Return the ``(n_basis, n_curves)`` coefficients of univariate curves."""
        coefs = _coefficients(fd)
        if coefs.ndim != 2:
            raise ValueError("FCCA needs univariate curves")
        return coefs

    def fit(self, X: FData, y: FData) -> FCCA:  # noqa: N803
        """Find the canonical weight functions of two samples of curves.

        Parameters
        ----------
        X, y : FData
            Two samples observed on the same curves, in that order.

        Returns
        -------
        FCCA
            The fitted estimator.

        Raises
        ------
        TypeError
            If either argument is not an :class:`~fdatools.core.FData`.
        ValueError
            If the samples hold different numbers of curves, are multivariate,
            or a parameter is out of range.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import Fourier, FData
        >>> from fdatools.decomposition import FCCA
        >>> rng = np.random.default_rng(8)
        >>> basis = Fourier(domain=(0.0, 1.0), n_basis=5)
        >>> x = rng.standard_normal((5, 40))
        >>> cca = FCCA(n=1, lam1=1e-3, lam2=1e-3).fit(FData(x, basis), FData(-x, basis))
        >>> cca.scores1.shape
        (40, 1)
        """
        if not isinstance(X, FData) or not isinstance(y, FData):
            raise TypeError("FCCA.fit needs two FData arguments")
        n_keep = _positive_int(self.n)
        if float(self.lam1) < 0.0 or float(self.lam2) < 0.0:
            raise ValueError("lam1 and lam2 must be non-negative")
        xp = _xp()
        cx, cy = self._univariate(X), self._univariate(y)
        if cx.shape[1] != cy.shape[1]:
            raise ValueError(
                f"the two samples hold {cx.shape[1]} and {cy.shape[1]} curves; they must match"
            )
        mean1 = xp.mean(cx, axis=1, keepdims=True) if self.center else xp.zeros_like(cx[:, :1])
        mean2 = xp.mean(cy, axis=1, keepdims=True) if self.center else xp.zeros_like(cy[:, :1])
        cx, cy = cx - mean1, cy - mean2
        basis1, basis2 = X.basis, y.basis
        n_curves = cx.shape[1]
        w1 = asarray(to_numpy(basis1.gram()), xp)
        w2 = asarray(to_numpy(basis2.gram()), xp)
        r1 = asarray(to_numpy(basis1.penalty(self.penalty)), xp)
        r2 = asarray(to_numpy(basis2.penalty(self.penalty)), xp)
        proj1 = xp.matmul(w1, cx)
        proj2 = xp.matmul(w2, cy)
        sxy = xp.matmul(proj1, xp.matrix_transpose(proj2)) / n_curves
        sxx = xp.matmul(proj1, xp.matrix_transpose(proj1)) / n_curves
        syy = xp.matmul(proj2, xp.matrix_transpose(proj2)) / n_curves

        p1, p2 = basis1.n_basis, basis2.n_basis
        zero12 = xp.zeros((p1, p2), dtype=xp.float64)
        left = xp.concat(
            [
                xp.concat([xp.zeros((p1, p1), dtype=xp.float64), sxy], axis=1),
                xp.concat([xp.matrix_transpose(sxy), xp.zeros((p2, p2), dtype=xp.float64)], axis=1),
            ],
            axis=0,
        )
        right = xp.concat(
            [
                xp.concat([sxx + float(self.lam1) * r1, zero12], axis=1),
                xp.concat([xp.matrix_transpose(zero12), syy + float(self.lam2) * r2], axis=1),
            ],
            axis=0,
        )
        right = 0.5 * (right + xp.matrix_transpose(right))
        rho, vec = _linalg.pencil_eigh(left, right)
        rho = xp.flip(asarray(rho, xp), axis=0)
        vec = xp.flip(asarray(vec, xp), axis=1)

        keep = min(n_keep, p1, p2)
        a = vec[:p1, :keep]
        b = vec[p1:, :keep]
        tiny = 1e-300
        norm1 = xp.sqrt(xp.clip(xp.sum(a * xp.matmul(w1, a), axis=0), tiny, None))
        norm2 = xp.sqrt(xp.clip(xp.sum(b * xp.matmul(w2, b), axis=0), tiny, None))
        signs = _positive_sum_signs(a)
        a = a / norm1 * signs
        b = b / norm2 * signs

        self.n_components_ = keep
        self.correlations_ = rho[: min(p1, p2)]
        self.weights1_ = FData(a, basis1)
        self.weights2_ = FData(b, basis2)
        self.scores1_ = xp.matmul(xp.matrix_transpose(cx), xp.matmul(w1, a))
        self.scores2_ = xp.matmul(xp.matrix_transpose(cy), xp.matmul(w2, b))
        self.grams_ = (w1, w2)
        self.means_ = (mean1, mean2)
        return self

    @property
    def weights1(self) -> FData:
        """Canonical weight functions of the first sample."""
        check_is_fitted(self)
        return self.weights1_

    @property
    def weights2(self) -> FData:
        """Canonical weight functions of the second sample."""
        check_is_fitted(self)
        return self.weights2_

    @property
    def correlations(self) -> Array:
        """The full canonical correlation spectrum, descending."""
        check_is_fitted(self)
        return self.correlations_

    @property
    def scores1(self) -> Array:
        """Canonical variate scores of the first sample."""
        check_is_fitted(self)
        return self.scores1_

    @property
    def scores2(self) -> Array:
        """Canonical variate scores of the second sample."""
        check_is_fitted(self)
        return self.scores2_

    def transform(self, X: FData, y: FData) -> tuple[Array, Array]:  # noqa: N803
        """Return the canonical variate scores of two new samples.

        Parameters
        ----------
        X, y : FData
            Curves to project, in the order given to :meth:`fit`.

        Returns
        -------
        tuple of array
            Scores of shape ``(n_curves, n)`` for each sample.

        Raises
        ------
        ValueError
            If a sample's basis size differs from the one seen during ``fit``.

        Examples
        --------
        >>> import numpy as np
        >>> from fdatools import Fourier, FData
        >>> from fdatools.decomposition import FCCA
        >>> rng = np.random.default_rng(9)
        >>> basis = Fourier(domain=(0.0, 1.0), n_basis=5)
        >>> x = FData(rng.standard_normal((5, 30)), basis)
        >>> y = FData(rng.standard_normal((5, 30)), basis)
        >>> cca = FCCA(n=2, lam1=1e-2, lam2=1e-2).fit(x, y)
        >>> bool(np.allclose(cca.transform(x, y)[0], cca.scores1))
        True
        """
        check_is_fitted(self)
        xp = _xp()
        w1, w2 = self.grams_
        m1, m2 = self.means_
        cx, cy = self._univariate(X), self._univariate(y)
        if cx.shape[0] != w1.shape[0] or cy.shape[0] != w2.shape[0]:
            raise ValueError("X and y must use the basis sizes seen during fit")
        a = _coefficients(self.weights1_)
        b = _coefficients(self.weights2_)
        return (
            xp.matmul(xp.matrix_transpose(cx - m1), xp.matmul(w1, a)),
            xp.matmul(xp.matrix_transpose(cy - m2), xp.matmul(w2, b)),
        )
