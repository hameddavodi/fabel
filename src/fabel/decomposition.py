r"""Functional principal components and canonical correlation.

Two estimators replace R ``fda``'s ``pca.fd``, ``varmx.pca.fd`` and ``cca.fd``.
Both work on the basis coefficients of an :class:`~fabel.core.FData` and both
follow the scikit-learn estimator API, so they drop straight into a
:class:`~sklearn.pipeline.Pipeline` behind :class:`~fabel.smoothing.Smoother`.

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

Examples
--------
>>> import numpy as np
>>> import fabel as fb
>>> rng = np.random.default_rng(0)
>>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=8)
>>> fd = fb.FData(rng.standard_normal((8, 30)), basis)
>>> pca = fb.FPCA(n=3).fit(fd)
>>> pca.scores.shape
(30, 3)
>>> bool(np.all(np.diff(pca.values) <= 1e-12))
True
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted, validate_data

from fabel import _linalg
from fabel._backend import to_numpy
from fabel._operator import LDO
from fabel.basis import Basis, BSpline
from fabel.core import FData

if TYPE_CHECKING:  # pragma: no cover - typing only
    NDArray = np.ndarray[Any, np.dtype[Any]]
else:
    NDArray = Any

__all__ = ["FCCA", "FPCA"]

#: Number of grid points R's ``varmx.pca.fd`` rotates on, and the value Fabel
#: uses so that a rotation is comparable with R's.
_ROTATION_GRID = 501

#: Angle below which the varimax sweep stops rotating a pair of columns.
_ROTATION_TOL = 1e-12

#: Maximum number of varimax sweeps.
_ROTATION_MAX_SWEEPS = 1000

#: Log-spaced grid searched when ``lam="gcv"``.
_LAMBDA_GRID = 10.0 ** np.linspace(-8.0, 8.0, 17)

#: Default B-spline order used when an estimator has to invent a basis.
_DEFAULT_ORDER = 4


def _default_basis(n_basis: int) -> Basis:
    """Return the fallback basis for a coefficient matrix with ``n_basis`` columns."""
    return BSpline(domain=(0.0, 1.0), n_basis=n_basis, order=min(_DEFAULT_ORDER, n_basis))


def _block_diagonal(matrix: NDArray, repeats: int) -> NDArray:
    """Return ``matrix`` repeated ``repeats`` times down a block diagonal."""
    if repeats == 1:
        return matrix
    return np.kron(np.eye(repeats), matrix)


def _stack_variables(coefs: NDArray) -> NDArray:
    """Flatten a ``(n_basis, n_curves[, n_vars])`` array to ``(n_basis * n_vars, n_curves)``."""
    if coefs.ndim == 2:
        return coefs
    return np.concatenate([coefs[:, :, k] for k in range(coefs.shape[2])], axis=0)


def _unstack_variables(stacked: NDArray, n_basis: int, n_vars: int) -> NDArray:
    """Undo :func:`_stack_variables` for a matrix of column vectors."""
    if n_vars == 1:
        return stacked
    parts = [stacked[k * n_basis : (k + 1) * n_basis, :] for k in range(n_vars)]
    return np.stack(parts, axis=2)


def _sign_align(vectors: NDArray) -> NDArray:
    """Fix the arbitrary sign of each column.

    The eigenvector sign is undetermined, so a deterministic public API needs a
    rule.  Fabel makes the coefficient of largest magnitude positive; ties are
    broken by the first such coefficient.
    """
    if vectors.size == 0:
        return vectors
    pivot = np.argmax(np.abs(vectors), axis=0)
    signs = np.sign(vectors[pivot, np.arange(vectors.shape[1])])
    signs[signs == 0.0] = 1.0
    return np.asarray(vectors * signs, dtype=np.float64)


def _varimax_rotation(values: NDArray) -> NDArray:
    """Return the varimax rotation of the columns of ``values``.

    Uses the classical pairwise-Jacobi sweep on the *raw* (un-normalised)
    loadings, maximising ``sum_j [ sum_i a_ij^4 - (sum_i a_ij^2)^2 / n ]``.

    Parameters
    ----------
    values : numpy.ndarray
        Loadings of shape ``(n_points, n_components)``.

    Returns
    -------
    numpy.ndarray
        Orthogonal ``(n_components, n_components)`` rotation ``T`` such that
        ``values @ T`` maximises the varimax criterion.
    """
    n_points, n_comp = values.shape
    rotation = np.eye(n_comp)
    loadings = np.array(values, dtype=np.float64, copy=True)
    if n_comp < 2:
        return rotation
    for _ in range(_ROTATION_MAX_SWEEPS):
        largest = 0.0
        for j in range(n_comp - 1):
            for k in range(j + 1, n_comp):
                x, y = loadings[:, j], loadings[:, k]
                u, v = x * x - y * y, 2.0 * x * y
                sum_u, sum_v = float(u.sum()), float(v.sum())
                num = 2.0 * (float((u * v).sum()) - sum_u * sum_v / n_points)
                den = float((u * u).sum()) - float((v * v).sum()) - (sum_u**2 - sum_v**2) / n_points
                angle = float(np.arctan2(num, den)) / 4.0
                if abs(angle) <= _ROTATION_TOL:
                    continue
                cos, sin = np.cos(angle), np.sin(angle)
                loadings[:, j], loadings[:, k] = cos * x + sin * y, cos * y - sin * x
                col_j, col_k = rotation[:, j].copy(), rotation[:, k].copy()
                rotation[:, j], rotation[:, k] = (
                    cos * col_j + sin * col_k,
                    cos * col_k - sin * col_j,
                )
                largest = max(largest, abs(angle))
        if largest <= _ROTATION_TOL:
            break
    return rotation


def _grid(basis: Basis, n_points: int = _ROTATION_GRID) -> NDArray:
    """Return an evenly spaced evaluation grid over the domain of ``basis``."""
    lower, upper = basis.domain
    return np.linspace(lower, upper, n_points)


class FPCA(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    r"""Functional principal component analysis, optionally with a roughness penalty.

    Replaces R's ``pca.fd`` and ``varmx.pca.fd``.

    Parameters
    ----------
    n : int, optional
        Number of harmonics to keep.  Default ``2``.
    lam : float or str, optional
        Roughness penalty on the harmonics.  ``0.0`` (default) is the
        unpenalised problem; ``"gcv"`` picks the penalty by leave-one-curve-out
        cross-validation of the ``n``-harmonic reconstruction error over the
        log grid ``10 ** linspace(-8, 8, 17)``.
    penalty : int or LDO, optional
        Roughness operator ``L``.  Default ``2`` (the second derivative).
    center : bool, optional
        Subtract the sample mean function before decomposing.  Default ``True``.
    basis : Basis, optional
        Basis to attach to a plain coefficient matrix.  Ignored when ``fit``
        receives an :class:`~fabel.core.FData`.  ``None`` builds a cubic
        B-spline on ``(0, 1)`` with one basis function per column, which lets
        the estimator run inside a generic pipeline.

    Attributes
    ----------
    harmonics : FData
        The ``n`` principal component functions.
    values : numpy.ndarray
        The **full** eigenvalue spectrum, descending -- one entry per basis
        function (times the number of variables), as R reports it.
    scores : numpy.ndarray
        Score matrix of shape ``(n_curves, n)``.  For multivariate curves the
        per-variable contributions are summed; see :attr:`scores_by_var`.
    varprop : numpy.ndarray
        Proportion of total variance explained by each retained harmonic.
    mean_fd : FData
        The sample mean function (a zero function when ``center=False``).
    rotation : numpy.ndarray or None
        The rotation applied by :meth:`rotate`, or ``None`` for an unrotated fit.

    Notes
    -----
    Eigenvector signs are mathematically arbitrary.  Fabel fixes them by making
    the largest-magnitude coefficient of every harmonic positive, so repeated
    fits on the same data give identical output.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> rng = np.random.default_rng(1)
    >>> fd = fb.FData(rng.standard_normal((7, 40)), fb.BSpline(domain=(0.0, 1.0), n_basis=7))
    >>> pca = fb.FPCA(n=2).fit(fd)
    >>> pca.harmonics.n_curves
    2
    >>> bool(pca.varprop.sum() <= 1.0 + 1e-12)
    True
    """

    def __init__(
        self,
        n: int = 2,
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
        return FData(np.asarray(data, dtype=np.float64).T, basis)

    # ------------------------------------------------------------- fitting

    def _decompose(self, coefs: NDArray, basis: Basis, lam: float) -> tuple[NDArray, NDArray]:
        """Return the eigenvalues and the ``B``-orthonormal harmonic coefficients."""
        n_curves = coefs.shape[1]
        stacked = _stack_variables(coefs)
        n_vars = 1 if coefs.ndim == 2 else coefs.shape[2]
        gram = _block_diagonal(np.asarray(basis.gram(), dtype=np.float64), n_vars)
        rough = _block_diagonal(np.asarray(basis.penalty(self.penalty), dtype=np.float64), n_vars)
        cov = stacked @ stacked.T / n_curves
        left = gram @ cov @ gram
        left = 0.5 * (left + left.T)
        right = gram + lam * rough
        right = 0.5 * (right + right.T)
        mu, vec = _linalg.pencil_eigh(left, right)
        return mu[::-1].copy(), vec[:, ::-1].copy()

    def _reconstruction_sse(self, coefs: NDArray, basis: Basis, lam: float) -> float:
        """Return the leave-one-curve-out reconstruction error at penalty ``lam``."""
        n_curves = coefs.shape[1]
        n_vars = 1 if coefs.ndim == 2 else coefs.shape[2]
        gram = _block_diagonal(np.asarray(basis.gram(), dtype=np.float64), n_vars)
        stacked = _stack_variables(coefs)
        keep = min(self.n, stacked.shape[0])
        total = 0.0
        for i in range(n_curves):
            rest = np.delete(coefs, i, axis=1)
            mean = rest.mean(axis=1, keepdims=True) if self.center else np.zeros_like(rest[:, :1])
            _, vec = self._decompose(rest - mean, basis, lam)
            harm = vec[:, :keep]
            residual = stacked[:, i : i + 1] - _stack_variables(mean)
            fitted = harm @ (harm.T @ gram @ residual)
            err = residual - fitted
            total += float(err.T @ gram @ err)
        return total

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
        """
        fd = self._as_fdata(X, reset=True)
        basis = fd.basis
        coefs = np.asarray(to_numpy(fd.coefs), dtype=np.float64)
        if coefs.ndim == 1:
            coefs = coefs[:, None]
        n_curves = coefs.shape[1]
        n_vars = 1 if coefs.ndim == 2 else coefs.shape[2]
        mean = coefs.mean(axis=1, keepdims=True) if self.center else np.zeros_like(coefs[:, :1])
        centred = coefs - mean

        lam = self.lam
        if isinstance(lam, str):
            if lam != "gcv":
                raise ValueError(f"lam must be a number or 'gcv', got {lam!r}")
            scores = [self._reconstruction_sse(centred, basis, float(c)) for c in _LAMBDA_GRID]
            lam = float(_LAMBDA_GRID[int(np.argmin(scores))])
        lam = float(lam)
        if lam < 0.0:
            raise ValueError(f"lam must be non-negative, got {lam}")

        values, vectors = self._decompose(centred, basis, lam)
        keep = min(self.n, vectors.shape[1])
        harmonics = _sign_align(vectors[:, :keep])

        gram = _block_diagonal(np.asarray(basis.gram(), dtype=np.float64), n_vars)
        stacked = _stack_variables(centred)
        per_var = np.stack(
            [
                stacked[k * basis.n_basis : (k + 1) * basis.n_basis, :].T
                @ np.asarray(basis.gram(), dtype=np.float64)
                @ harmonics[k * basis.n_basis : (k + 1) * basis.n_basis, :]
                for k in range(n_vars)
            ],
            axis=2,
        )

        self.lam_ = lam
        self.n_components_ = keep
        self.n_curves_ = n_curves
        self.values_ = values
        self.harmonics_ = FData(_unstack_variables(harmonics, basis.n_basis, n_vars), basis)
        self.scores_by_var_ = per_var
        self.scores_ = per_var.sum(axis=2)
        total = float(values.sum())
        self.varprop_ = values[:keep] / total if total != 0.0 else np.zeros(keep)
        self.mean_fd_ = FData(_unstack_variables(mean, basis.n_basis, n_vars), basis)
        self.rotation_: NDArray | None = None
        self._gram = gram
        return self

    # ----------------------------------------------------------- accessors

    @property
    def harmonics(self) -> FData:
        """The principal component functions."""
        check_is_fitted(self)
        return self.harmonics_

    @property
    def values(self) -> NDArray:
        """The full eigenvalue spectrum, descending."""
        check_is_fitted(self)
        return self.values_

    @property
    def scores(self) -> NDArray:
        """Score matrix of shape ``(n_curves, n)``."""
        check_is_fitted(self)
        return self.scores_

    @property
    def scores_by_var(self) -> NDArray:
        """Score contributions of shape ``(n_curves, n, n_vars)``."""
        check_is_fitted(self)
        return self.scores_by_var_

    @property
    def varprop(self) -> NDArray:
        """Proportion of total variance carried by each retained harmonic."""
        check_is_fitted(self)
        return self.varprop_

    @property
    def mean_fd(self) -> FData:
        """The sample mean function."""
        check_is_fitted(self)
        return self.mean_fd_

    @property
    def rotation(self) -> NDArray | None:
        """The rotation applied by :meth:`rotate`, or ``None``."""
        check_is_fitted(self)
        return self.rotation_

    # ---------------------------------------------------------- transforms

    def transform(self, X: Any) -> NDArray:  # noqa: N803
        """Return the principal component scores of ``X``.

        Parameters
        ----------
        X : FData or array of shape (n_samples, n_basis)
            Curves to project onto the fitted harmonics.

        Returns
        -------
        numpy.ndarray
            Scores of shape ``(n_samples, n)``.
        """
        check_is_fitted(self)
        fd = self._as_fdata(X, reset=False)
        basis = self.harmonics_.basis
        coefs = np.asarray(to_numpy(fd.coefs), dtype=np.float64)
        if coefs.ndim == 1:
            coefs = coefs[:, None]
        mean = np.asarray(to_numpy(self.mean_fd_.coefs), dtype=np.float64)
        if mean.ndim == 1:
            mean = mean[:, None]
        centred = _stack_variables(coefs - mean)
        harm = _stack_variables(np.asarray(to_numpy(self.harmonics_.coefs), dtype=np.float64))
        gram = self._gram
        if centred.shape[0] != gram.shape[0]:
            raise ValueError("X does not carry the number of variables seen during fit")
        del basis
        return np.asarray(centred.T @ gram @ harm, dtype=np.float64)

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
        """
        check_is_fitted(self)
        scores = np.atleast_2d(np.asarray(X, dtype=np.float64))
        harm = np.asarray(to_numpy(self.harmonics_.coefs), dtype=np.float64)
        mean = np.asarray(to_numpy(self.mean_fd_.coefs), dtype=np.float64)
        if harm.ndim == 2:
            coefs = harm @ scores.T + mean
        else:
            coefs = np.einsum("bkv,nk->bnv", harm, scores) + mean
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
            and ``varprop`` describe the rotated solution.  ``mean_fd`` and the
            eigenvalue spectrum of the unrotated fit are unchanged.
        """
        check_is_fitted(self)
        if method != "varimax":
            raise ValueError(f"unknown rotation {method!r}; only 'varimax' is available")
        basis = self.harmonics_.basis
        grid = _grid(basis)
        harm = np.asarray(to_numpy(self.harmonics_.coefs), dtype=np.float64)
        design = np.asarray(basis(grid), dtype=np.float64)
        if harm.ndim == 2:
            loadings = design @ harm
        else:
            loadings = np.concatenate(
                [design @ harm[:, :, k] for k in range(harm.shape[2])], axis=0
            )
        rot = _varimax_rotation(loadings)

        rotated = self.__class__(**self.get_params())
        rotated.__dict__.update({k: v for k, v in self.__dict__.items() if k.endswith("_")})
        rotated._gram = self._gram
        if harm.ndim == 2:
            rotated.harmonics_ = FData(harm @ rot, basis)
        else:
            rotated.harmonics_ = FData(np.einsum("bkv,kj->bjv", harm, rot), basis)
        rotated.scores_by_var_ = np.einsum("nkv,kj->njv", self.scores_by_var_, rot)
        rotated.scores_ = rotated.scores_by_var_.sum(axis=2)
        rotated.values_ = (rotated.scores_**2).sum(axis=0) / self.n_curves_
        total = float(self.values_.sum())
        rotated.varprop_ = (
            rotated.values_ / total if total != 0.0 else np.zeros_like(rotated.values_)
        )
        rotated.rotation_ = rot
        return rotated

    # ------------------------------------------------------------ plotting

    def plot(self, ax: Any = None, *, n_points: int = 201) -> Any:
        """Plot each harmonic as a perturbation of the mean.

        Every panel shows the mean function and the mean plus and minus a
        multiple of the harmonic, the display R's ``plot.pca.fd`` produces.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or sequence of Axes, optional
            Axes to draw on.  ``None`` creates one row of panels.
        n_points : int, optional
            Number of evaluation points.  Default ``201``.

        Returns
        -------
        numpy.ndarray
            The axes drawn on.
        """
        check_is_fitted(self)
        import matplotlib.pyplot as plt

        basis = self.harmonics_.basis
        grid = _grid(basis, n_points)
        keep = self.n_components_
        if ax is None:
            _, axes = plt.subplots(1, keep, figsize=(4.0 * keep, 3.0), squeeze=False)
            axes = axes.ravel()
        else:
            axes = np.atleast_1d(np.asarray(ax, dtype=object)).ravel()
        mean = np.asarray(to_numpy(self.mean_fd_(grid)), dtype=np.float64).reshape(n_points, -1)
        harm = np.asarray(to_numpy(self.harmonics_(grid)), dtype=np.float64)
        harm = harm.reshape(n_points, keep, -1)
        for j in range(keep):
            size = float(np.sqrt(max(self.values_[j], 0.0)))
            effect = size * harm[:, j, :]
            axes[j].plot(grid, mean, color="black", label="mean")
            axes[j].plot(grid, mean + effect, color="tab:blue", linestyle="--", label="+")
            axes[j].plot(grid, mean - effect, color="tab:red", linestyle=":", label="-")
            axes[j].set_title(f"PC {j + 1} ({100 * self.varprop_[j]:.1f}%)")
        return axes


class FCCA(BaseEstimator):  # type: ignore[misc]
    r"""Functional canonical correlation analysis.

    Replaces R's ``cca.fd``.

    Parameters
    ----------
    n : int, optional
        Number of canonical variate pairs to keep.  Default ``2``.
    lam1, lam2 : float, optional
        Roughness penalties on the first and second set of canonical weights.
    penalty : int or LDO, optional
        Roughness operator ``L``.  Default ``2``.
    center : bool, optional
        Subtract each sample mean before decomposing.  Default ``True``.

    Attributes
    ----------
    weights1, weights2 : FData
        Canonical weight functions, normalised to unit :math:`L^2` norm.
    correlations : numpy.ndarray
        The **full** canonical correlation spectrum, descending, as R reports it.
    scores1, scores2 : numpy.ndarray
        Canonical variate scores of shape ``(n_curves, n)``.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> rng = np.random.default_rng(2)
    >>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=6)
    >>> a = rng.standard_normal((6, 25))
    >>> fd1 = fb.FData(a, basis)
    >>> fd2 = fb.FData(a + 0.1 * rng.standard_normal((6, 25)), basis)
    >>> cca = fb.FCCA(n=2, lam1=1e-6, lam2=1e-6).fit(fd1, fd2)
    >>> bool(cca.correlations[0] > 0.9)
    True
    """

    def __init__(
        self,
        n: int = 2,
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
    def _centred(fd: FData, center: bool) -> tuple[NDArray, Basis]:
        """Return the (optionally centred) coefficient matrix of ``fd``."""
        coefs = np.asarray(to_numpy(fd.coefs), dtype=np.float64)
        if coefs.ndim == 1:
            coefs = coefs[:, None]
        if coefs.ndim != 2:
            raise ValueError("FCCA needs univariate curves")
        if center:
            coefs = coefs - coefs.mean(axis=1, keepdims=True)
        return coefs, fd.basis

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
        """
        if not isinstance(X, FData) or not isinstance(y, FData):
            raise TypeError("FCCA.fit needs two FData arguments")
        cx, basis1 = self._centred(X, self.center)
        cy, basis2 = self._centred(y, self.center)
        if cx.shape[1] != cy.shape[1]:
            raise ValueError(
                f"the two samples hold {cx.shape[1]} and {cy.shape[1]} curves; they must match"
            )
        n_curves = cx.shape[1]
        w1 = np.asarray(basis1.gram(), dtype=np.float64)
        w2 = np.asarray(basis2.gram(), dtype=np.float64)
        r1 = np.asarray(basis1.penalty(self.penalty), dtype=np.float64)
        r2 = np.asarray(basis2.penalty(self.penalty), dtype=np.float64)
        sxy = w1 @ cx @ cy.T @ w2 / n_curves
        sxx = w1 @ cx @ cx.T @ w1 / n_curves
        syy = w2 @ cy @ cy.T @ w2 / n_curves

        p1, p2 = basis1.n_basis, basis2.n_basis
        left = np.zeros((p1 + p2, p1 + p2))
        left[:p1, p1:] = sxy
        left[p1:, :p1] = sxy.T
        right = np.zeros_like(left)
        right[:p1, :p1] = sxx + float(self.lam1) * r1
        right[p1:, p1:] = syy + float(self.lam2) * r2
        left = 0.5 * (left + left.T)
        right = 0.5 * (right + right.T)
        rho, vec = _linalg.pencil_eigh(left, right)
        rho, vec = rho[::-1].copy(), vec[:, ::-1].copy()

        keep = min(self.n, p1, p2)
        spectrum = rho[: min(p1, p2)]
        a = vec[:p1, :keep]
        b = vec[p1:, :keep]
        a = a / np.sqrt(np.maximum(np.einsum("ij,ik,kj->j", a, w1, a), np.finfo(float).tiny))
        b = b / np.sqrt(np.maximum(np.einsum("ij,ik,kj->j", b, w2, b), np.finfo(float).tiny))
        a, b = _sign_align(a), _sign_align(b)

        self.n_components_ = keep
        self.correlations_ = spectrum
        self.weights1_ = FData(a, basis1)
        self.weights2_ = FData(b, basis2)
        self.scores1_ = cx.T @ w1 @ a
        self.scores2_ = cy.T @ w2 @ b
        self._maps = (w1, w2)
        self._means = (
            np.asarray(to_numpy(X.mean().coefs), dtype=np.float64).reshape(-1, 1)
            if self.center
            else np.zeros((p1, 1)),
            np.asarray(to_numpy(y.mean().coefs), dtype=np.float64).reshape(-1, 1)
            if self.center
            else np.zeros((p2, 1)),
        )
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
    def correlations(self) -> NDArray:
        """The full canonical correlation spectrum, descending."""
        check_is_fitted(self)
        return self.correlations_

    @property
    def scores1(self) -> NDArray:
        """Canonical variate scores of the first sample."""
        check_is_fitted(self)
        return cast("NDArray", self.scores1_)

    @property
    def scores2(self) -> NDArray:
        """Canonical variate scores of the second sample."""
        check_is_fitted(self)
        return cast("NDArray", self.scores2_)

    def transform(self, X: FData, y: FData) -> tuple[NDArray, NDArray]:  # noqa: N803
        """Return the canonical variate scores of two new samples.

        Parameters
        ----------
        X, y : FData
            Curves to project, in the order given to :meth:`fit`.

        Returns
        -------
        tuple of numpy.ndarray
            Scores of shape ``(n_curves, n)`` for each sample.
        """
        check_is_fitted(self)
        w1, w2 = self._maps
        m1, m2 = self._means
        cx = np.asarray(to_numpy(X.coefs), dtype=np.float64).reshape(w1.shape[0], -1) - m1
        cy = np.asarray(to_numpy(y.coefs), dtype=np.float64).reshape(w2.shape[0], -1) - m2
        a = np.asarray(to_numpy(self.weights1_.coefs), dtype=np.float64)
        b = np.asarray(to_numpy(self.weights2_.coefs), dtype=np.float64)
        return cx.T @ w1 @ a, cy.T @ w2 @ b
