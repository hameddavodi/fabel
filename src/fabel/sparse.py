r"""Functional PCA for sparse, irregular longitudinal data (PACE).

When every curve is seen at only a handful of irregular times, smoothing each
curve on its own fails: a curve with four points cannot carry a ten-function
basis.  PACE -- *principal components analysis through conditional
expectation* (Yao, Müller & Wang, 2005) -- pools the data of all curves
instead.  It replaces R ``fda``'s ``smooth.sparse.mean``, ``covPACE``,
``pcaPACE`` and ``scoresPACE``.

The input is the irregular per-curve form :func:`fabel.smoothing.smooth`
accepts: a sequence ``t`` of observation times, one array per curve, and a
matching sequence ``y`` of observed values.

1. **Mean.**  One penalised least-squares fit to the pooled points
   ``(t_ij, y_ij)``:

   .. math:: \min_c \sum_{ij} (y_{ij} - \varphi(t_{ij})^{T} c)^2
             + \lambda_\mu\, c^{T} R\, c .

2. **Covariance.**  Every pair of distinct observations of the same curve gives
   a raw covariance ``r_ij r_ik`` at ``(t_ij, t_ik)``, with
   ``r_ij = y_ij - μ(t_ij)``.  A tensor-product surface
   ``G(s, t) = ψ(s)ᵀ C ψ(t)`` is fitted to them by penalised least squares,
   with the roughness penalty ``λ_G vec(C)ᵀ (W ⊗ P + P ⊗ W) vec(C)``
   (``W`` the Gram matrix, ``P`` the roughness matrix of ``ψ``).  The diagonal
   pairs ``j = k`` are left out: they carry the measurement-error variance
   ``σ²`` on top of ``G(t, t)``.

3. **Error variance.**  The left-out diagonal ``r_ij²`` is smoothed on its own
   into ``V(t)``, and ``σ²`` is the mean gap between ``V`` and ``G`` over the
   middle half of the domain, where both are best determined:

   .. math:: \sigma^2 = \frac{2}{|T|} \int_{T_1} \{V(t) - G(t, t)\}\,dt ,
             \qquad T_1 = [a + |T|/4,\; b - |T|/4] .

4. **Eigenfunctions.**  The harmonics ``ξ = θ(t)ᵀ b`` in a basis ``θ`` solve
   the penalised eigenproblem of :class:`~fabel.decomposition.FPCA` with the
   covariance operator of ``G``:

   .. math:: J C J^{T} b = \mu\,(W_\theta + \lambda R_\theta)\,b ,
             \qquad J = \int \theta\,\psi^{T} .

5. **Scores.**  The best linear unbiased predictor of the scores of curve
   ``i`` given its own observations is the conditional expectation

   .. math:: \hat\xi_i = \Lambda \Xi_i^{T} \Sigma_i^{-1} (y_i - \mu_i),
             \qquad \Sigma_i = \Xi_i \Lambda \Xi_i^{T} + \sigma^2 I ,

   where ``Ξ_i`` holds the harmonics at the times of curve ``i`` and ``Λ`` the
   retained eigenvalues.  It works with any number of points per curve.

Every integral is exact (Gauss-Legendre on the basis break points), where R
uses the Romberg approximation of ``inprod()``; see
``tests/parity/test_pace.py`` for the measured consequences.  R's
``scoresPACE`` does not compute the conditional expectation above (it
evaluates the harmonics at a single, wrongly indexed point per curve); Fabel
implements the published estimator.

References
----------
Yao, F., Müller, H.-G. and Wang, J.-L. (2005).  Functional data analysis for
sparse longitudinal data.  *Journal of the American Statistical Association*,
100, 577-590.

Examples
--------
>>> import numpy as np
>>> from fabel import BSpline
>>> from fabel.sparse import PACE
>>> rng = np.random.default_rng(0)
>>> t = [np.sort(rng.uniform(0.0, 1.0, 5)) for _ in range(200)]
>>> a = rng.normal(0.0, 1.0, 200)
>>> y = [ai * np.sin(np.pi * ti) + 0.1 * rng.normal(size=5) for ai, ti in zip(a, t)]
>>> pace = PACE(n=1, basis=BSpline(domain=(0.0, 1.0), n_basis=4)).fit(y, t=t)
>>> pace.scores.shape
(200, 1)
>>> bool(abs(np.corrcoef(pace.scores[:, 0], a)[0, 1]) > 0.95)
True
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from itertools import pairwise
from types import ModuleType
from typing import Any, SupportsIndex

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from fabel import _linalg
from fabel._backend import asarray, default_namespace, to_numpy
from fabel._operator import LDO
from fabel.basis import Basis, BSpline
from fabel.core import BiFData, FData, inprod
from fabel.decomposition import _positive_int, _positive_sum_signs

__all__ = ["PACE", "SparseCov", "sparse_cov", "sparse_mean"]

Array = Any

#: Basis size used for the mean when none is given (capped by the data).
_DEFAULT_MEAN_SIZE = 10

#: Basis size used for the covariance surface when none is given.
_DEFAULT_COV_SIZE = 6

#: Default B-spline order of an invented basis.
_DEFAULT_ORDER = 4

#: Smallest eigenvalue, relative to the largest, of a normal matrix that is
#: still treated as non-singular.
_SINGULAR_TOL = 1e-13

#: A non-positive ``σ²`` estimate is replaced by this fraction of the mean
#: squared residual, so that every ``Σ_i`` stays positive definite.
_SIGMA2_FLOOR = 1e-6

#: Gauss-Legendre nodes per panel for the ``σ²`` integral.
_QUAD_NODES = 12

#: Panels each basis interval is split into for the ``σ²`` integral.
_QUAD_SPLIT = 4


def _xp() -> ModuleType:
    """Return the NumPy array namespace every computation here runs in."""
    return default_namespace()


# --------------------------------------------------------------------------- #
# input handling
# --------------------------------------------------------------------------- #


def _curves(y: Any, t: Any) -> tuple[list[Array], list[Array]]:
    """Return ``(times, values)`` as lists of 1-D float arrays, one per curve.

    ``t=None`` means ``y`` is R's list form: one ``(n_i, 2)`` array per curve
    holding the times in column 0 and the values in column 1.
    """
    xp = _xp()
    if isinstance(y, (str, bytes)) or not hasattr(y, "__len__"):
        raise TypeError("y must be a sequence with one array per curve")
    if len(y) == 0:
        raise ValueError("y must hold at least one curve")
    times: list[Array] = []
    values: list[Array] = []
    if t is None:
        for entry in y:
            pair = asarray(entry, xp)
            if pair.ndim != 2 or pair.shape[1] != 2:
                raise ValueError(
                    "without t, each curve must be an (n_i, 2) array of (time, value) rows"
                )
            times.append(xp.asarray(pair[:, 0], copy=True))
            values.append(xp.asarray(pair[:, 1], copy=True))
    else:
        if isinstance(t, (str, bytes)) or len(t) != len(y):
            raise ValueError("t must be a sequence with one time array per curve of y")
        for points, obs in zip(t, y, strict=True):
            times.append(xp.reshape(asarray(points, xp), (-1,)))
            values.append(xp.reshape(asarray(obs, xp), (-1,)))
    for points, obs in zip(times, values, strict=True):
        if points.shape[0] != obs.shape[0]:
            raise ValueError("each curve needs as many values as observation times")
        if points.shape[0] == 0:
            raise ValueError("every curve needs at least one observation")
        if not bool(xp.all(xp.isfinite(points))) or not bool(xp.all(xp.isfinite(obs))):
            raise ValueError("observation times and values must be finite")
    return times, values


def _default_basis(times: list[Array], size: int) -> Basis:
    """Return a B-spline basis over the observed range with at most ``size`` functions."""
    xp = _xp()
    pooled = xp.concat(times)
    lower, upper = float(xp.min(pooled)), float(xp.max(pooled))
    if not upper > lower:
        raise ValueError("the observation times must span a non-degenerate interval")
    n_basis = min(size, int(xp.unique_values(pooled).shape[0]))
    return BSpline(domain=(lower, upper), n_basis=n_basis, order=min(_DEFAULT_ORDER, n_basis))


def _check_lambda(value: float, name: str) -> float:
    """Return ``value`` as a float, rejecting negative or non-finite penalties."""
    lam = float(value)
    if not lam >= 0.0 or lam == float("inf"):
        raise ValueError(f"{name} must be a finite non-negative number, got {value!r}")
    return lam


def _as_operator(penalty: int | LDO) -> LDO:
    """Return ``penalty`` as a linear differential operator."""
    return penalty if isinstance(penalty, LDO) else LDO(int(penalty))


def _solve_checked(matrix: Array, rhs: Array, what: str) -> Array:
    """Solve a symmetric normal system, raising if it is numerically singular."""
    xp = _xp()
    sym = 0.5 * (matrix + xp.matrix_transpose(matrix))
    eig = xp.linalg.eigvalsh(sym)
    top = float(xp.max(xp.abs(eig)))
    if top == 0.0 or float(xp.min(eig)) <= _SINGULAR_TOL * top:
        raise ValueError(
            f"the {what} is not determined by the data: the normal equations are "
            "singular; use a smaller basis or a positive penalty"
        )
    return asarray(_linalg.solve_spd(sym, rhs), xp)


# --------------------------------------------------------------------------- #
# mean
# --------------------------------------------------------------------------- #


def sparse_mean(
    y: Any,
    t: Any = None,
    basis: Basis | None = None,
    *,
    lam: float = 0.0,
    penalty: int | LDO = 2,
) -> FData:
    r"""Estimate the mean function of sparse, irregularly observed curves.

    Replaces R's ``smooth.sparse.mean``.  All observations are pooled into one
    penalised least-squares fit

    .. math:: (\Phi^{T}\Phi + \lambda R)\,c = \Phi^{T} y ,

    where ``Φ`` is the basis evaluated at every observation time of every curve.

    Parameters
    ----------
    y : sequence of array_like
        Observed values, one 1-D array per curve.  With ``t=None``, one
        ``(n_i, 2)`` array per curve holding ``(time, value)`` rows instead
        (R's list form).
    t : sequence of array_like, optional
        Observation times, one 1-D array per curve, matching ``y``.
    basis : Basis, optional
        Basis of the mean.  Default: a cubic B-spline over the observed range
        with ``min(10, number of distinct times)`` functions.
    lam : float, optional
        Roughness penalty.  Default ``0.0``.
    penalty : int or LDO, optional
        Roughness operator.  Default ``2`` (the second derivative), as in R.

    Returns
    -------
    FData
        The mean as a single curve (coefficients of shape ``(n_basis, 1)``).

    Raises
    ------
    ValueError
        If the input is malformed, ``lam`` is negative, or the basis is not
        determined by the pooled data.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.sparse import sparse_mean
    >>> t = [np.array([0.0, 0.4, 0.9]), np.array([0.1, 0.5, 1.0]), np.array([0.2, 0.7])]
    >>> y = [2.0 * ti + 1.0 for ti in t]
    >>> mean = sparse_mean(y, t, BSpline(domain=(0.0, 1.0), n_basis=4))
    >>> np.round(mean(np.array([0.0, 0.5, 1.0]))[:, 0], 10).tolist()
    [1.0, 2.0, 3.0]
    """
    times, values = _curves(y, t)
    lam_value = _check_lambda(lam, "lam")
    if basis is None:
        basis = _default_basis(times, _DEFAULT_MEAN_SIZE)
    return _fit_mean(times, values, basis, lam_value, _as_operator(penalty))


def _fit_mean(times: list[Array], values: list[Array], basis: Basis, lam: float, op: LDO) -> FData:
    """Return the pooled penalised least-squares mean."""
    xp = _xp()
    phi = basis(xp.concat(times))
    target = xp.concat(values)
    normal = xp.matmul(xp.matrix_transpose(phi), phi)
    if lam > 0.0:
        normal = normal + lam * asarray(basis.penalty(op, xp=xp), xp)
    rhs = xp.matmul(xp.matrix_transpose(phi), target[:, None])
    return FData(_solve_checked(normal, rhs, "mean"), basis)


# --------------------------------------------------------------------------- #
# covariance
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SparseCov:
    """Covariance estimate of sparse longitudinal data.

    Returned by :func:`sparse_cov`; replaces the list R's ``covPACE`` returns.

    Attributes
    ----------
    cov : BiFData
        The smoothed covariance surface ``G(s, t)``, fitted to the off-diagonal
        raw covariances only.
    mean : FData
        The mean function the residuals were taken from.
    sigma2 : float
        Estimated measurement-error variance ``σ²``.
    variance : FData
        Smooth ``V(t)`` of the diagonal raw covariances ``r_ij²``, i.e. the
        variance of an observation, ``G(t, t) + σ²``.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.sparse import sparse_cov
    >>> rng = np.random.default_rng(2)
    >>> t = [np.sort(rng.uniform(0.0, 1.0, 6)) for _ in range(40)]
    >>> y = [rng.normal() + 0.2 * rng.normal(size=6) for _ in t]
    >>> est = sparse_cov(y, t, basis=BSpline(domain=(0.0, 1.0), n_basis=3, order=2))
    >>> est.cov.coefs.shape
    (3, 3)
    """

    cov: BiFData
    mean: FData
    sigma2: float
    variance: FData


def _pair_normal_equations(
    times: list[Array], residuals: list[Array], basis: Basis
) -> tuple[Array, Array]:
    """Return ``(AᵀA, Aᵀz)`` of the off-diagonal raw-covariance regression.

    Row ``(i, j, k)``, ``j ≠ k``, of ``A`` is ``ψ(t_ij) ⊗ ψ(t_ik)`` and
    ``z = r_ij r_ik``.  Summed over all ordered pairs of one curve, including
    ``j = k``, the products factor as ``G_i ⊗ G_i`` and ``u_i ⊗ u_i`` with
    ``G_i = Ψ_iᵀ Ψ_i`` and ``u_i = Ψ_iᵀ r_i``; the diagonal terms are then
    subtracted.  No pair is ever formed explicitly.
    """
    xp = _xp()
    k = basis.n_basis
    normal = xp.zeros((k * k, k * k), dtype=xp.float64)
    rhs = xp.zeros((k * k,), dtype=xp.float64)
    for points, resid in zip(times, residuals, strict=True):
        if points.shape[0] < 2:
            continue
        psi = basis(points)
        gram = xp.matmul(xp.matrix_transpose(psi), psi)
        proj = xp.matmul(xp.matrix_transpose(psi), resid)
        rows = xp.reshape(psi[:, :, None] * psi[:, None, :], (psi.shape[0], k * k))
        normal = normal + _kron(gram, gram) - xp.matmul(xp.matrix_transpose(rows), rows)
        rhs = (
            rhs
            + xp.reshape(proj[:, None] * proj[None, :], (k * k,))
            - xp.matmul(xp.matrix_transpose(rows), resid * resid)
        )
    return normal, rhs


def _kron(a: Array, b: Array) -> Array:
    """Return the Kronecker product of two matrices."""
    xp = _xp()
    rows, cols = a.shape[0] * b.shape[0], a.shape[1] * b.shape[1]
    return xp.reshape(a[:, None, :, None] * b[None, :, None, :], (rows, cols))


def _fit_cov(
    times: list[Array],
    residuals: list[Array],
    basis: Basis,
    lam: float,
    op: LDO,
    gram: Array | None = None,
) -> Array:
    """Return the symmetric coefficient matrix ``C`` of the covariance surface.

    ``gram`` overrides the exact Gram matrix in the penalty (used by the parity
    tests to substitute R's ``inprod()`` approximation).
    """
    xp = _xp()
    if not any(points.shape[0] >= 2 for points in times):
        raise ValueError("the covariance needs at least one curve with two or more observations")
    k = basis.n_basis
    normal, rhs = _pair_normal_equations(times, residuals, basis)
    if lam > 0.0:
        w = asarray(basis.gram(), xp) if gram is None else asarray(gram, xp)
        rough = asarray(basis.penalty(op, xp=xp), xp)
        normal = normal + lam * (_kron(w, rough) + _kron(rough, w))
    coefs = xp.reshape(_solve_checked(normal, rhs, "covariance surface"), (k, k))
    return 0.5 * (coefs + xp.matrix_transpose(coefs))


def _diagonal_variance(
    times: list[Array], residuals: list[Array], basis: Basis, lam: float, op: LDO
) -> Array:
    """Return the coefficients of the smooth ``V(t)`` of the squared residuals."""
    return _fit_mean(times, [r * r for r in residuals], basis, lam, op).coefs


def _middle_half_quadrature(basis: Basis) -> tuple[Array, Array]:
    """Return Gauss-Legendre nodes and weights on the middle half of the domain."""
    xp = _xp()
    lower, upper = basis.domain
    quarter = 0.25 * (upper - lower)
    left, right = lower + quarter, upper - quarter
    inner = [b for b in basis._natural_breaks() if left < b < right]
    edges = [left, *inner, right]
    panels: list[float] = []
    for a, b in pairwise(edges):
        panels.extend(a + (b - a) * j / _QUAD_SPLIT for j in range(_QUAD_SPLIT))
    panels.append(right)
    nodes, weights = _linalg.composite_gauss_legendre(
        asarray(to_numpy(xp.asarray(panels, dtype=xp.float64)), xp), _QUAD_NODES
    )
    return asarray(nodes, xp), asarray(weights, xp)


def _estimate_sigma2(
    cov: Array, variance: Array, basis: Basis, residuals: list[Array], *, warn: bool = True
) -> float:
    """Return ``σ²`` from the gap between ``V(t)`` and ``G(t, t)``.

    A non-positive estimate (possible when the data carry almost no noise) is
    replaced by a small positive floor, with a :class:`RuntimeWarning` when
    ``warn`` is true (the caller uses the estimate).
    """
    xp = _xp()
    nodes, weights = _middle_half_quadrature(basis)
    psi = basis(nodes)
    diag = xp.sum(xp.matmul(psi, cov) * psi, axis=1)
    gap = xp.matmul(psi, variance)[:, 0] - diag
    lower, upper = basis.domain
    sigma2 = 2.0 * float(xp.sum(weights * gap)) / (upper - lower)
    if sigma2 > 0.0:
        return sigma2
    pooled = xp.concat(residuals)
    floor = max(_SIGMA2_FLOOR * float(xp.mean(pooled * pooled)), 1e-12)
    if not warn:
        return floor
    warnings.warn(
        f"the estimated measurement-error variance is not positive ({sigma2:.3g}); "
        f"using {floor:.3g} instead",
        RuntimeWarning,
        stacklevel=3,
    )
    return floor


def _residuals(times: list[Array], values: list[Array], mean: FData) -> list[Array]:
    """Return ``y_ij - μ(t_ij)`` for every curve."""
    return [obs - mean(points)[:, 0] for points, obs in zip(times, values, strict=True)]


def _cov_estimate(
    times: list[Array],
    values: list[Array],
    mean: FData,
    basis: Basis,
    lam: float,
    op: LDO,
    *,
    warn: bool = True,
) -> SparseCov:
    """Return the full covariance estimate for validated input."""
    if mean.n_curves != 1:
        raise ValueError(f"mean must be a single curve, got {mean.n_curves}")
    residuals = _residuals(times, values, mean)
    coefs = _fit_cov(times, residuals, basis, lam, op)
    variance = _diagonal_variance(times, residuals, basis, lam, op)
    sigma2 = _estimate_sigma2(coefs, variance, basis, residuals, warn=warn)
    return SparseCov(
        cov=BiFData(coefs, basis, basis),
        mean=mean,
        sigma2=sigma2,
        variance=FData(variance, basis),
    )


def sparse_cov(
    y: Any,
    t: Any = None,
    *,
    mean: FData | None = None,
    basis: Basis | None = None,
    lam: float = 0.0,
    penalty: int | LDO = 2,
) -> SparseCov:
    r"""Estimate the covariance surface of sparse, irregularly observed curves.

    Replaces R's ``covPACE``.  The raw covariances ``r_ij r_ik`` of every pair
    of *distinct* observations of the same curve are fitted by the tensor
    product surface ``G(s, t) = ψ(s)ᵀ C ψ(t)``, minimising

    .. math:: \sum_i \sum_{j \ne k} \{r_{ij} r_{ik} - \psi(t_{ij})^{T} C\,
              \psi(t_{ik})\}^2 + \lambda\, \mathrm{vec}(C)^{T}
              (W \otimes P + P \otimes W)\, \mathrm{vec}(C) .

    The diagonal ``j = k`` is excluded because ``E r_ij² = G(t, t) + σ²``;
    it is smoothed separately and gives the measurement-error variance
    ``σ²`` (see :mod:`fabel.sparse`).

    Parameters
    ----------
    y : sequence of array_like
        Observed values, one 1-D array per curve, or R's list form of
        ``(n_i, 2)`` ``(time, value)`` arrays when ``t`` is ``None``.
    t : sequence of array_like, optional
        Observation times, one 1-D array per curve.
    mean : FData, optional
        Mean function to subtract.  Default: :func:`sparse_mean` of the data
        with its default basis.
    basis : Basis, optional
        Basis ``ψ`` of both arguments of the surface.  Default: a cubic
        B-spline over the observed range with ``min(6, number of distinct
        times)`` functions.
    lam : float, optional
        Roughness penalty.  Default ``0.0``.
    penalty : int or LDO, optional
        Roughness operator.  Default ``2``.

    Returns
    -------
    SparseCov
        The surface, the mean, ``σ²`` and the smoothed diagonal.

    Raises
    ------
    ValueError
        If the input is malformed, no curve has two observations, or the
        surface is not determined by the data (use a smaller basis or
        ``lam > 0``).

    Warns
    -----
    RuntimeWarning
        If the ``σ²`` estimate is not positive and is replaced by a small floor.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.sparse import sparse_cov
    >>> rng = np.random.default_rng(3)
    >>> t = [np.sort(rng.uniform(0.0, 1.0, 5)) for _ in range(200)]
    >>> y = [rng.normal(0.0, 2.0) + rng.normal(0.0, 0.5, 5) for _ in t]
    >>> est = sparse_cov(y, t, basis=BSpline(domain=(0.0, 1.0), n_basis=1, order=1))
    >>> bool(abs(float(est.cov.coefs[0, 0]) - 4.0) < 1.0)
    True
    >>> bool(abs(est.sigma2 - 0.25) < 0.1)
    True
    """
    times, values = _curves(y, t)
    lam_value = _check_lambda(lam, "lam")
    if mean is None:
        mean = _fit_mean(
            times, values, _default_basis(times, _DEFAULT_MEAN_SIZE), 0.0, _as_operator(2)
        )
    if basis is None:
        basis = _default_basis(times, _DEFAULT_COV_SIZE)
    return _cov_estimate(times, values, mean, basis, lam_value, _as_operator(penalty))


# --------------------------------------------------------------------------- #
# eigenfunctions and scores
# --------------------------------------------------------------------------- #


def _harmonics(
    cov: Array,
    cov_basis: Basis,
    harm_basis: Basis,
    lam: float,
    op: LDO,
    n: int,
    cross_gram: Array | None = None,
) -> tuple[Array, Array]:
    """Return ``(values, coefs)`` of the leading ``n`` eigenfunctions of ``G``.

    ``cross_gram`` overrides the exact ``J = ∫ θ ψᵀ`` (used by the parity tests
    to substitute R's ``inprod()`` approximation).
    """
    xp = _xp()
    j = asarray(to_numpy(inprod(harm_basis, cov_basis)), xp) if cross_gram is None else cross_gram
    j = asarray(j, xp)
    left = xp.matmul(j, xp.matmul(cov, xp.matrix_transpose(j)))
    left = 0.5 * (left + xp.matrix_transpose(left))
    right = asarray(harm_basis.gram(), xp)
    if lam > 0.0:
        right = right + lam * asarray(harm_basis.penalty(op, xp=xp), xp)
    right = 0.5 * (right + xp.matrix_transpose(right))
    mu, vec = _linalg.pencil_eigh(left, right)
    values = xp.flip(asarray(mu, xp), axis=0)
    vectors = xp.flip(asarray(vec, xp), axis=1)
    keep = min(n, int(vectors.shape[1]))
    harm = vectors[:, :keep]
    return values[:keep], harm * _whitened_sum_signs(harm, right)


def _whitened_sum_signs(harm: Array, metric: Array) -> Array:
    """Return the ``±1`` per column that makes ``1ᵀ U b`` positive.

    ``U`` is the upper Cholesky factor of the metric ``W + λR = UᵀU``, so
    ``U b`` is the harmonic in orthonormal coordinates.  This is R's
    ``pcaPACE`` rule, measured on 90 harmonics of random sparse weather
    subsamples (100 % agreement; a positive plain coefficient sum agreed on
    only 76-79 %).  A zero sum keeps the sign it has.
    """
    xp = _xp()
    upper = xp.matrix_transpose(xp.linalg.cholesky(metric))
    return _positive_sum_signs(xp.matmul(upper, harm))


def _blup_scores(
    times: list[Array],
    values: list[Array],
    mean: FData,
    harmonics: FData,
    eigenvalues: Array,
    sigma2: float,
) -> Array:
    """Return the conditional-expectation scores, one row per curve."""
    xp = _xp()
    lam = xp.where(eigenvalues > 0.0, eigenvalues, 0.0)
    rows: list[Array] = []
    for points, obs in zip(times, values, strict=True):
        xi = harmonics(points)
        resid = obs - mean(points)[:, 0]
        weighted = xi * lam
        sigma = xp.matmul(weighted, xp.matrix_transpose(xi))
        sigma = sigma + sigma2 * xp.eye(points.shape[0], dtype=xp.float64)
        rows.append(xp.matmul(xp.matrix_transpose(weighted), _linalg.solve_spd(sigma, resid)))
    return xp.stack([asarray(row, xp) for row in rows], axis=0)


class PACE(TransformerMixin, BaseEstimator):  # type: ignore[misc]
    r"""Functional PCA for sparse longitudinal data by conditional expectation.

    Replaces R's ``pcaPACE`` and ``scoresPACE`` (and runs ``smooth.sparse.mean``
    and ``covPACE`` on the way); see :mod:`fabel.sparse` for the five steps.

    Parameters
    ----------
    n : int or SupportsIndex, optional
        Number of harmonics to keep.  Default ``2``.
    basis : Basis, optional
        Basis of the covariance surface.  Default: a cubic B-spline over the
        observed range with ``min(6, number of distinct times)`` functions.
    mean_basis : Basis, optional
        Basis of the mean.  Default: as ``basis`` but with up to 10 functions.
    harmonic_basis : Basis, optional
        Basis of the harmonics.  Default: the covariance basis.
    lam_mean, lam_cov, lam : float, optional
        Roughness penalties on the mean, the covariance surface and the
        harmonics.  All default to ``0.0``.
    penalty : int or LDO, optional
        Roughness operator used by all three penalties.  Default ``2``.
    sigma2 : float, optional
        Measurement-error variance to use in the scores.  ``None`` (default)
        estimates it from the data; a given value must be positive.

    Attributes
    ----------
    mean_fd : FData
        The mean function.
    cov : BiFData
        The covariance surface, fitted without the diagonal.
    sigma2_ : float
        The measurement-error variance used for the scores.
    harmonics : FData
        The ``n`` eigenfunctions, ``(W + λR)``-orthonormal.  Each is signed
        so that its coefficients in ``(W + λR)``-orthonormal coordinates
        (``U b`` with ``W + λR = UᵀU``) sum to a positive number, R's rule.
    values : array
        The ``n`` leading eigenvalues of the covariance operator.
    varprop : array
        ``values / values.sum()``, as R's ``pcaPACE`` reports it.
    scores : array
        Conditional-expectation scores of the fitted curves, ``(n_curves, n)``.

    Notes
    -----
    ``fit`` and ``transform`` take ragged input (a different number of points
    per curve), which scikit-learn's generic estimator checks cannot generate,
    so :func:`sklearn.utils.estimator_checks.check_estimator` does not apply.
    Parameters, ``clone`` and pickling follow the estimator API.  Negative
    eigenvalues of an indefinite surface estimate are kept in ``values`` but
    count as zero in the scores.

    ``σ²`` is the difference of two estimated surfaces on their diagonal, so
    with few curves it is poorly determined and can come out negative; it is
    then floored (with a :class:`RuntimeWarning`) and the scores approach the
    least-squares fit of the harmonics to each curve.  Pass ``sigma2`` when
    the measurement error is known, or use more curves or smoother surfaces.
    With ``sigma2`` given, the estimate is still kept in ``cov_estimate_`` but
    is not used, so no warning is raised.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel import BSpline
    >>> from fabel.sparse import PACE
    >>> rng = np.random.default_rng(1)
    >>> t = [np.sort(rng.uniform(0.0, 1.0, 4)) for _ in range(150)]
    >>> y = [1.0 + rng.normal() * np.cos(np.pi * ti) + 0.05 * rng.normal(size=4) for ti in t]
    >>> pace = PACE(n=2, basis=BSpline(domain=(0.0, 1.0), n_basis=4)).fit(y, t=t)
    >>> pace.harmonics.n_curves
    2
    >>> bool(pace.varprop[0] > 0.8)
    True
    """

    def __init__(
        self,
        n: SupportsIndex = 2,
        *,
        basis: Basis | None = None,
        mean_basis: Basis | None = None,
        harmonic_basis: Basis | None = None,
        lam_mean: float = 0.0,
        lam_cov: float = 0.0,
        lam: float = 0.0,
        penalty: int | LDO = 2,
        sigma2: float | None = None,
    ) -> None:
        self.n = n
        self.basis = basis
        self.mean_basis = mean_basis
        self.harmonic_basis = harmonic_basis
        self.lam_mean = lam_mean
        self.lam_cov = lam_cov
        self.lam = lam
        self.penalty = penalty
        self.sigma2 = sigma2

    def fit(self, X: Any, y: Any = None, *, t: Any = None) -> PACE:  # noqa: N803
        """Estimate mean, covariance, harmonics and scores from sparse curves.

        Parameters
        ----------
        X : sequence of array_like
            Observed values, one 1-D array per curve; or, with ``t=None``,
            R's list form of ``(n_i, 2)`` ``(time, value)`` arrays.
        y : ignored
            Present for API compatibility.
        t : sequence of array_like, optional
            Observation times, one 1-D array per curve.

        Returns
        -------
        PACE
            The fitted estimator.

        Raises
        ------
        ValueError
            If the input or a parameter is invalid, or a basis is not
            determined by the data.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.sparse import PACE
        >>> rng = np.random.default_rng(5)
        >>> t = [np.sort(rng.uniform(0.0, 10.0, 5)) for _ in range(30)]
        >>> y = [rng.normal() * ti / 10.0 + 0.1 * rng.normal(size=5) for ti in t]
        >>> PACE(n=1).fit(y, t=t).values.shape
        (1,)
        """
        n_keep = _positive_int(self.n)
        lam_mean = _check_lambda(self.lam_mean, "lam_mean")
        lam_cov = _check_lambda(self.lam_cov, "lam_cov")
        lam = _check_lambda(self.lam, "lam")
        if self.sigma2 is not None and not float(self.sigma2) > 0.0:
            raise ValueError(f"sigma2 must be positive, got {self.sigma2!r}")
        op = _as_operator(self.penalty)
        times, values = _curves(X, t)
        cov_basis = (
            self.basis if self.basis is not None else _default_basis(times, _DEFAULT_COV_SIZE)
        )
        mean_basis = (
            self.mean_basis
            if self.mean_basis is not None
            else _default_basis(times, _DEFAULT_MEAN_SIZE)
        )
        harm_basis = self.harmonic_basis if self.harmonic_basis is not None else cov_basis
        mean = _fit_mean(times, values, mean_basis, lam_mean, op)
        # A given sigma2 replaces the estimate, so a floored estimate is not worth a warning.
        estimate = _cov_estimate(
            times, values, mean, cov_basis, lam_cov, op, warn=self.sigma2 is None
        )
        eigenvalues, coefs = _harmonics(estimate.cov.coefs, cov_basis, harm_basis, lam, op, n_keep)
        xp = _xp()
        total = float(xp.sum(eigenvalues))
        self.mean_fd_ = mean
        self.cov_ = estimate.cov
        self.cov_estimate_ = estimate
        self.sigma2_ = float(self.sigma2) if self.sigma2 is not None else estimate.sigma2
        self.harmonics_ = FData(coefs, harm_basis)
        self.values_ = eigenvalues
        self.varprop_ = eigenvalues / total if total != 0.0 else xp.zeros_like(eigenvalues)
        self.n_components_ = int(eigenvalues.shape[0])
        self.n_curves_ = len(times)
        self.scores_ = _blup_scores(times, values, mean, self.harmonics_, eigenvalues, self.sigma2_)
        return self

    # ----------------------------------------------------------- accessors

    @property
    def mean_fd(self) -> FData:
        """The mean function."""
        check_is_fitted(self)
        return self.mean_fd_

    @property
    def cov(self) -> BiFData:
        """The covariance surface."""
        check_is_fitted(self)
        return self.cov_

    @property
    def harmonics(self) -> FData:
        """The retained eigenfunctions."""
        check_is_fitted(self)
        return self.harmonics_

    @property
    def values(self) -> Array:
        """The retained eigenvalues, descending."""
        check_is_fitted(self)
        return self.values_

    @property
    def varprop(self) -> Array:
        """Share of each retained eigenvalue in their sum."""
        check_is_fitted(self)
        return self.varprop_

    @property
    def scores(self) -> Array:
        """Conditional-expectation scores of the fitted curves."""
        check_is_fitted(self)
        return self.scores_

    # ---------------------------------------------------------- transforms

    def transform(self, X: Any, t: Any = None) -> Array:  # noqa: N803
        """Return the conditional-expectation scores of new sparse curves.

        Parameters
        ----------
        X : sequence of array_like
            Observed values, one 1-D array per curve (or R's list form when
            ``t`` is ``None``).
        t : sequence of array_like, optional
            Observation times, one 1-D array per curve.

        Returns
        -------
        array
            Scores of shape ``(n_curves, n)``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.sparse import PACE
        >>> rng = np.random.default_rng(6)
        >>> t = [np.sort(rng.uniform(0.0, 1.0, 5)) for _ in range(40)]
        >>> y = [rng.normal() * np.sin(np.pi * ti) + 0.1 * rng.normal(size=5) for ti in t]
        >>> pace = PACE(n=2).fit(y, t=t)
        >>> bool(np.allclose(pace.transform(y, t), pace.scores))
        True
        """
        check_is_fitted(self)
        times, values = _curves(X, t)
        return _blup_scores(
            times, values, self.mean_fd_, self.harmonics_, self.values_, self.sigma2_
        )

    def fit_transform(self, X: Any, y: Any = None, *, t: Any = None) -> Array:  # noqa: N803
        """Fit to ``X`` and return the scores of the same curves.

        Parameters
        ----------
        X : sequence of array_like
            Observed values, one array per curve.
        y : ignored
            Present for API compatibility.
        t : sequence of array_like, optional
            Observation times, one array per curve.

        Returns
        -------
        array
            Scores of shape ``(n_curves, n)``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.sparse import PACE
        >>> rng = np.random.default_rng(7)
        >>> t = [np.sort(rng.uniform(0.0, 1.0, 5)) for _ in range(40)]
        >>> y = [rng.normal() + 0.1 * rng.normal(size=5) for ti in t]
        >>> PACE(n=1).fit_transform(y, t=t).shape
        (40, 1)
        """
        return self.fit(X, y, t=t).scores_

    def inverse_transform(self, X: Any, t: Any) -> Array:  # noqa: N803
        """Evaluate the curves rebuilt from scores, ``μ(t) + Σ_k ξ_k φ_k(t)``.

        Parameters
        ----------
        X : array_like of shape (n_samples, n)
            Scores, e.g. from :meth:`transform`.
        t : array_like
            Evaluation points (one shared grid).

        Returns
        -------
        array
            Values of shape ``(n_samples, len(t))``.

        Raises
        ------
        ValueError
            If ``X`` does not have one column per harmonic.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.sparse import PACE
        >>> rng = np.random.default_rng(8)
        >>> t = [np.sort(rng.uniform(0.0, 1.0, 6)) for _ in range(50)]
        >>> y = [rng.normal() * ti + 0.05 * rng.normal(size=6) for ti in t]
        >>> pace = PACE(n=1).fit(y, t=t)
        >>> pace.inverse_transform(pace.scores, np.linspace(0.1, 0.9, 11)).shape
        (50, 11)
        """
        check_is_fitted(self)
        xp = _xp()
        scores = asarray(X, xp)
        if scores.ndim == 1:
            scores = scores[None, :]
        if scores.ndim != 2 or scores.shape[1] != self.n_components_:
            raise ValueError(
                f"X must have shape (n_samples, {self.n_components_}), got {tuple(scores.shape)}"
            )
        grid = xp.reshape(asarray(t, xp), (-1,))
        # Harmonics are (W + λR)-orthonormal, and the scores are the weights of
        # the expansion x = μ + Σ ξ_k φ_k used by the conditional expectation.
        curves = xp.matmul(scores, xp.matrix_transpose(self.harmonics_(grid)))
        return curves + self.mean_fd_(grid)[:, 0][None, :]
