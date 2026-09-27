r"""Functional linear regression: :func:`fregress` and the :class:`FRegress` estimator.

One entry point replaces R ``fda``'s ``fRegress`` family (``fRegress``,
``fRegress.formula``, ``fRegress.fd``, ``fRegress.double``, ``predict.fRegress``,
``fRegress.stderr`` and ``fRegress.CV``).  The kind of model is read off the
arguments:

* **scalar response** (``y`` a vector of ``n`` numbers).  Each covariate enters
  as a scalar ``z_ij β_j`` or, when it is a curve, as ``∫ x_ij(t) β_j(t) dt``:

  .. math:: y_i = \sum_j \int x_{ij}(t)\,\beta_j(t)\,dt + \varepsilon_i .

* **functional response** (``y`` an :class:`~fabel.core.FData` of ``n``
  curves).  Scalar covariates give the function-on-scalar model and curves give
  the concurrent model; the two mix freely:

  .. math:: y_i(t) = \sum_j x_{ij}(t)\,\beta_j(t) + \varepsilon_i(t) .

Writing every coefficient function in its own basis, ``β_j = θ_jᵀ b_j``, both
models are penalised least-squares problems (Ramsay, Hooker & Graves 2009,
chapters 9-10) whose normal equations ``C b = D`` stack one block per term:

.. math::

    C_{jk} = \sum_i w_i \int x_{ij} x_{ik}\,\theta_j \theta_k^{T}\,dt
             + \delta_{jk}\,\lambda_j R_j ,
    \qquad
    D_j = \sum_i w_i \int x_{ij}\,y_i\,\theta_j\,dt

for a functional response, and ``C = ZᵀWZ + diag(λ_j R_j)``, ``D = ZᵀWy`` with
``Z_i = (∫ x_ij θ_j)_j`` for a scalar one.  ``R_j = ∫ (Lθ_j)(Lθ_j)ᵀ`` is the
roughness penalty of the ``j``-th coefficient.  Every integral is computed
exactly (Gauss-Legendre on the break points of the bases involved), whereas R
approximates them -- see ``tests/parity/test_regression.py`` for the measured
consequences.

Examples
--------
>>> import numpy as np
>>> from fabel.regression import fregress
>>> rng = np.random.default_rng(0)
>>> z = rng.standard_normal(40)
>>> y = 1.0 + 2.0 * z + 0.01 * rng.standard_normal(40)
>>> model = fregress(y, [1.0, z])
>>> np.round([float(b.coefs[0, 0]) for b in model.beta], 2).tolist()
[1.0, 2.0]
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from math import isfinite
from types import ModuleType
from typing import Any

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted, column_or_1d, validate_data

from fabel._backend import array_namespace, asarray, default_namespace
from fabel._operator import LDO
from fabel.basis import Basis, Constant, _same_domain
from fabel.core import FData, _project, _quadrature, inprod
from fabel.smoothing import SmoothResult

__all__ = [
    "FRegress",
    "FRegressCV",
    "FRegressResult",
    "FRegressStderr",
    "fregress",
]

Array = Any

#: Name R's formula interface gives the intercept term; Fabel keeps it.
INTERCEPT = "const"

_FORMULA = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_.]*)\s*~(.*)$")
_TERM = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")


# --------------------------------------------------------------------------- #
# terms
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Term:
    """One covariate together with the basis and penalty of its coefficient.

    ``values`` is an :class:`FData` of ``n`` curves for a functional covariate
    and a length-``n`` vector for a scalar one.
    """

    name: str
    values: Any
    basis: Basis
    lam: float
    penalty: LDO

    @property
    def functional(self) -> bool:
        """Return whether the covariate is a curve rather than a number."""
        return isinstance(self.values, FData)


def _as_operator(penalty: int | LDO) -> LDO:
    """Coerce an integer derivative order to an :class:`LDO`."""
    return penalty if isinstance(penalty, LDO) else LDO(int(penalty))


def _is_number(value: Any) -> bool:
    """Return whether ``value`` is a plain real number (not an array)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _scalar_vector(value: Any, n: int, xp: ModuleType, name: str) -> Array:
    """Return a scalar covariate as a length-``n`` vector, broadcasting a number."""
    if _is_number(value):
        return float(value) * xp.ones((n,), dtype=xp.float64)
    vector = asarray(value, xp)
    if len(vector.shape) == 2 and 1 in vector.shape:
        vector = xp.reshape(vector, (-1,))
    if len(vector.shape) != 1 or vector.shape[0] != n:
        raise ValueError(
            f"covariate {name!r} must be a number or hold {n} values, "
            f"got shape {tuple(vector.shape)}"
        )
    return vector


def _covariate_items(x: Any) -> list[tuple[str, Any]]:
    """Return ``(name, covariate)`` pairs for a mapping, sequence or single covariate."""
    if isinstance(x, Mapping):
        return [(str(key), value) for key, value in x.items()]
    if isinstance(x, (list, tuple)):
        return [(f"x{j}", value) for j, value in enumerate(x)]
    return [("x0", x)]


def _beta_spec(beta: Any, name: str, source: str, position: int) -> Any:
    """Look up the coefficient specification of one term.

    A mapping is searched by the term's own name first and then by the variable
    it was expanded from (so ``{"region": basis}`` covers every region dummy); a
    sequence is read by position.
    """
    if beta is None:
        return None
    if isinstance(beta, Mapping):
        if name in beta:
            return beta[name]
        return beta.get(source)
    if isinstance(beta, (list, tuple)):
        if len(beta) <= position:
            raise ValueError(f"beta has {len(beta)} entries but there are more covariates")
        return beta[position]
    return beta


def _resolve_spec(
    spec: Any, default_basis: Basis, lam: float, penalty: int | LDO
) -> tuple[Basis, float, LDO]:
    """Unpack ``None``, a basis, ``(basis, lam)`` or ``(basis, lam, penalty)``."""
    if spec is None:
        return default_basis, float(lam), _as_operator(penalty)
    if isinstance(spec, Basis):
        return spec, float(lam), _as_operator(penalty)
    if isinstance(spec, tuple) and 1 <= len(spec) <= 3 and isinstance(spec[0], Basis):
        chosen_lam = float(spec[1]) if len(spec) > 1 else float(lam)
        chosen_pen = spec[2] if len(spec) > 2 else penalty
        return spec[0], chosen_lam, _as_operator(chosen_pen)
    raise TypeError(
        f"a beta entry must be None, a Basis, (Basis, lam) or (Basis, lam, penalty); got {spec!r}"
    )


# --------------------------------------------------------------------------- #
# formula interface
# --------------------------------------------------------------------------- #


def _is_categorical(value: Any) -> bool:
    """Return whether ``value`` is a sequence of labels (strings)."""
    if isinstance(value, (str, FData)) or _is_number(value):
        return False
    try:
        entries = list(value)
    except TypeError:
        return False
    return bool(entries) and all(isinstance(entry, str) for entry in entries)


def _parse_formula(formula: str) -> tuple[str, list[str], bool]:
    """Split ``"y ~ a + b - 1"`` into the response, the terms and the intercept flag."""
    match = _FORMULA.match(formula)
    if match is None:
        raise ValueError(f"formula must look like 'response ~ term + term', got {formula!r}")
    response, rhs = match.group(1), match.group(2)
    intercept = True
    terms: list[str] = []
    for sign, word in re.findall(r"([+-]?)\s*([^+\-\s]+)", rhs):
        if word in {"0", "1"}:
            intercept = word == "1" and sign != "-"
            continue
        if sign == "-" or _TERM.match(word) is None:
            raise ValueError(f"unsupported formula term {sign}{word!r} in {formula!r}")
        if word not in terms:
            terms.append(word)
    if not terms and not intercept:
        raise ValueError(f"formula {formula!r} has no terms")
    return response, terms, intercept


def _expand_formula(
    formula: str, data: Mapping[str, Any]
) -> tuple[Any, list[tuple[str, str, Any]]]:
    """Return the response and ``(name, source, covariate)`` triples of a formula.

    The intercept is the constant ``1`` named ``const``.  A covariate given as a
    sequence of labels becomes treatment-coded indicators ``<name>.<level>``
    over its sorted levels, dropping the first level when there is an intercept
    -- the coding of R's ``fRegress.formula``.
    """
    response, terms, intercept = _parse_formula(formula)
    missing = [name for name in (response, *terms) if name not in data]
    if missing:
        raise KeyError(f"formula variables not found in data: {missing}")
    items: list[tuple[str, str, Any]] = []
    if intercept:
        items.append((INTERCEPT, INTERCEPT, 1.0))
    for name in terms:
        value = data[name]
        if not _is_categorical(value):
            items.append((name, name, value))
            continue
        labels = [str(entry) for entry in value]
        levels = sorted(set(labels))
        kept = levels[1:] if intercept else levels
        items.extend(
            (f"{name}.{level}", name, [1.0 if label == level else 0.0 for label in labels])
            for level in kept
        )
    return data[response], items


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class FRegressStderr:
    """Sampling variability of a fitted functional regression.

    Replaces the output of R's ``fRegress.stderr``.

    Attributes
    ----------
    beta : tuple of FData
        Pointwise standard error ``sqrt(θ_j(t)ᵀ V_jj θ_j(t))`` of each
        coefficient function, projected (in L2) onto the coefficient's basis.
    fitted : array or FData
        Standard error of the fitted values: a vector for a scalar response,
        curves in the response basis for a functional one.
    cov : array
        Covariance ``V`` of the stacked coefficient vector ``b`` (R's ``bvar``).

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.regression import fregress
    >>> z = np.linspace(-1.0, 1.0, 20)
    >>> y = 3.0 * z + np.sin(7.0 * z)
    >>> se = fregress(y, [1.0, z]).stderr()
    >>> se.cov.shape
    (2, 2)
    """

    beta: tuple[FData, ...]
    fitted: Any
    cov: Array


@dataclass(frozen=True)
class FRegressCV:
    """Leave-one-out cross-validation of a functional regression.

    Replaces the output of R's ``fRegress.CV``.

    Attributes
    ----------
    sse : float
        Sum over observations of the squared leave-one-out error (integrated
        over the domain for a functional response).
    errors : array or FData
        ``y_i - ŷ_i^{(-i)}``, the error of the model fitted without observation
        ``i``: a vector for a scalar response, curves in the response basis for a
        functional one.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.regression import fregress
    >>> z = np.linspace(-1.0, 1.0, 20)
    >>> cv = fregress(2.0 * z, [1.0, z]).cv()
    >>> bool(cv.sse < 1e-20)
    True
    """

    sse: float
    errors: Any


@dataclass(frozen=True, eq=False)
class FRegressResult:
    """A fitted functional linear model, as returned by :func:`fregress`.

    Attributes
    ----------
    beta : tuple of FData
        Estimated coefficient function of every term, one curve each (R's
        ``betaestlist``).  The coefficient of a scalar term in a scalar-response
        model lives in a :class:`~fabel.basis.Constant` basis.
    fitted : array or FData
        Fitted values ``ŷ``: a vector for a scalar response, curves expressed in
        the response basis for a functional one (R's ``yhatfdobj``).
    names : tuple of str
        Term names, in model order.
    y : array or FData
        The response.
    terms : tuple
        The covariates with their coefficient bases, smoothing parameters and
        penalty operators, as used in the fit.
    weights : array
        Observation weights.
    cmat, dmat : array
        The normal equations ``C b = D`` (R's ``Cmat`` and ``Dmat``).
    df : float or None
        Degrees of freedom ``tr H`` of the hat matrix (scalar response only).
    gcv : float or None
        ``SSE / (n - df)²`` (scalar response only).
    ocv : float or None
        Ordinary cross-validation score, the sum of squared leave-one-out
        residuals (scalar response only).
    y2c_map : array or None
        Data-to-coefficient map of the response smooth, kept when ``y`` was a
        :class:`~fabel.smoothing.SmoothResult`; used by :meth:`stderr`.

    Examples
    --------
    >>> import numpy as np
    >>> import fabel as fb
    >>> from fabel.regression import fregress
    >>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=6)
    >>> rng = np.random.default_rng(1)
    >>> y = fb.FData(rng.standard_normal((6, 12)), basis)
    >>> group = np.repeat([0.0, 1.0], 6)
    >>> model = fregress(y, {"const": 1.0, "group": group})
    >>> model.names
    ('const', 'group')
    >>> model.fitted.n_curves
    12
    """

    beta: tuple[FData, ...]
    fitted: Any
    names: tuple[str, ...]
    y: Any
    terms: tuple[_Term, ...]
    weights: Array
    cmat: Array
    dmat: Array
    df: float | None
    gcv: float | None
    ocv: float | None
    y2c_map: Array | None = None

    # ----------------------------------------------------------- properties

    @property
    def functional_response(self) -> bool:
        """Whether the response is a set of curves.

        Returns
        -------
        bool
            ``True`` for a functional response, ``False`` for a scalar one.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import fregress
        >>> fregress(np.arange(5.0), [1.0]).functional_response
        False
        """
        return isinstance(self.y, FData)

    @property
    def coefficients(self) -> Array:
        """The stacked coefficient vector ``b`` solving ``C b = D``.

        Returns
        -------
        array
            Vector of length ``sum_j n_basis_j``.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import fregress
        >>> fregress(np.array([1.0, 2.0, 3.0]), [1.0]).coefficients.tolist()
        [2.0]
        """
        xp = array_namespace(self.cmat)
        return xp.concat([xp.reshape(b.coefs, (-1,)) for b in self.beta])

    # -------------------------------------------------------------- methods

    def predict(self, x: Any = None) -> Any:
        """Predict the response for new covariates (R's ``predict.fRegress``).

        Parameters
        ----------
        x : mapping, sequence or covariate, optional
            New covariates, laid out like those passed to :func:`fregress`: a
            mapping keyed by term name, or a sequence in term order.  Numbers
            are broadcast.  ``None`` returns the fitted values.

        Returns
        -------
        array or FData
            Predictions: a vector for a scalar response, curves in the response
            basis for a functional one.

        Raises
        ------
        ValueError
            If ``x`` does not supply one covariate per term, or a covariate has
            the wrong kind or size.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import fregress
        >>> z = np.linspace(0.0, 1.0, 10)
        >>> model = fregress(1.0 + 2.0 * z, [1.0, z])
        >>> np.round(model.predict([1.0, np.array([0.0, 3.0])]), 8).tolist()
        [1.0, 7.0]
        """
        if x is None:
            return self.fitted
        values = self._new_covariates(x)
        if self.functional_response:
            return _functional_fit(self.y.basis, self.terms, values, self.beta)
        return _scalar_design(self.terms, values) @ self.coefficients

    def stderr(self, sigma_e: Any = None, y2c_map: Any = None) -> FRegressStderr:
        """Return the standard errors of the coefficients (R's ``fRegress.stderr``).

        The coefficients are linear in the data, ``b = M y``, so their covariance
        is ``V = M Σ Mᵀ``.  For a functional response ``y`` enters through the
        coefficients of its smooth, ``c_i = S y_i`` (``S`` the ``y2c_map`` of the
        smoothing step) and the curves are independent, so
        ``V = C⁻¹ (Σ_i M_i S Σ Sᵀ M_iᵀ) C⁻¹``.

        Parameters
        ----------
        sigma_e : float or array, optional
            Residual covariance ``Σ``.  A number means ``σ² I``.  For a scalar
            response it defaults to ``SSE / (n - df) · I``; for a functional
            response it is the ``(n_t, n_t)`` covariance of the residuals at the
            observation points and must be given.
        y2c_map : array, optional
            Map from observations to response coefficients.  The identity for a
            scalar response.  For a functional response it defaults to the map
            stored when ``y`` was a :class:`~fabel.smoothing.SmoothResult`.

        Returns
        -------
        FRegressStderr
            Pointwise standard-error curves, fitted-value standard errors and the
            coefficient covariance.

        Raises
        ------
        ValueError
            If a functional response lacks ``sigma_e`` or ``y2c_map``, or the
            shapes do not match.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import fregress
        >>> z = np.linspace(-1.0, 1.0, 30)
        >>> model = fregress(z + np.cos(9.0 * z), [1.0, z])
        >>> se = model.stderr(sigma_e=0.25)
        >>> bool(np.allclose(np.diag(se.cov), 0.25 * np.diag(np.linalg.inv(model.cmat))))
        True
        """
        if self.functional_response:
            return self._functional_stderr(sigma_e, y2c_map)
        return self._scalar_stderr(sigma_e, y2c_map)

    def cv(self) -> FRegressCV:
        """Leave-one-out cross-validation (R's ``fRegress.CV``).

        Each observation is removed in turn, the model is refitted on the rest
        with the same bases and smoothing parameters, and the held-out
        observation is predicted.

        Returns
        -------
        FRegressCV
            The summed squared errors and the individual errors.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import fregress
        >>> z = np.linspace(0.0, 1.0, 12)
        >>> y = z + 0.1 * np.sin(20.0 * z)
        >>> model = fregress(y, [1.0, z])
        >>> bool(np.isclose(model.cv().sse, model.ocv))
        True
        """
        if self.functional_response:
            return self._functional_cv()
        return self._scalar_cv()

    # ------------------------------------------------------------ internals

    def _new_covariates(self, x: Any) -> list[Any]:
        """Normalise new covariates to one entry per term, in term order."""
        if isinstance(x, Mapping):
            missing = [term.name for term in self.terms if term.name not in x]
            if missing:
                raise ValueError(f"new covariates are missing terms {missing}")
            items = [x[term.name] for term in self.terms]
        else:
            items = [value for _, value in _covariate_items(x)]
        if len(items) != len(self.terms):
            raise ValueError(f"expected {len(self.terms)} covariates, got {len(items)}")
        sizes = {v.n_curves for v in items if isinstance(v, FData)} | {
            int(asarray(v).shape[0])
            for v in items
            if not isinstance(v, FData) and not _is_number(v)
        }
        if len(sizes) > 1:
            raise ValueError(f"new covariates disagree on the number of observations: {sizes}")
        n = sizes.pop() if sizes else 1
        xp = default_namespace()
        out: list[Any] = []
        for term, value in zip(self.terms, items, strict=True):
            if term.functional != isinstance(value, FData):
                kind = "a curve" if term.functional else "a number per observation"
                raise ValueError(f"covariate {term.name!r} must be {kind}")
            out.append(
                value if isinstance(value, FData) else _scalar_vector(value, n, xp, term.name)
            )
        return out

    def _scalar_stderr(self, sigma_e: Any, y2c_map: Any) -> FRegressStderr:
        xp = array_namespace(self.cmat)
        values = [term.values for term in self.terms]
        design = _scalar_design(self.terms, values)
        n = int(design.shape[0])
        cinv = xp.linalg.inv(self.cmat)
        c2b = cinv @ (xp.matrix_transpose(design) * self.weights[None, :])
        if sigma_e is None:
            resid = self.y - self.fitted
            assert self.df is not None
            sigma_e = float(xp.sum(self.weights * resid**2)) / (n - self.df)
        sigma = _covariance(sigma_e, n, xp)
        smap = xp.eye(n, dtype=xp.float64) if y2c_map is None else asarray(y2c_map, xp)
        if tuple(smap.shape) != (n, n):
            raise ValueError(f"y2c_map must be ({n}, {n}) for a scalar response")
        mapped = c2b @ smap
        cov = mapped @ sigma @ xp.matrix_transpose(mapped)
        fitted = xp.sqrt(xp.sum((design @ cov) * design, axis=1))
        return FRegressStderr(beta=_beta_stderr(self.terms, cov), fitted=fitted, cov=cov)

    def _functional_stderr(self, sigma_e: Any, y2c_map: Any) -> FRegressStderr:
        smap_in = self.y2c_map if y2c_map is None else y2c_map
        if smap_in is None:
            raise ValueError(
                "a functional response needs y2c_map: pass it, or fit from a SmoothResult"
            )
        if sigma_e is None:
            raise ValueError("a functional response needs sigma_e, the residual covariance")
        xp = array_namespace(self.cmat)
        basis = self.y.basis
        smap = asarray(smap_in, xp)
        if len(smap.shape) != 2 or smap.shape[0] != basis.n_basis:
            raise ValueError(f"y2c_map must have {basis.n_basis} rows, got {tuple(smap.shape)}")
        sigma = _covariance(sigma_e, int(smap.shape[1]), xp)
        coef_cov = smap @ sigma @ xp.matrix_transpose(smap)
        quad = _Quadrature(basis, self.terms, xp)
        values = [term.values for term in self.terms]
        blocks = quad.response_blocks(values, self.weights)
        middle = xp.zeros(self.cmat.shape, dtype=xp.float64)
        for block in blocks:
            middle = middle + block @ coef_cov @ xp.matrix_transpose(block)
        cinv = xp.linalg.inv(self.cmat)
        cov = cinv @ middle @ cinv
        fitted = quad.fitted_stderr(values, cov)
        return FRegressStderr(beta=_beta_stderr(self.terms, cov), fitted=fitted, cov=cov)

    def _scalar_cv(self) -> FRegressCV:
        xp = array_namespace(self.cmat)
        values = [term.values for term in self.terms]
        design = _scalar_design(self.terms, values)
        y = self.y
        errors = []
        for i in range(int(design.shape[0])):
            row = design[i, :]
            weight = self.weights[i]
            cmat = self.cmat - weight * (row[:, None] * row[None, :])
            dmat = self.dmat - weight * y[i] * row
            coef = xp.linalg.solve(cmat, dmat[:, None])[:, 0]
            errors.append(y[i] - xp.sum(row * coef))
        errs = xp.stack(errors)
        return FRegressCV(sse=float(xp.sum(errs**2)), errors=errs)

    def _functional_cv(self) -> FRegressCV:
        xp = array_namespace(self.cmat)
        basis = self.y.basis
        quad = _Quadrature(basis, self.terms, xp)
        values = [term.values for term in self.terms]
        ycoefs = asarray(self.y.coefs, xp)
        columns = []
        for i in range(self.y.n_curves):
            one = [v[i : i + 1] for v in values]
            cpart, dpart = quad.normal_equations(one, ycoefs[:, i : i + 1], self.weights[i : i + 1])
            coef = xp.linalg.solve(self.cmat - cpart, (self.dmat - dpart)[:, None])[:, 0]
            betas = _split(self.terms, coef, xp)
            fit = quad.fitted_coefs(one, betas)
            columns.append(ycoefs[:, i] - fit[:, 0])
        err = xp.stack(columns, axis=1)
        gram = asarray(basis.gram(), xp)
        sse = float(xp.sum(err * (gram @ err)))
        return FRegressCV(sse=sse, errors=FData(err, basis))


def _covariance(sigma_e: Any, size: int, xp: ModuleType) -> Array:
    """Return ``sigma_e`` as a ``(size, size)`` matrix, expanding a number to ``σ² I``."""
    if _is_number(sigma_e):
        return float(sigma_e) * xp.eye(size, dtype=xp.float64)
    sigma = asarray(sigma_e, xp)
    if tuple(sigma.shape) != (size, size):
        raise ValueError(f"sigma_e must be a number or a ({size}, {size}) matrix")
    return sigma


def _split(terms: Sequence[_Term], coef: Array, xp: ModuleType) -> list[Array]:
    """Split the stacked coefficient vector into one block per term."""
    out = []
    start = 0
    for term in terms:
        stop = start + term.basis.n_basis
        out.append(coef[start:stop])
        start = stop
    return out


def _beta_stderr(terms: Sequence[_Term], cov: Array) -> tuple[FData, ...]:
    """Project the pointwise standard error of every coefficient onto its basis."""
    xp = array_namespace(cov)
    out = []
    start = 0
    for term in terms:
        size = term.basis.n_basis
        block = cov[start : start + size, start : start + size]
        start += size
        nodes, weights = _quadrature(term.basis)
        theta = asarray(term.basis(nodes), xp)
        variance = xp.sum((theta @ block) * theta, axis=1)
        pointwise = xp.sqrt(xp.clip(variance, 0.0, None))
        coefs = _project(term.basis, nodes, weights, pointwise[:, None])
        out.append(FData(coefs, term.basis))
    return tuple(out)


# --------------------------------------------------------------------------- #
# scalar response
# --------------------------------------------------------------------------- #


def _scalar_design(terms: Sequence[_Term], values: Sequence[Any]) -> Array:
    """Return ``Z`` with ``Z_i = (∫ x_ij θ_j)_j`` (``z_ij`` for a scalar term)."""
    xp = default_namespace()
    blocks = []
    for term, value in zip(terms, values, strict=True):
        if isinstance(value, FData):
            blocks.append(asarray(inprod(value, term.basis), xp))
        else:
            blocks.append(asarray(value, xp)[:, None])
    return xp.concat(blocks, axis=1)


def _penalty_blocks(terms: Sequence[_Term], xp: ModuleType) -> Array:
    """Return the block-diagonal penalty ``diag(λ_j R_j)``."""
    size = sum(term.basis.n_basis for term in terms)
    out = xp.zeros((size, size), dtype=xp.float64)
    start = 0
    for term in terms:
        stop = start + term.basis.n_basis
        if term.lam != 0.0:
            block = term.lam * asarray(term.basis.penalty(term.penalty), xp)
            out[start:stop, start:stop] = block
        start = stop
    return out


def _fit_scalar(
    y: Array, terms: tuple[_Term, ...], weights: Array
) -> tuple[list[Array], Array, Array, Array, float, float, float]:
    """Solve the scalar-response normal equations and return the fit summaries."""
    xp = default_namespace()
    design = _scalar_design(terms, [term.values for term in terms])
    n = int(design.shape[0])
    weighted = xp.matrix_transpose(design) * weights[None, :]
    cmat = weighted @ design + _penalty_blocks(terms, xp)
    dmat = weighted @ y
    cinv = xp.linalg.inv(cmat)
    coef = cinv @ dmat
    fitted = design @ coef
    hat_diag = xp.sum((design @ cinv) * xp.matrix_transpose(weighted), axis=1)
    df = float(xp.sum(hat_diag))
    resid = y - fitted
    sse = float(xp.sum(weights * resid**2))
    gcv = sse / (n - df) ** 2 if n > df else float("inf")
    ocv = float(xp.sum((resid / (1.0 - hat_diag)) ** 2))
    return _split(terms, coef, xp), fitted, cmat, dmat, df, gcv, ocv


# --------------------------------------------------------------------------- #
# functional response
# --------------------------------------------------------------------------- #


class _Quadrature:
    """Quadrature shared by every integral of a functional-response model.

    The rule is composite Gauss-Legendre on the union of the break points of the
    response basis, the covariate bases and the coefficient bases, refined until
    it resolves the products ``x_j x_k θ_j θ_k`` exactly for splines and to
    rounding error for Fourier series.
    """

    def __init__(
        self,
        response: Basis,
        terms: Sequence[_Term],
        xp: ModuleType,
        values: Sequence[Any] | None = None,
    ) -> None:
        covariates = [term.values for term in terms] if values is None else values
        cov_bases = {value.basis for value in covariates if isinstance(value, FData)}
        beta_bases = {term.basis for term in terms}
        factors = [response, *cov_bases, *cov_bases, *beta_bases, *beta_bases]
        nodes, weights = _quadrature(*factors)
        self.xp = xp
        self.nodes = nodes
        self.weights = asarray(weights, xp)
        self.response = response
        self.phi = asarray(response(nodes), xp)
        self.theta = [asarray(term.basis(nodes), xp) for term in terms]
        self.terms = tuple(terms)

    def covariate(self, value: Any) -> Array:
        """Return a covariate on the nodes: ``(Q, n)`` for curves, ``(n,)`` for numbers."""
        if isinstance(value, FData):
            return asarray(value(self.nodes), self.xp)
        return asarray(value, self.xp)

    def _product(self, left: Array, right: Array, obs_weights: Array) -> Array:
        """Return ``Σ_i w_i x_ij(t) x_ik(t)`` on the nodes (a number if both are scalar)."""
        xp = self.xp
        if len(left.shape) == 1 and len(right.shape) == 1:
            return xp.sum(obs_weights * left * right)
        a = left if len(left.shape) == 2 else left[None, :]
        b = right if len(right.shape) == 2 else right[None, :]
        return xp.sum((a * b) * obs_weights[None, :], axis=1)

    def normal_equations(
        self, values: Sequence[Any], ycoefs: Array, obs_weights: Array
    ) -> tuple[Array, Array]:
        """Return ``C`` (without penalty) and ``D`` for the given observations."""
        xp = self.xp
        covs = [self.covariate(value) for value in values]
        yvals = self.phi @ ycoefs
        rows = []
        dparts = []
        for j, theta_j in enumerate(self.theta):
            row = []
            for k, theta_k in enumerate(self.theta):
                s = self._product(covs[j], covs[k], obs_weights)
                row.append(xp.matrix_transpose(theta_j) @ ((self.weights * s)[:, None] * theta_k))
            rows.append(xp.concat(row, axis=1))
            xy = self._product(covs[j], yvals, obs_weights)
            dparts.append(xp.matrix_transpose(theta_j) @ (self.weights * xy))
        return xp.concat(rows, axis=0), xp.concat(dparts)

    def response_blocks(self, values: Sequence[Any], obs_weights: Array) -> list[Array]:
        """Return ``M_i`` with ``D = Σ_i M_i c_i``, one ``(P, K_y)`` block per curve."""
        xp = self.xp
        covs = [self.covariate(value) for value in values]
        n = int(obs_weights.shape[0])
        shared = [
            xp.matrix_transpose(theta) @ (self.weights[:, None] * self.phi) for theta in self.theta
        ]
        blocks = []
        for i in range(n):
            parts = []
            for j, theta in enumerate(self.theta):
                cov = covs[j]
                if len(cov.shape) == 1:
                    parts.append(obs_weights[i] * cov[i] * shared[j])
                else:
                    scaled = (self.weights * cov[:, i])[:, None] * self.phi
                    parts.append(obs_weights[i] * (xp.matrix_transpose(theta) @ scaled))
            blocks.append(xp.concat(parts, axis=0))
        return blocks

    def _fitted_values(self, values: Sequence[Any], betas: Sequence[Array]) -> Array:
        """Return ``Σ_j x_ij(t) β_j(t)`` on the nodes, shape ``(Q, n)``."""
        xp = self.xp
        total = None
        for theta, value, coef in zip(self.theta, values, betas, strict=True):
            curve = theta @ xp.reshape(coef, (-1,))
            cov = self.covariate(value)
            part = curve[:, None] * cov if len(cov.shape) == 2 else curve[:, None] * cov[None, :]
            total = part if total is None else total + part
        assert total is not None
        return total

    def fitted_coefs(self, values: Sequence[Any], betas: Sequence[Array]) -> Array:
        """Return the response-basis coefficients of the fitted curves."""
        return _project(self.response, self.nodes, self.weights, self._fitted_values(values, betas))

    def fitted_stderr(self, values: Sequence[Any], cov: Array) -> FData:
        """Return the pointwise standard error of every fitted curve, projected."""
        xp = self.xp
        covs = [self.covariate(value) for value in values]
        n = max(int(c.shape[-1]) for c in covs)
        sizes = [theta.shape[1] for theta in self.theta]
        offsets = [sum(sizes[:j]) for j in range(len(sizes))]
        variance = xp.zeros((self.nodes.shape[0], n), dtype=xp.float64)
        for j, theta_j in enumerate(self.theta):
            for k, theta_k in enumerate(self.theta):
                block = cov[offsets[j] : offsets[j] + sizes[j], offsets[k] : offsets[k] + sizes[k]]
                g = xp.sum((theta_j @ block) * theta_k, axis=1)
                xj = covs[j] if len(covs[j].shape) == 2 else covs[j][None, :]
                xk = covs[k] if len(covs[k].shape) == 2 else covs[k][None, :]
                variance = variance + g[:, None] * (xj * xk)
        pointwise = xp.sqrt(xp.clip(variance, 0.0, None))
        return FData(_project(self.response, self.nodes, self.weights, pointwise), self.response)


def _functional_fit(
    basis: Basis, terms: tuple[_Term, ...], values: Sequence[Any], betas: Sequence[FData]
) -> FData:
    """Return ``Σ_j x_j β_j`` for new covariates, in the response basis."""
    xp = default_namespace()
    quad = _Quadrature(basis, terms, xp, values)
    return FData(quad.fitted_coefs(values, [beta.coefs for beta in betas]), basis)


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #


def _normalise_response(y: Any) -> tuple[Any, Array | None, int]:
    """Return the response, the smooth's ``y2c_map`` (if any) and the sample size."""
    xp = default_namespace()
    y2c = None
    if isinstance(y, SmoothResult):
        if y.constraint is not None:
            raise ValueError("a constrained smooth is not linear in the data")
        y2c = None if isinstance(y.y2c_map, tuple) else asarray(y.y2c_map, xp)
        y = y.fd
    if isinstance(y, FData):
        if len(y.coefs.shape) != 2:
            raise ValueError("a functional response must not carry a variable axis")
        return FData(asarray(y.coefs, xp), y.basis), y2c, y.n_curves
    vector = asarray(y, xp)
    if len(vector.shape) == 2 and vector.shape[1] == 1:
        vector = vector[:, 0]
    if len(vector.shape) != 1:
        raise ValueError(f"a scalar response must be a vector, got shape {tuple(vector.shape)}")
    return vector, None, int(vector.shape[0])


def _build_terms(
    items: Sequence[tuple[str, str, Any]],
    response: Any,
    n: int,
    beta: Any,
    lam: float,
    penalty: int | LDO,
) -> tuple[_Term, ...]:
    """Validate the covariates and attach a coefficient basis to each."""
    xp = default_namespace()
    functional_y = isinstance(response, FData)
    curves = [value for _, _, value in items if isinstance(value, FData)]
    if functional_y:
        domain = response.domain
    elif curves:
        domain = curves[0].domain
    else:
        domain = (0.0, 1.0)
    terms = []
    for position, (name, source, value) in enumerate(items):
        if isinstance(value, FData):
            if len(value.coefs.shape) != 2:
                raise ValueError(f"covariate {name!r} must not carry a variable axis")
            if value.n_curves != n:
                raise ValueError(f"covariate {name!r} has {value.n_curves} curves, expected {n}")
            if not _same_domain(value.domain, domain):
                raise ValueError(f"covariate {name!r} lives on {value.domain}, not {domain}")
            covariate: Any = FData(asarray(value.coefs, xp), value.basis)
            default = response.basis if functional_y else value.basis
        else:
            covariate = _scalar_vector(value, n, xp, name)
            default = response.basis if functional_y else Constant(domain)
        spec = _beta_spec(beta, name, source, position)
        basis, term_lam, operator = _resolve_spec(spec, default, lam, penalty)
        if not _same_domain(basis.domain, domain):
            raise ValueError(f"the basis of beta {name!r} lives on {basis.domain}, not {domain}")
        if not functional_y and not isinstance(value, FData) and not isinstance(basis, Constant):
            raise ValueError(
                f"covariate {name!r} is scalar and so is the response: its coefficient is a "
                "number and needs a Constant basis"
            )
        if term_lam < 0.0 or not isfinite(term_lam):
            raise ValueError(f"lam for {name!r} must be finite and non-negative, got {term_lam}")
        terms.append(_Term(name, covariate, basis, term_lam, operator))
    return tuple(terms)


def fregress(
    y: Any,
    x: Any = None,
    beta: Any = None,
    *,
    lam: float = 0.0,
    penalty: int | LDO = 2,
    weights: Any = None,
) -> FRegressResult:
    r"""Fit a functional linear model; the model type follows from the arguments.

    Replaces R's ``fRegress`` (including ``fRegress.formula``, ``fRegress.fd``
    and ``fRegress.double``).

    Parameters
    ----------
    y : array, FData, SmoothResult or str
        The response: ``n`` numbers (scalar response), an :class:`FData` of
        ``n`` curves (functional response), or a
        :class:`~fabel.smoothing.SmoothResult` whose ``y2c_map`` is then kept for
        :meth:`FRegressResult.stderr`.  A string is a formula
        ``"response ~ a + b"``, evaluated against the mapping ``x``: the
        intercept ``const`` is included unless the formula has ``- 1`` or
        ``+ 0``, and a covariate given as labels (strings) is expanded into
        treatment-coded indicators ``name.level``, first sorted level dropped,
        as in R.
    x : mapping, sequence or covariate
        The covariates, each either an :class:`FData` of ``n`` curves, ``n``
        numbers, or a single number broadcast to all observations (``1.0`` is
        an intercept).  A mapping names the terms; a sequence names them
        ``x0, x1, ...``.  The data mapping of a formula.
    beta : mapping or sequence, optional
        Coefficient specification per term (by name or position): ``None``, a
        :class:`~fabel.basis.Basis`, ``(basis, lam)`` or
        ``(basis, lam, penalty)``.  A mapping may also be keyed by the variable a
        formula term came from.  The default basis is the response basis for a
        functional response, and for a scalar response the covariate's own basis
        (curves) or a :class:`~fabel.basis.Constant` (numbers).
    lam : float, optional
        Smoothing parameter applied to every term that does not give its own.
        Default ``0``.
    penalty : int or LDO, optional
        Roughness operator for every term that does not give its own; an integer
        means ``D^penalty``.  Default ``2``.
    weights : array, optional
        Positive observation weights, one per observation.

    Returns
    -------
    FRegressResult
        Coefficient functions, fitted values and fit diagnostics, with
        :meth:`~FRegressResult.predict`, :meth:`~FRegressResult.stderr` and
        :meth:`~FRegressResult.cv`.

    Raises
    ------
    ValueError
        If the covariates disagree with the response in size or domain, a
        scalar-on-scalar term is given a non-constant basis, or the formula is
        malformed.
    KeyError
        If a formula names a variable missing from the data.

    Examples
    --------
    Scalar response on a functional covariate:

    >>> import numpy as np
    >>> import fabel as fb
    >>> from fabel.regression import fregress
    >>> rng = np.random.default_rng(3)
    >>> basis = fb.BSpline(domain=(0.0, 1.0), n_basis=7)
    >>> x = fb.FData(rng.standard_normal((7, 40)), basis)
    >>> y = np.asarray(fb.inprod(x, fb.FData(np.linspace(-1, 1, 7), basis)))[:, 0]
    >>> model = fregress(y, [1.0, x])
    >>> np.round(model.beta[1].coefs[:, 0], 6).tolist()
    [-1.0, -0.666667, -0.333333, 0.0, 0.333333, 0.666667, 1.0]

    The formula interface with a categorical covariate:

    >>> region = ["north", "south", "north", "east", "south", "east"]
    >>> yf = fb.FData(rng.standard_normal((7, 6)), basis)
    >>> fregress("temp ~ region", {"temp": yf, "region": region}).names
    ('const', 'region.north', 'region.south')
    """
    if isinstance(y, str):
        if not isinstance(x, Mapping):
            raise TypeError("a formula needs the data as a mapping in x")
        response_raw, items = _expand_formula(y, x)
    else:
        if x is None:
            raise ValueError("fregress needs covariates x")
        response_raw = y
        items = [(name, name, value) for name, value in _covariate_items(x)]
    if not items:
        raise ValueError("fregress needs at least one covariate")
    response, y2c, n = _normalise_response(response_raw)
    xp = default_namespace()
    obs_weights = (
        xp.ones((n,), dtype=xp.float64) if weights is None else _scalar_vector(weights, n, xp, "w")
    )
    if bool(xp.any(obs_weights <= 0.0)):
        raise ValueError("weights must be positive")
    terms = _build_terms(items, response, n, beta, lam, penalty)
    names = tuple(term.name for term in terms)
    if len(set(names)) != len(names):
        raise ValueError(f"term names must be unique, got {names}")
    if isinstance(response, FData):
        quad = _Quadrature(response.basis, terms, xp)
        values = [term.values for term in terms]
        cpart, dmat = quad.normal_equations(values, response.coefs, obs_weights)
        cmat = cpart + _penalty_blocks(terms, xp)
        coef = xp.linalg.solve(cmat, dmat[:, None])[:, 0]
        blocks = _split(terms, coef, xp)
        fitted: Any = FData(quad.fitted_coefs(values, blocks), response.basis)
        df = gcv = ocv = None
    else:
        blocks, fitted, cmat, dmat, df, gcv, ocv = _fit_scalar(response, terms, obs_weights)
    beta_fds = tuple(
        FData(block[:, None], term.basis) for term, block in zip(terms, blocks, strict=True)
    )
    return FRegressResult(
        beta=beta_fds,
        fitted=fitted,
        names=names,
        y=response,
        terms=terms,
        weights=obs_weights,
        cmat=cmat,
        dmat=dmat,
        df=df,
        gcv=gcv,
        ocv=ocv,
        y2c_map=y2c,
    )


# --------------------------------------------------------------------------- #
# scikit-learn estimator
# --------------------------------------------------------------------------- #


def _is_matrix(X: Any) -> bool:  # noqa: N803
    """Whether ``X`` is a plain covariate matrix rather than curves."""
    if isinstance(X, FData):
        return False
    return not (isinstance(X, (list, tuple)) and any(isinstance(v, FData) for v in X))


class FRegress(RegressorMixin, BaseEstimator):  # type: ignore[misc]
    """Functional linear regression as a scikit-learn regressor.

    ``X`` is either a numeric matrix of scalar covariates (one column per
    covariate), a single :class:`FData` covariate, or a list of covariates as
    accepted by :func:`fregress`.  ``y`` is a vector (scalar response) or an
    :class:`FData` (functional response).

    Parameters
    ----------
    beta : mapping or sequence, optional
        Coefficient specification per covariate, as in :func:`fregress`.  With
        ``fit_intercept`` the intercept is the first term, named ``const``.
    lam : float, optional
        Default smoothing parameter.
    penalty : int or LDO, optional
        Default roughness operator.
    fit_intercept : bool, optional
        Prepend the constant covariate ``1``.  Default ``True``.

    Attributes
    ----------
    result_ : FRegressResult
        The fitted model.
    coef_ : tuple of FData
        The estimated coefficient functions (``result_.beta``).
    n_features_in_ : int
        Number of columns seen during ``fit`` when ``X`` was a matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from fabel.regression import FRegress
    >>> rng = np.random.default_rng(0)
    >>> X = rng.standard_normal((50, 2))
    >>> y = 0.5 + X @ np.array([1.0, -2.0])
    >>> model = FRegress().fit(X, y)
    >>> np.round(model.predict(X[:2]) - y[:2], 10).tolist()
    [0.0, 0.0]
    """

    def __init__(
        self,
        beta: Any = None,
        *,
        lam: float = 0.0,
        penalty: int | LDO = 2,
        fit_intercept: bool = True,
    ) -> None:
        self.beta = beta
        self.lam = lam
        self.penalty = penalty
        self.fit_intercept = fit_intercept

    def _covariates(self, X: Any, *, reset: bool, y: Any = None) -> tuple[list[Any], Any]:  # noqa: N803
        """Return the covariate list for ``X`` and the validated response ``y``.

        A numeric matrix goes through scikit-learn's ``validate_data`` (together
        with a numeric ``y`` when one is given); curves are passed on as they are.
        """
        if not _is_matrix(X):
            covariates: list[Any] = [X] if isinstance(X, FData) else list(X)
        elif y is None or isinstance(y, (FData, SmoothResult)):
            data = validate_data(self, X, reset=reset, ensure_min_samples=1)
            covariates = [data[:, j] for j in range(data.shape[1])]
        else:
            data, y = validate_data(self, X, y, reset=reset, ensure_min_samples=1, y_numeric=True)
            covariates = [data[:, j] for j in range(data.shape[1])]
        if self.fit_intercept:
            covariates = [1.0, *covariates]
        return covariates, y

    def fit(self, X: Any, y: Any) -> FRegress:  # noqa: N803
        """Fit the model.

        Parameters
        ----------
        X : array of shape (n_samples, n_features), FData or list
            Covariates.
        y : array of shape (n_samples,), FData or SmoothResult
            Response.  A numeric ``y`` must be finite; a column vector is
            flattened with a ``DataConversionWarning``.

        Returns
        -------
        FRegress
            The fitted estimator.

        Raises
        ------
        ValueError
            If ``y`` is missing or not finite, ``X`` and ``y`` hold different
            numbers of samples, or there is only one sample.

        Examples
        --------
        >>> import numpy as np
        >>> from fabel.regression import FRegress
        >>> X = np.array([[0.0], [1.0], [2.0], [3.0]])
        >>> model = FRegress().fit(X, 1.0 + 2.0 * X[:, 0])
        >>> [round(float(b.coefs[0, 0]), 10) for b in model.coef_]
        [1.0, 2.0]
        """
        if y is None:
            raise ValueError(
                f"{type(self).__name__} requires y to be passed, but the target y is None."
            )
        if not isinstance(y, (FData, SmoothResult)) and not _is_matrix(X):
            y = column_or_1d(
                check_array(y, ensure_2d=False, dtype="numeric", input_name="y"), warn=True
            )
        covariates, y = self._covariates(X, reset=True, y=y)
        if isinstance(y, SmoothResult):
            n = y.fd.n_curves
        elif isinstance(y, FData):
            n = y.n_curves
        else:
            n = len(y)
        if n < 2:
            raise ValueError(f"{type(self).__name__} needs at least 2 samples, got n_samples = {n}")
        for value in covariates:
            n_value = value.n_curves if isinstance(value, FData) else None
            if n_value is None and not _is_number(value):
                n_value = len(value)
            if n_value is not None and n_value != n:
                raise ValueError(f"X has {n_value} samples but y has {n}")
        self.result_ = fregress(y, covariates, self.beta, lam=self.lam, penalty=self.penalty)
        self.coef_ = self.result_.beta
        return self

    def predict(self, X: Any) -> Any:  # noqa: N803
        """Predict the response for new covariates.

        Parameters
        ----------
        X : array of shape (n_samples, n_features), FData or list
            Covariates laid out as in :meth:`fit`.

        Returns
        -------
        array or FData
            Predicted response.
        """
        check_is_fitted(self)
        covariates, _ = self._covariates(X, reset=False)
        if not any(isinstance(v, FData) for v in covariates):
            n = max(len(v) for v in covariates if not _is_number(v))
            covariates = [
                float(v) * default_namespace().ones((n,)) if _is_number(v) else v
                for v in covariates
            ]
        return self.result_.predict(covariates)
