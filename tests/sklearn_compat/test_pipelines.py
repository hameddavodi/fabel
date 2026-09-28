"""fdatools estimators inside scikit-learn tooling (SPEC 5.1).

Covers the SPEC pipeline (``Smoother -> FPCA -> LogisticRegression``) under
``GridSearchCV`` and, for all five estimators including the two-view ``FCCA``,
the parameter API, ``clone`` and ``pickle`` round trips.
"""

from __future__ import annotations

import inspect
import pickle
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
import pytest
from sklearn.base import clone
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.utils.validation import check_is_fitted

from fdatools.basis import BSpline
from fdatools.core import FData
from fdatools.decomposition import FCCA, FPCA
from fdatools.registration import Registrator
from fdatools.regression import FRegress
from fdatools.smoothing import Smoother, smooth

GRID = np.linspace(0.0, 1.0, 31)
BASIS = BSpline(domain=(0.0, 1.0), n_basis=9)


def labelled_curves(n: int = 60, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Two classes of noisy curves that differ in the shape of one bump."""
    rng = np.random.default_rng(seed)
    labels = np.arange(n) % 2
    amp = 1.0 + 0.2 * rng.standard_normal((n, 1))
    centre = np.where(labels == 1, 0.35, 0.65)[:, None]
    curves = amp * np.exp(-(((GRID[None, :] - centre) / 0.12) ** 2))
    return curves + 0.05 * rng.standard_normal((n, GRID.size)), labels


def bump_fdata(centres: list[float]) -> FData:
    """Smooth bumps centred at ``centres``, one curve per centre."""
    basis = BSpline(domain=(0.0, 1.0), n_basis=12)
    values = np.exp(-(((GRID[:, None] - np.asarray(centres)[None, :]) / 0.15) ** 2))
    coefs, *_ = np.linalg.lstsq(np.asarray(basis(GRID)), values, rcond=None)
    return FData(coefs, basis)


def two_views(seed: int = 5) -> tuple[FData, FData]:
    """Two correlated samples of curves for FCCA."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((BASIS.n_basis, 30))
    return FData(a, BASIS), FData(a + 0.2 * rng.standard_normal(a.shape), BASIS)


# --------------------------------------------------------------------------- #
# SPEC 5.1 pipeline
# --------------------------------------------------------------------------- #


def spec_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("smooth", Smoother(BASIS, t=GRID, lam=1e-4)),
            ("fpca", FPCA(n=5, basis=BASIS)),
            ("clf", LogisticRegression()),
        ]
    )


def test_spec_pipeline_grid_search() -> None:
    x, y = labelled_curves()
    search = GridSearchCV(spec_pipeline(), {"fpca__n": [3, 5, 8]}, cv=3).fit(x, y)
    assert search.best_params_["fpca__n"] in (3, 5, 8)
    assert search.best_score_ > 0.9
    assert list(search.cv_results_["param_fpca__n"]) == [3, 5, 8]
    x_new, y_new = labelled_curves(20, seed=1)
    assert search.score(x_new, y_new) > 0.9
    fitted_fpca = search.best_estimator_.named_steps["fpca"]
    assert fitted_fpca.scores.shape == (60, search.best_params_["fpca__n"])


def test_spec_pipeline_pickles() -> None:
    x, y = labelled_curves()
    pipe = spec_pipeline().fit(x, y)
    again = pickle.loads(pickle.dumps(pipe))
    np.testing.assert_array_equal(again.predict_proba(x), pipe.predict_proba(x))


# --------------------------------------------------------------------------- #
# parameter API, clone and pickle for all five estimators
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Case:
    """An estimator, how to fit it, how to read its output, and a parameter to set."""

    make: Callable[[], Any]
    fit: Callable[[Any], Any]
    output: Callable[[Any], Any]
    param: tuple[str, Any]


def _smoother_case() -> Case:
    x, _ = labelled_curves(12)
    return Case(
        make=lambda: Smoother(BASIS, t=GRID, lam=1e-3),
        fit=lambda est: est.fit(x),
        output=lambda est: est.transform(x),
        param=("lam", 1e-2),
    )


def _fpca_case() -> Case:
    fd, _ = two_views()
    return Case(
        make=lambda: FPCA(n=3, lam=1e-4),
        fit=lambda est: est.fit(fd),
        output=lambda est: est.transform(fd),
        param=("n", 4),
    )


def _fcca_case() -> Case:
    fd1, fd2 = two_views()
    return Case(
        make=lambda: FCCA(n=2, lam1=1e-3, lam2=1e-3),
        fit=lambda est: est.fit(fd1, fd2),
        output=lambda est: np.concatenate(est.transform(fd1, fd2), axis=1),
        param=("lam2", 1e-2),
    )


def _fregress_case() -> Case:
    rng = np.random.default_rng(7)
    x = rng.standard_normal((25, 3))
    y = 1.0 + x @ np.array([0.5, -1.0, 2.0]) + 0.1 * rng.standard_normal(25)
    return Case(
        make=lambda: FRegress(lam=0.0),
        fit=lambda est: est.fit(x, y),
        output=lambda est: est.predict(x),
        param=("fit_intercept", False),
    )


def _registrator_case() -> Case:
    fd = bump_fdata([0.4, 0.5, 0.6])
    return Case(
        make=lambda: Registrator(BSpline(domain=(0.0, 1.0), n_basis=4), lam=1e-4),
        fit=lambda est: est.fit(fd),
        output=lambda est: est.transform(fd),
        param=("criterion", "least_squares"),
    )


CASES = {
    "Smoother": _smoother_case,
    "FPCA": _fpca_case,
    "FCCA": _fcca_case,
    "FRegress": _fregress_case,
    "Registrator": _registrator_case,
}


def _same_params(left: dict[str, Any], right: dict[str, Any]) -> None:
    assert left.keys() == right.keys()
    for key, value in left.items():
        assert repr(value) == repr(right[key]), key


@pytest.fixture(params=list(CASES), ids=list(CASES))
def case(request: pytest.FixtureRequest) -> Case:
    return CASES[request.param]()


def test_get_params_matches_constructor(case: Case) -> None:
    est = case.make()
    names = [name for name in inspect.signature(type(est).__init__).parameters if name != "self"]
    assert sorted(est.get_params(deep=False)) == sorted(names)
    assert est.get_params(deep=True).keys() == est.get_params(deep=False).keys()


def test_set_params_returns_self_and_updates(case: Case) -> None:
    est = case.make()
    name, value = case.param
    assert est.set_params(**{name: value}) is est
    assert est.get_params()[name] == value
    with pytest.raises(ValueError, match="Invalid parameter"):
        est.set_params(not_a_parameter=1)


def test_clone_is_unfitted_copy(case: Case) -> None:
    est = case.fit(case.make())
    copy = clone(est)
    assert copy is not est
    _same_params(copy.get_params(), est.get_params())
    with pytest.raises(NotFittedError):
        check_is_fitted(copy)
    np.testing.assert_allclose(case.output(case.fit(copy)), case.output(est), rtol=1e-12)


def test_fit_returns_self_and_keeps_params(case: Case) -> None:
    est = case.make()
    before = est.get_params()
    assert case.fit(est) is est
    _same_params(est.get_params(), before)


def test_pickle_unfitted(case: Case) -> None:
    est = case.make()
    again = pickle.loads(pickle.dumps(est))
    _same_params(again.get_params(), est.get_params())
    with pytest.raises(NotFittedError):
        check_is_fitted(again)


def test_pickle_fitted_output_unchanged(case: Case) -> None:
    est = case.fit(case.make())
    expected = case.output(est)
    again = pickle.loads(pickle.dumps(est))
    np.testing.assert_array_equal(case.output(again), expected)
    np.testing.assert_array_equal(case.output(est), expected)


# --------------------------------------------------------------------------- #
# FCCA directly (two-view fit)
# --------------------------------------------------------------------------- #


def test_fcca_fit_transform_reproduces_training_scores() -> None:
    fd1, fd2 = two_views()
    cca = FCCA(n=2, lam1=1e-3, lam2=1e-3).fit(fd1, fd2)
    s1, s2 = cca.transform(fd1, fd2)
    np.testing.assert_allclose(s1, cca.scores1, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(s2, cca.scores2, rtol=1e-12, atol=1e-12)
    assert cca.correlations[0] > 0.9
    assert cca.get_params() == {
        "n": 2,
        "lam1": 1e-3,
        "lam2": 1e-3,
        "penalty": 2,
        "center": True,
    }


def test_fcca_refit_is_idempotent() -> None:
    fd1, fd2 = two_views()
    cca = FCCA(n=2, lam1=1e-3, lam2=1e-3)
    first = cca.fit(fd1, fd2).scores1.copy()
    np.testing.assert_array_equal(cca.fit(fd1, fd2).scores1, first)


# --------------------------------------------------------------------------- #
# estimator-specific input handling fixed for the checks
# --------------------------------------------------------------------------- #


def test_fregress_rejects_missing_or_bad_targets() -> None:
    x = np.random.default_rng(0).standard_normal((10, 2))
    with pytest.raises(ValueError, match="requires y to be passed"):
        FRegress().fit(x, None)
    y = np.ones(10)
    y[3] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        FRegress().fit(x, y)
    with pytest.raises(ValueError, match="n_samples = 1"):
        FRegress().fit(x[:1], np.ones(1))


def test_fregress_functional_covariate_validates_y() -> None:
    fd, _ = two_views()
    good = np.asarray(fd.coefs).sum(axis=0)
    with pytest.raises(ValueError, match="requires y to be passed"):
        FRegress().fit(fd, None)
    bad = good.copy()
    bad[0] = np.inf
    with pytest.raises(ValueError, match="infinity"):
        FRegress().fit(fd, bad)
    with pytest.warns(Warning, match="column-vector"):
        est = FRegress(beta=[None, (BASIS, 1e-3)]).fit(fd, good[:, None])
    assert est.predict(fd).shape == (fd.n_curves,)
    with pytest.raises(ValueError, match="n_samples = 1"):
        FRegress().fit(fd[:1], good[:1])


def test_fregress_functional_response_on_scalar_matrix() -> None:
    fd, _ = two_views()
    x = np.random.default_rng(1).standard_normal((fd.n_curves, 1))
    est = FRegress(lam=1e-3).fit(x, fd)
    assert est.n_features_in_ == 1
    assert est.predict(x).n_curves == fd.n_curves
    with pytest.raises(ValueError, match="samples"):
        FRegress().fit(x[:5], fd)


def test_fregress_smooth_result_response() -> None:
    curves, _ = labelled_curves(12)
    response = smooth(curves.T, GRID, basis=BASIS, lam=1e-4)
    x = np.random.default_rng(3).standard_normal((12, 1))
    est = FRegress(lam=1e-3).fit(x, response)
    direct = FRegress(lam=1e-3).fit(x, response.fd)
    np.testing.assert_allclose(
        np.asarray(est.predict(x).coefs), np.asarray(direct.predict(x).coefs), rtol=1e-12
    )


def test_registrator_accepts_few_coefficients() -> None:
    # Three coefficients give a quadratic B-spline curve; same shape, different size.
    rows = np.outer([1.0, 1.5, 2.0, 2.5, 3.0], [0.0, 1.0, 0.0])
    est = Registrator().fit(rows)
    basis = est.target_.basis
    assert isinstance(basis, BSpline)
    assert basis.n_basis == 3
    assert basis.order == 3
    assert est.n_iter_ >= 0
    np.testing.assert_allclose(est.transform(rows), rows, atol=1e-8)
