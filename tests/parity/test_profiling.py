"""Parity of ``fabel.profiling`` against R fda 6.3.0's CSTR functions.

Golden file ``tests/golden/profiling.json`` (``tools/golden_r/profiling.R``):

* ``quadset_*``        -- :func:`simpson_rule` vs ``quadset``;
* ``cstr2in_*``        -- :func:`cstr_inputs` vs ``CSTR2in``;
* ``cstr2_*``          -- :func:`cstr_model` right-hand side vs ``CSTR2``;
* ``cstr_fitls_*``     -- :meth:`ProfiledODE.residuals` vs ``CSTRfitLS``
  (data residuals, concentration equation residual, full Jacobian);
* ``cstr_fn_errors``   -- R's ``CSTRfn``/``CSTRsse``/``CSTRres`` fail in 6.3.0
  (``CSTRfitLS`` drops the temperature equation residual), recorded;
* ``cstr_inner_*``     -- :meth:`ProfiledODE.fit_states` vs a Gauss-Newton fit
  built in R from ``CSTRfitLS`` + ``CSTR2``;
* ``cstr_profile_*``   -- :meth:`ProfiledODE.fit` vs ``nls`` on that profiled
  criterion.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fabel import BSpline
from fabel.profiling import (
    CSTR_PARAMETERS,
    ProfiledODE,
    cstr_inputs,
    cstr_model,
    simpson_rule,
)

from .conftest import case_rtol, compare, golden_cases

pytestmark = pytest.mark.parity

MODULE = "profiling"
CASES = {c["name"]: c for c in golden_cases(MODULE)}
DATA = CASES["cstr_data"]["output"]
PAR_NAMES = tuple(CSTR_PARAMETERS)


def _cases(prefix: str) -> list[dict[str, Any]]:
    return [c for c in CASES.values() if c["name"].startswith(prefix)]


def _ids(cases: list[dict[str, Any]]) -> list[str]:
    return [c["name"] for c in cases]


def _basis() -> BSpline:
    spec = DATA["basis"]
    breaks = [spec["rangeval"][0], *spec["params"], spec["rangeval"][1]]
    return BSpline(domain=tuple(spec["rangeval"]), order=4, breaks=breaks)


def _problem(
    fit: list[int], lam: list[float], estimate: tuple[str, ...] = PAR_NAMES
) -> ProfiledODE:
    y = np.array(DATA["y"], dtype=float)
    for j, observed in enumerate(fit):
        if not observed:
            y[:, j] = np.nan
    return ProfiledODE(
        cstr_model("all.cool.step", estimate=estimate),
        np.array(DATA["t"]),
        y,
        _basis(),
        lam=lam,
        state_weights=[1.0 / DATA["Cwt"], 1.0 / DATA["Twt"]],
    )


def _colmajor(matrix: Any) -> np.ndarray:
    """Flatten an R matrix (stored row-major in JSON) in R's column-major order."""
    a = np.array(matrix, dtype=float)
    if a.ndim == 1:
        return a
    return a.T.ravel()


def test_data_quadrature_is_the_default_rule() -> None:
    problem = _problem([1, 1], [100.0, 100.0])
    compare(problem.nodes, DATA["quad_nodes"], 1e-12)
    compare(problem.weights, DATA["quad_weights"], 1e-12)


@pytest.mark.parametrize("case", _cases("quadset_"), ids=_ids(_cases("quadset_")))
def test_quadset(case: dict[str, Any]) -> None:
    nodes, weights = simpson_rule(case["input"]["breaks"], case["input"]["n_quad"])
    rtol = case_rtol(MODULE, case)
    compare(nodes, case["output"]["nodes"], rtol)
    compare(weights, case["output"]["weights"], rtol)


@pytest.mark.parametrize("case", _cases("cstr2in_"), ids=_ids(_cases("cstr2in_")))
def test_cstr2in(case: dict[str, Any]) -> None:
    got = cstr_inputs(case["input"]["t"], case["input"]["condition"])
    np.testing.assert_array_equal(got, np.array(case["output"]["inputs"], dtype=float))


@pytest.mark.parametrize("case", _cases("cstr2_"), ids=_ids(_cases("cstr2_")))
def test_cstr2_rhs(case: dict[str, Any]) -> None:
    inp = case["input"]
    model = cstr_model(inp["condition"], **inp["constants"])
    got = model(np.array(inp["x"]), np.array(inp["t"]), np.array(inp["theta"]))
    compare(got, case["output"]["dx"], case_rtol(MODULE, case))


@pytest.mark.parametrize("case", _cases("cstr_fitls_"), ids=_ids(_cases("cstr_fitls_")))
def test_cstr_fitls(case: dict[str, Any]) -> None:
    inp, out = case["input"], case["output"]
    rtol = case_rtol(MODULE, case)
    problem = _problem(inp["fit"], inp["lambda"])
    coef = np.array(inp["coef"], dtype=float)
    res = problem.residuals([coef[:, 0], coef[:, 1]], inp["theta"])
    n_data, q = res.data.shape[0], problem.nodes.shape[0]
    compare(res.data, _colmajor(out["Sres"]), rtol)
    # R returns only the concentration equation residual (see cstr_fn_errors)
    compare(res.equation[:q], _colmajor(out["Lres"]), rtol)
    compare(res.jac_coefs[:n_data], out["DSres"], rtol)
    compare(res.jac_coefs[n_data:], out["DLres"], rtol)


def test_cstr_fn_fails_in_r() -> None:
    out = CASES["cstr_fn_errors"]["output"]
    for name in ("CSTRfn", "CSTRfn_gradwrd_false", "CSTRsse", "CSTRres"):
        assert out[name] == "non-conformable arguments"
    assert out["Lres_columns"] == 1
    problem = _problem([1, 1], [100.0, 100.0])
    coef = np.array(DATA["coef_smooth"], dtype=float)
    res = problem.residuals([coef[:, 0], coef[:, 1]], list(CSTR_PARAMETERS.values()))
    assert res.equation.shape == (2 * problem.nodes.shape[0],)


@pytest.mark.parametrize("case", _cases("cstr_inner_"), ids=_ids(_cases("cstr_inner_")))
def test_cstr_inner_fit(case: dict[str, Any]) -> None:
    inp, out = case["input"], case["output"]
    rtol = case_rtol(MODULE, case)
    inner = _problem(inp["fit"], inp["lambda"]).fit_states(inp["theta"])
    assert inner.converged
    coef = np.array(out["coef"], dtype=float)
    compare(inner.coefs[0], coef[:, 0], rtol)
    compare(inner.coefs[1], coef[:, 1], rtol)
    compare(inner.criterion, out["criterion"], rtol)
    compare(inner.sse, out["sse"], rtol)


@pytest.mark.parametrize("case", _cases("cstr_profile_"), ids=_ids(_cases("cstr_profile_")))
def test_cstr_profile(case: dict[str, Any]) -> None:
    inp, out = case["input"], case["output"]
    rtol = case_rtol(MODULE, case)
    estimate = tuple(np.atleast_1d(inp["estimate"]).tolist())
    fit = _problem(inp["fit"], inp["lambda"], estimate).fit(inp["start"])
    assert fit.converged
    compare(fit.theta, out["theta"], rtol)
    compare(fit.inner.sse, out["sse"], rtol)
    compare(fit.cov, out["vcov"], rtol)
    coef = np.array(out["coef"], dtype=float)
    compare(fit.inner.coefs[0], coef[:, 0], rtol)
    compare(fit.inner.coefs[1], coef[:, 1], rtol)
