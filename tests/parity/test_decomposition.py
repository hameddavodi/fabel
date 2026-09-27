"""Parity of :mod:`fabel.decomposition` against golden output from R ``fda`` 6.3.0.

Cases are parametrised by ``(case, field)``, as in ``test_smoothing.py``: a
defect in R's scores must not hide agreement in its eigenvalues.

Inputs
------
The golden file records only the outputs of ``pca.fd``/``varmx.pca.fd``/
``cca.fd``.  The smoothed curves they were computed from are read from
``tests/parity/data/fd_inputs.json`` (R's own ``smooth.basis`` coefficients),
so Fabel decomposes bit-for-bit the same data as R did.

Signs
-----
An eigenvector is only defined up to sign.  Measured on all seven ``pca.fd``
cases (18 harmonics), R returns every harmonic with a positive coefficient sum;
Fabel adopts that rule, so unrotated harmonics and scores are compared as they
are.  ``varmx.pca.fd`` does not re-sign or re-order its rotated harmonics, and
``cca.fd`` follows no sign rule at all, so for those the Fabel columns are first
matched to R's by a signed permutation (varimax) or by one sign per canonical
pair (CCA, the same sign on both weight functions so the pair stays positively
correlated).  Only the arbitrary part is aligned; magnitudes are compared at
the golden ``rtol``.
"""

from __future__ import annotations

import json
from functools import cache
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from fabel import LDO, FData
from fabel.decomposition import FCCA, FPCA

from .conftest import build_basis, case_rtol, golden_cases

pytestmark = pytest.mark.parity

MODULE = "decomposition"

INPUTS = Path(__file__).resolve().parent / "data" / "fd_inputs.json"

# --------------------------------------------------------------------------- #
# R defects
# --------------------------------------------------------------------------- #

_INPROD_GRAM_REASON = (
    "R's pca.fd builds the covariance operator W V W with the Gram matrix W "
    "from inprod(), whose Romberg quadrature is only good to 4-5 digits, while "
    "it normalises with the exact Gram matrix.  inprod(basis, basis) differs "
    "from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the "
    "synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline "
    "and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on "
    "the highest-frequency pair.  Substituting R's inprod() Gram into Fabel's "
    "eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) "
    "and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic "
    "accelerator penalty error, see test_smoothing.py).  With the exact Gram "
    "the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait "
    "eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 "
    "(weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient "
    "(gait harmonics, spread over the coupled hip/knee blocks).  Fabel uses "
    "the exact Gram matrix."
)

_INPROD_SCORES_REASON = (
    "R's pca.fd scores are inprod(centred curves, harmonics), a Romberg "
    "quadrature good to 4-5 digits, not the exact C' W h.  Measured inside R: "
    "inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on "
    "weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7.  "
    "Fabel's scores are the exact inner products; the relative gaps asserted "
    "here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 "
    "(synthetic) and 6.5e-4 (gait)."
)

_VARIMAX_REASON = (
    "R's varmx() stops before the varimax criterion is stationary and returns "
    "a rotation that is not orthogonal.  On R's own 501-point harmonic values "
    "the gradient of the criterion on the rotation group (the skew part of "
    "L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 "
    "/ weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at "
    "Fabel's, and max|T'T - I| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T.  On the "
    "weather cases R's criterion is also lower (0.04072609184 against "
    "0.04072610195; 0.04092609258 against 0.04092669536).  Growth inherits the "
    "inprod() Gram error of its source pca.fd case as well."
)

_PCA_CASES = (
    "pca_fd_weather_nharm2_harmlambda0",
    "pca_fd_weather_nharm2_harmlambda10000",
    "pca_fd_weather_nharm4_harmlambda0",
    "pca_fd_weather_nharm4_harmlambda10000",
    "pca_fd_growth_hgtm_nharm3",
    "pca_fd_synthetic_bspline_n10curves",
    "pca_fd_gait_multivariate_nharm3",
)

R_FDA_DEFECTS: dict[tuple[str, str], str] = {
    ("pca_fd_synthetic_bspline_n10curves", "values"): _INPROD_GRAM_REASON,
    ("pca_fd_synthetic_bspline_n10curves", "varprop"): _INPROD_GRAM_REASON,
    ("pca_fd_synthetic_bspline_n10curves", "harmonics"): _INPROD_GRAM_REASON,
    ("pca_fd_growth_hgtm_nharm3", "harmonics"): _INPROD_GRAM_REASON,
    ("pca_fd_gait_multivariate_nharm3", "values"): _INPROD_GRAM_REASON,
    ("pca_fd_gait_multivariate_nharm3", "harmonics"): _INPROD_GRAM_REASON,
}
for _case in _PCA_CASES:
    R_FDA_DEFECTS[(_case, "scores")] = _INPROD_SCORES_REASON
    if "weather" in _case:
        R_FDA_DEFECTS[(_case, "harmonics")] = _INPROD_GRAM_REASON
for _case in (
    "varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0",
    "varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000",
    "varmx_pca_fd_growth_hgtm_nharm3",
):
    for _field in ("values", "harmonics", "scores", "varprop", "rotmat"):
        R_FDA_DEFECTS[(_case, _field)] = _VARIMAX_REASON


# --------------------------------------------------------------------------- #
# case inputs
# --------------------------------------------------------------------------- #


@cache
def _inputs() -> dict[str, Any]:
    with INPUTS.open(encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    return data


def _coefs(key: str) -> np.ndarray:
    return np.asarray(_inputs()[key], dtype=float)


def pca_inputs(name: str, case: dict[str, Any]) -> tuple[FData, float, int | LDO]:
    """Return the curves, harmonic penalty and operator of one ``pca.fd`` case."""
    inp = case["input"]
    basis = build_basis(inp["basis"])
    lam = float(inp.get("harm_lambda", 0.0))
    if "weather" in name:
        # harmfdPar = fdPar(basis, harmaccelLfd, lambda) with the 365-day period.
        return FData(_coefs("temp_fd_coefs"), basis), lam, LDO.harmonic(365.0)
    if "growth" in name:
        return FData(_coefs("hgtm_fd_coefs"), basis), lam, 2
    if "synthetic" in name:
        return FData(np.asarray(inp["coefs"], dtype=float), basis), lam, 2
    return FData(_coefs("gait_fd_coefs"), basis), lam, LDO.harmonic(20.0)


_FITS: dict[str, Any] = {}


def fit_for(name: str, case: dict[str, Any]) -> Any:
    """Run (and memoise) the Fabel fit that replays one golden case."""
    if name in _FITS:
        return _FITS[name]
    if name.startswith("pca"):
        fd, lam, operator = pca_inputs(name, case)
        _FITS[name] = FPCA(n=int(case["input"]["nharm"]), lam=lam, penalty=operator).fit(fd)
    elif name.startswith("varmx"):
        source = case["input"]["source_case"]
        _FITS[name] = fit_for(source, CASES_BY_NAME[source]).rotate("varimax")
    else:
        inp = case["input"]
        basis = build_basis(inp["basis"])
        lam = float(inp["lambda"])
        x = FData(_coefs("temp2_fd_coefs"), basis)
        y = FData(_coefs("logprecip_fd_coefs"), basis)
        _FITS[name] = FCCA(n=int(inp["ncan"]), lam1=lam, lam2=lam, penalty=2).fit(x, y)
    return _FITS[name]


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #


def compare(actual: Any, expected: Any, rtol: float) -> None:
    """Compare against R, reshaping to R's layout and scaling ``atol``."""
    want = np.asarray(expected, dtype=float)
    got = np.reshape(np.asarray(actual, dtype=float), want.shape)
    scale = float(np.max(np.abs(want))) if want.size else 1.0
    np.testing.assert_allclose(got, want, rtol=rtol, atol=1e-12 * max(1.0, scale))


def signed_permutation(fabel: np.ndarray, r: np.ndarray) -> tuple[list[int], np.ndarray]:
    """Return the column order and signs that best map ``fabel``'s columns onto ``r``'s."""
    n = fabel.shape[1]
    best: tuple[float, list[int]] = (-1.0, list(range(n)))
    for candidate in permutations(range(n)):
        score = sum(abs(float(np.dot(fabel[:, j], r[:, i]))) for i, j in enumerate(candidate))
        if score > best[0]:
            best = (score, list(candidate))
    order = best[1]
    signs = np.array([np.sign(np.dot(fabel[:, j], r[:, i])) or 1.0 for i, j in enumerate(order)])
    return order, signs


def check_pca_field(name: str, case: dict[str, Any], field: str, rtol: float) -> None:
    """Compare one field of a ``pca.fd`` or ``varmx.pca.fd`` case."""
    fit = fit_for(name, case)
    out = case["output"]
    harmonics = np.asarray(fit.harmonics.coefs, dtype=float)
    multivariate = harmonics.ndim == 3
    scores = np.asarray(fit.scores_by_var if multivariate else fit.scores, dtype=float)
    values = np.asarray(fit.values, dtype=float)
    varprop = np.asarray(fit.varprop, dtype=float)
    rotation = None if fit.rotation is None else np.asarray(fit.rotation, dtype=float)
    if name.startswith("varmx"):
        order, signs = signed_permutation(harmonics, np.asarray(out["harmonics_coefs"]))
        harmonics = harmonics[:, order] * signs
        scores = scores[:, order] * signs
        values, varprop = values[order], varprop[order]
        assert rotation is not None
        rotation = rotation[:, order] * signs
    if field == "values":
        compare(values, out["values"], rtol)
    elif field == "harmonics":
        compare(harmonics, out["harmonics_coefs"], rtol)
    elif field == "scores":
        compare(scores, out["scores"], rtol)
    elif field == "varprop":
        compare(varprop, out["varprop"], rtol)
    elif field == "meanfd":
        compare(fit.mean_fd.coefs, out["meanfd_coefs"], rtol)
    elif field == "rotmat":
        compare(rotation, out["rotmat"], rtol)
    else:
        raise AssertionError(f"unknown pca field {field!r}")


def check_cca_field(name: str, case: dict[str, Any], field: str, rtol: float) -> None:
    """Compare one field of a ``cca.fd`` case."""
    fit = fit_for(name, case)
    out = case["output"]
    w1 = np.asarray(fit.weights1.coefs, dtype=float)
    signs = np.sign(np.sum(w1 * np.asarray(out["ccawtfd1_coefs"]), axis=0))
    if field == "ccacorr":
        compare(fit.correlations, out["ccacorr"], rtol)
    elif field == "ccawtfd1":
        compare(w1 * signs, out["ccawtfd1_coefs"], rtol)
    elif field == "ccawtfd2":
        compare(np.asarray(fit.weights2.coefs) * signs, out["ccawtfd2_coefs"], rtol)
    elif field == "ccavar1":
        compare(np.asarray(fit.scores1) * signs, out["ccavar1"], rtol)
    elif field == "ccavar2":
        compare(np.asarray(fit.scores2) * signs, out["ccavar2"], rtol)
    else:
        raise AssertionError(f"unknown cca field {field!r}")


def fields_of(case: dict[str, Any]) -> tuple[str, ...]:
    """Return the comparable output fields of one golden case."""
    if case["name"].startswith("cca"):
        return ("ccacorr", "ccawtfd1", "ccawtfd2", "ccavar1", "ccavar2")
    base = ("values", "harmonics", "scores", "varprop", "meanfd")
    return (*base, "rotmat") if "rotmat" in case["output"] else base


CASES = golden_cases(MODULE)
CASES_BY_NAME = {case["name"]: case for case in CASES}
PARAMS = [(case, field) for case in CASES for field in fields_of(case)]


def test_every_golden_case_is_replayed() -> None:
    """All twelve golden cases are exercised, none silently dropped."""
    assert len(CASES) == 12
    assert {case["name"] for case, _ in PARAMS} == set(CASES_BY_NAME)


@pytest.mark.parametrize(("case", "field"), PARAMS, ids=[f"{c['name']}-{f}" for c, f in PARAMS])
def test_decomposition_matches_r_fda(
    case: dict[str, Any], field: str, request: pytest.FixtureRequest
) -> None:
    """Replay one field of one golden decomposition case against R ``fda`` 6.3.0."""
    name = case["name"]
    reason = R_FDA_DEFECTS.get((name, field))
    if reason is not None:
        request.applymarker(pytest.mark.xfail(strict=True, reason=reason))
    rtol = case_rtol(MODULE, case)
    if name.startswith("cca"):
        check_cca_field(name, case, field, rtol)
    else:
        check_pca_field(name, case, field, rtol)


def test_pca_harmonics_follow_r_sign_rule() -> None:
    """Every unrotated harmonic has a positive coefficient sum, as R reports them."""
    for name in _PCA_CASES:
        harmonics = np.asarray(fit_for(name, CASES_BY_NAME[name]).harmonics.coefs, dtype=float)
        stacked = harmonics if harmonics.ndim == 2 else np.concatenate(list(harmonics.T), axis=1).T
        want = np.asarray(CASES_BY_NAME[name]["output"]["harmonics_coefs"], dtype=float)
        want_stacked = want if want.ndim == 2 else np.concatenate(list(want.T), axis=1).T
        assert np.all(stacked.sum(axis=0) > 0.0), name
        assert np.all(want_stacked.sum(axis=0) > 0.0), name
