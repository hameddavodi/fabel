# Phase 1 review — `_backend`, `_linalg`, `_operator`, `basis`, `core`, `_plot`

Date: 2026-09-06. Reviewer: fable (read-only). Scope per task brief.
Scratch scripts: `/private/tmp/claude-501/-Users-hamed-work-fda/52508ebb-3a0a-404d-8463-dc23ca129132/scratchpad/{verify_basis.py,verify_core.py,core_r.R,probe.py}`.

## Gate results (raw)

```
.venv/bin/pytest -q            -> 545 passed, 24 xfailed in 5.59s   (exit 0)
ruff check . ; ruff format --check . ; mypy --strict src/fabel
                               -> All checks passed! / 20 files already formatted / Success: no issues found in 7 source files
coverage (unit+parity)         -> basis 96%, core 96%, _linalg 93%, _backend 93%, _operator 100%, _plot 97%, TOTAL 96%
```

## 1. `R_FDA_DEFECTS` verdict table

Method: exact rational integrals via sympy (Bernstein Gram), `scipy.integrate.quad` at 1e-13/1e-14 tolerances (Fourier / power / inprod), symbolic derivatives (monomial / power), central finite differences (polygonal, order-6 D²), and pointwise products of factor curves evaluated by **both** Fabel and R `eval.fd` on a 2001-point grid (`times.fd`, `^.fd`). R live values were regenerated with `Rscript` and match the golden files to 0.0 relative.

"rel" = max|a − b| / max|b|.

| Case | Verdict | Numbers |
|---|---|---|
| `bspline_penalty_k4_n4_dom0_1_L0` | **R wrong** (not a convention difference) | exact Bernstein Gram [0,0] = 1/7; Fabel rel 7.8e-16; R rel 6.0. R's matrix is exactly the monomial Gram (Hilbert, rel 0.0). R's own `eval.basis` on this basis returns Bernstein values (0.421875 = 0.75³ at t = 0.25) and R's own `inprod(b, b)` returns 1/7 in [0,0] — R's `eval.penalty` is inconsistent with R's own basis functions. |
| `..._L1` | **R wrong** | exact [0,0] = 9/5; Fabel rel 8.6e-16; R rel 1.33 (first row zero). |
| `..._L2` | **R wrong** | exact [0,0] = 12; Fabel rel 7.9e-16; R rel 1.0. |
| `fourier_periodmismatch_penalty_L0` | **R wrong (quadrature error)** | Fabel rel 4.4e-16; R rel 1.38e-6; R[7,7] = 0.499999308 vs exact 0.5. |
| `..._L1` | **R wrong** | Fabel rel 3.6e-16; R rel 1.38e-6 (max abs 1.1e-4). |
| `..._L2` | **R wrong** | Fabel rel 4.4e-16; R rel 1.38e-6 (max abs 1.7e-2; R[7,7] 12468.3464 vs 12468.3637). |
| `monomial_eval_n3_d2` | **R wrong** | D²tᵉ = e(e−1)tᵉ⁻²; Fabel rel 0; R rel 1.0 (R gives 0 for e = 2). |
| `monomial_eval_n4_d2` / `n5` / `n6` | **R wrong** | Fabel rel 0; R rel 0.5 / 0.33 / 0.25 (R coefficient e(e−2) confirmed: e = 5 gives 15 instead of 20). |
| `monomial_customexp_eval_d2` | **R wrong** | Fabel rel 0; R rel 0.25. |
| `power_cfg1_eval_d1` | **R wrong** | Fabel rel 0; R rel 1.0 (R zeroes e = 0.5 and e = 1, returns 4t for e = 2). |
| `power_cfg2_eval_d1` | **R wrong** | Fabel rel 0; R rel 1.0 (R returns all zeros). |
| `power_cfg2_penalty_L0` | **R wrong** | exact [0,1] = ln 6 = 1.791759469; Fabel matches to 2e-16; R NaN at [0,1],[1,0]; R's other entries rel 2e-16. |
| `polygonal_n5_eval_d1` / `n11` / `n20` | **R wrong** | R d1 output == d0 values (rel 0 / 1e-16 / 1e-16); Fabel d1 vs central FD of d0 rel 3e-11; R d1 vs FD rel ≈ 1.1. R's own `eval.penalty(polyg, 1)` uses true slopes (±4 for h = 0.25). |
| `fd_mul_bspline` | **R wrong** | truth = R `eval.fd(f1)·eval.fd(f2)` (Fabel factor evals agree with R to 6e-16). Fabel product (order 7, n = 31, multiplicity-4 knots) rel 7.5e-15. R `times.fd` (order 7, n = 13) rel 1.28e-1 (max abs 0.254, scale 1.99). Least-squares projection of truth onto R's own basis is 6.1e-2 away from R's coefficients — R is not even the projection. |
| `fd_power2` | **R wrong** | Fabel rel 1.9e-14; R rel 5.67e-5. |
| `inprod_fourier_L0_0` / `L1_1` / `L2_2` | **R wrong (Romberg)** | vs quad: Fabel 7e-16 / 4e-16 / 1e-16; R golden 4.3e-5 / 1.2e-4 / 2.0e-4. R's own `t(C) %*% eval.penalty(b, L) %*% C` matches quad to 4e-16 — R's `inprod` disagrees with R's own exact Gram. |
| `inprod_basis_bspline_x_fourier` | **R wrong (Romberg)** | Fabel rel 9.7e-16; R rel 5.7e-5. |
| `deriv_fd_bspline_order6_L2` | **R wrong** | Fabel D² (order-4 basis) vs central FD rel 3.5e-7 (h = 1e-4); vs Fabel `fd(t, 2)` 3.9e-16; vs R `eval.fd(t, f, 2)` 3.9e-16. R `deriv.fd` vs R's own `eval.fd(t, f, 2)` rel 1.29e-2 (max abs 60, scale 4655). |

All 24 `xfail(strict=True)` entries are justified. No convention difference found; in every case R's stated result is inconsistent with another R entry point on the same object. Reason strings are accurate (the two R-side numbers quoted — 12.8% / 0.254 / 1.99, 1.3% / 59.9 / 4655, 5.7e-5, ~1.4e-6, ~1.3e-4 — all reproduce).

One nit in a reason string: `power_cfg1_eval_d1` says R "returns 4t" for e = 2; observed values are 4t (0.4 at t = 0.1) — correct.

## 2. Confirmed bugs (with reproducers)

### B1 — `FData[-1]` silently returns an empty object (High)
`src/fabel/core.py:270` — `self.coefs[:, index : index + 1, ...]` with `index = -1` is `slice(-1, 0)` → shape `(n_basis, 0)`. No error, downstream ops return empty arrays.
```python
fd = FData(np.arange(24.0).reshape(8, 3), BSpline(n_basis=8)); fd[-1].coefs.shape  # (8, 0)
```
Fix: normalise `index = range(self.n_curves)[index]` (raises IndexError when out of range), then slice. Also accept `numpy.integer` (`fd[np.int64(1)]` currently raises `TypeError: 'numpy.int64' object is not iterable`): test `isinstance(index, (int, np.integer))` via `operator.index`.

### B2 — `derivative()` crashes on splines with repeated interior knots (High)
`src/fabel/basis.py:504-517` — `_derivative_map` builds `BSpline(order=order-n, breaks=self.breaks)`; `_clean_breaks` rejects multiplicity `> order-n`. Fabel's own exact product (`fd * fd` → order 7, multiplicity 4) therefore cannot be differentiated 4 or more times, and a discontinuous spline (multiplicity = order) cannot be differentiated at all.
```python
fd = FData(np.arange(24.0).reshape(8, 3), BSpline(n_basis=8))
(fd * fd).derivative(4)   # ValueError: break 0.2 repeats 4 times, which exceeds the order 3
FData(np.ones(4), BSpline(order=2, breaks=[0, .5, .5, 1])).derivative(1)  # ValueError
BSpline(order=4, breaks=[0, .5, .5, .5, 1])._derivative_map(2)           # ValueError
```
Fix: the derivative of a multiplicity-m knot in order k is multiplicity min(m, k−n) in order k−n (the delta from a jump is dropped); the rows of `operator` to drop are the ones sitting on the now-degenerate spans, which is not `operator[n:-n]` when interior spans also collapse — derive the kept rows from the reduced knot vector instead of assuming only the end rows vanish. Add a unit test that differentiates a product basis down to order 1 and compares with finite differences.

### B3 — Autograd is cut in `inprod`, `fd * fd`, `fd ** p`, `fd.std()` (High vs SPEC §5.2 "Gradients flow through … inner products")
`core.py:648-652` (`to_numpy(first.coefs)`), `core.py:558` (`to_numpy(self(nodes))`), `core.py:616`, `core.py:462`. Verified: torch coefs with `requires_grad=True` → `inprod(fd, fd).requires_grad is False`; `(fd*fd).coefs.requires_grad is False`; `fd.std().coefs` is a NumPy array.
```python
c = torch.randn(8, 2, dtype=torch.float64, requires_grad=True); fd = FData(c, BSpline(n_basis=8))
inprod(fd, fd).requires_grad   # False
```
Fix: `_cross_gram` is a constant → keep it NumPy, but do the coefficient contractions in `result_namespace(first.coefs, second.coefs)` via `as_backend(matrix, like)` **before** multiplying. For products: evaluate `self(nodes)` in the coefs' namespace (nodes converted with `asarray(nodes, xp)`), keep `_project` backend-agnostic (`solve_spd` already dispatches to torch). Also `core.py:543 float(other)` on a tensor scalar warns and drops that scalar's gradient — multiply directly when `other` is an array.

### B4 — `fd + scalar` / `fd - scalar` silently wrong when the basis cannot represent a constant (High)
`core.py:503-509` `_constant` L2-projects the constant; on `Monomial(exponents=[1,2])`, `Power` without 0, `Exponential` without rate 0 the result is a different function with no warning.
```python
m = FData(np.ones(2), Monomial(exponents=[1, 2])); s = m + 1.0
np.max(np.abs(s(g) - (m(g) + 1)))   # 1.0   (Exponential(rates=[1,2]) + 1 → 0.153)
```
Fix: check the projection residual (`values - basis(nodes) @ coefs`) against a tolerance and raise `ValueError("basis cannot represent a constant")`; or short-circuit for bases with an exact constant (B-spline: partition of unity → coefs = value; Fourier: `value*sqrt(P)` on const; Monomial/Power with exponent 0; Exponential rate 0).

### B5 — Fallback product basis is a 1–6 % approximation presented as an operator result (Medium)
`basis.py:1509-1510` — non-matching families fall back to `BSpline(order=8, n_basis=max(n1+n2, 8))`. Measured: `BSpline(10) * Fourier(9)` product rel error **1.3e-2**; `Fourier(9, period=1) * Fourier(9, period=2)` rel error **5.9e-2**. CLAUDE.md "parity is law, 1e-8" — an operator that returns 1 % error with no signal violates that. Fix: either refine the fallback until the projection residual is below tolerance (residual check in `_project`), or raise/warn for non-closed products. Same residual check fixes B4.

### B6 — `cov()` / `std()` with 3-D coefficients silently mix variables into replications (Medium)
`core.py:486-489` flattens `(n_basis, n_curves, n_vars)` to `(n_basis, n_curves*n_vars)` and divides by `n_curves-1` → returns a single `(n_basis, n_basis)` sum over variables, not a covariance. `std()` returns `(n_basis, n_vars)` which `FData` then reads as `n_curves = n_vars`.
```python
fd3 = FData(np.random.randn(8, 5, 2), BSpline(n_basis=8)); fd3.cov().coefs.shape  # (8, 8); fd3.std().n_curves  # 2
```
Fix: either raise for `n_vars > 1` (as `inprod` does) or return the `(n_basis, n_basis, n_vars, n_vars)` cross-covariance R's `var.fd` returns.

### B7 — No torch tests for anything except `fd(t)` (Medium; eval-coverage finding)
`tests/torch/` is empty; the only torch test is `tests/unit/test_core.py:88`. B3 would have been caught by a 5-line test on `inprod`/`*`/`std`. Add `tests/torch/test_autograd.py` covering `basis(t)`, `fd(t)`, `derivative`, `inprod`, `*`, `**`, `mean/center/std`, `penalty(xp=torch)`, and a numpy-coefs/torch-t mixed case.

## 3. Suggestions ranked by severity

1. (High) B1, B2, B3, B4 above.
2. (Medium) B5, B6, B7.
3. (Medium) `core.py:641` `inprod` compares domains with exact `!=`, while `Basis.__mul__` uses `_same_domain` with `_TOL`; `BSpline(domain=(0, 1+1e-12))` vs `(0, 1)` raises in `inprod` but multiplies fine. Use one helper.
4. (Low) `core.py:645-650` the `len(coefs.shape) != 2` guard runs after `_cross_gram` — move it first (wasted quadrature before the error).
5. (Low) `basis.py:240` `Basis.__call__` silently flattens any 2-D `t` — a `(n, m)` matrix of points becomes `n*m` rows with no error; accept only 1-D or `(n, 1)`.
6. (Low) `core.py:543` `float(other)` in `__mul__`: rejects a per-curve scalar vector and, for tensor scalars, emits a warning and detaches. Branch on `is_array_api_obj`.
7. (Low) `core.py:93` `_refined_spline` does a function-local `from fabel.basis import BSpline` although `fabel.basis` is already imported at module top — no cycle exists; hoist it.
8. (Low) `_linalg.cached_gram` eviction pops the oldest inserted key (FIFO) but the docstring calls it "memoise"; fine, but note it is process-global and unbounded in key size (an `LDO` with an `FData` weight is a distinct identity key per object → cache pollution on repeated `penalty(LDO(weights=[fd]))` calls). Consider excluding non-hashable-by-value operators from the cache.
9. (Low) `basis.py:1490-1494` `_product_basis` for `Fourier × Fourier` with matching period gives `n = 2*max(n1, n2) − 1`, exact — fine; but a `Fourier` whose domain is not a whole number of periods multiplied by a `BSpline` hits the fallback (B5).
10. (Nit) `_bspline_matrix` is dense in `(n_t, n_knots)` through every recursion level; fine for Phase 1, but SPEC §6 promises banded/LRU-cached evaluation — flag for the benchmark task.

## 4. Things checked and found correct

- No `numpy` import in `basis.py`, `core.py` (typing block only), `_operator.py`, `_plot.py`; torch imported only lazily in `to_torch`. Conventions satisfied.
- Immutability: all public objects `frozen=True`; every method returns a new object.
- Order-1 splines: partition of unity holds, `deriv=1` returns zeros. `deriv >= order` returns zeros. `deriv = order−1` at the right endpoint returns 0 — matches R (`eval.basis(1, b, 3)` row is 0 in R 4.6.1 / fda 6.3.0).
- Repeated interior knots: evaluation sums to 1 to 2e-16; product of two cubics on a common break yields multiplicity 4 in order 7 and reproduces the pointwise product to 7.5e-15.
- Extrapolation: B-spline/polygonal raise `ValueError`; Fourier/monomial/power/exponential evaluate (consistent with R).
- Gradient does flow through `basis(t)` w.r.t. `t` (max error vs finite differences 6.6e-11), through `fd(t)` w.r.t. both `t` and coefficients, and through `derivative()` / `mean()` — including numpy coefficients with a torch `t`.
- `LDO` with functional weights works in both `fd(t, LDO)` and `basis.penalty(LDO)`.
- `Fourier` even `n_basis` bumped to odd, matching R.
- Parity harness (`compare`) scales `atol` by matrix norm — sensible and documented.
