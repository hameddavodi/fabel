# Parity report

Fabel 1.0.0 against R `fda` 6.3.0 (R version 4.6.1 (2026-06-24)).
Generated 2026-09-27 by `tools/parity_report.py` from a live run of `tests/parity`.

Errors are measured on every `assert_allclose` a parity test makes:
absolute error is `max |fabel - R|`, relative error is the normwise
`max |fabel - R| / max |R|` (entrywise ratios are meaningless on matrix entries
that are zero up to rounding). The summary columns show
the worst value over **passing** cases only; strict xfails (cases where R is
demonstrably the less accurate side) are listed with their measured reason below.

**Totals:** 486 parity checks -- 414 pass, 72 strict xfail (R defect), 0 fail, 0 skipped.

## By public symbol

| Symbol | R counterpart | Golden cases | Checks | Max abs err | Max rel err | Status |
|---|---|---:|---:|---:|---:|---|
| `FData` | fd, eval.fd, mean.fd, sd.fd, center.fd, deriv.fd, +.fd | 17 | 16 | 3.41e-12 | 1.77e-15 | pass, 3 xfail (R defect) |
| `BiFData` | bifd, eval.bifd | 1 | 1 | 0 | 0 | pass |
| `LDO` | Lfd, int2Lfd, vec2Lfd | 3 | 3 | 4.55e-13 | 5.31e-15 | pass |
| `inprod` | inprod, inprod.bspline | 8 | 8 | 9.95e-14 | 1.63e-15 | pass, 4 xfail (R defect) |
| `Basis` | basisfd | 0 | 0 | — | — | no golden cases |
| `BSpline` | create.bspline.basis, bsplineS, bsplinepen | 129 | 129 | 6.40e-09 | 5.93e-14 | pass, 3 xfail (R defect) |
| `Fourier` | create.fourier.basis, fourier, fourierpen | 43 | 43 | 4.77e-07 | 2.84e-14 | pass, 3 xfail (R defect) |
| `Monomial` | create.monomial.basis, monomial | 42 | 42 | 7.11e-15 | 1.78e-16 | pass, 5 xfail (R defect) |
| `Exponential` | create.exponential.basis, expon | 18 | 18 | 9.09e-13 | 1.67e-16 | pass |
| `Power` | create.power.basis, powerbasis, powerpen | 8 | 8 | 0 | 0 | pass, 3 xfail (R defect) |
| `Constant` | create.constant.basis | 4 | 4 | 0 | 0 | pass |
| `Polygonal` | create.polygonal.basis, polyg, polygpen | 12 | 12 | 7.11e-15 | 2.77e-15 | pass, 3 xfail (R defect) |
| `smooth` | smooth.basis, smooth.basisPar, Data2fd, smooth.monotone, smooth.pos | 30 | 176 | 5.00e-05 | 3.02e-06 | pass, 43 xfail (R defect) |
| `SmoothResult` | smooth.basis return list | 0 | 0 | — | — | no golden cases |
| `Smoother` | (new: sklearn estimator) | 0 | 0 | — | — | no golden cases |
| `gcv_curve` | lambda2gcv | 6 | 7 | 4.98e-10 | 4.01e-10 | pass, 2 xfail (R defect) |
| `lambda_to_df` | lambda2df | 3 | 3 | 6.32e-13 | 1.20e-14 | pass |
| `df_to_lambda` | df2lambda | 3 | 3 | — | — | all 3 checks xfail (R defect) |
| `FPCA` | pca.fd, varmx.pca.fd | 10 | 0 | — | — | golden cases exist, no parity test yet |
| `FCCA` | cca.fd | 2 | 0 | — | — | golden cases exist, no parity test yet |
| `fregress` | fRegress, predict.fRegress, fRegress.stderr, fRegress.CV | 9 | 0 | — | — | not built on this branch |
| `FRegress` | fRegress (estimator form) | 0 | 0 | — | — | not built on this branch |
| `register` | register.fd, AmpPhaseDecomp | 3 | 0 | — | — | not built on this branch |
| `landmark_register` | landmarkreg | 1 | 0 | — | — | not built on this branch |
| `Registrator` | register.fd (estimator form) | 0 | 0 | — | — | not built on this branch |
| `PDA` | pda.fd, pda.overlay | 3 | 0 | — | — | not built on this branch |
| `phase_plane` | phaseplanePlot | 0 | 0 | — | — | not built on this branch |
| `stats.cov` | var.fd | 1 | 0 | — | — | not built on this branch |
| `stats.cor` | cor.fd | 1 | 0 | — | — | not built on this branch |
| `stats.depth` | fdepth | 1 | 0 | — | — | not built on this branch |
| `stats.boxplot` | fbplot, boxplot.fd | 1 | 0 | — | — | not built on this branch |
| `stats.f_test` | Fperm.fd | 1 | 0 | — | — | not built on this branch |
| `stats.t_test` | tperm.fd | 1 | 0 | — | — | not built on this branch |
| `datasets.load_*` | data(package='fda') | 11 | 11 | 3.21e-09 | 5.94e-14 | pass |
| `nn.BasisLayer` | (new: PyTorch layer) | 0 | 0 | — | — | not built on this branch |
| `nn.FDataDataset` | (new: PyTorch dataset) | 0 | 0 | — | — | not built on this branch |
| `from_pandas` | (new) | 0 | 0 | — | — | no golden cases |
| `to_pandas` | (new) | 0 | 0 | — | — | no golden cases |
| `to_xarray` | (new) | 0 | 0 | — | — | no golden cases |
| `read_rds` | readRDS on fd / bifd / basisfd objects | 0 | 2 | 2.22e-16 | 9.98e-17 | pass |

## Strict xfails (R fda is the less accurate side)

- `bspline_penalty_k4_n4_dom0_1_L0` (basis, rel err 8.57e-01): With no interior knots R returns the monomial Gram (the Hilbert matrix) instead of the Bernstein-basis Gram; Fabel's [0, 0] entry is 1/7.
- `bspline_penalty_k4_n4_dom0_1_L1` (basis, rel err 1.33e+00): Same zero-interior-knot path as L0: R's first rows are identically zero, which is impossible for a basis with nonzero first derivatives.
- `bspline_penalty_k4_n4_dom0_1_L2` (basis, rel err 3.00e+00): R returns [[0,0,0,0],[0,0,0,0],[0,0,4,6],[0,0,6,12]]. D^2 B_0 = 6(1-t), so entry [0, 0] is the integral of 36 (1-t)^2 = 12, not 0.
- `fourier_periodmismatch_penalty_L0` (basis, rel err 1.38e-06): Period 2 over domain (0, 1) leaves R on a numerical quadrature accurate only to ~1.4e-6 (its [7, 7] is 0.499999); the exact value is 1/2.
- `fourier_periodmismatch_penalty_L1` (basis, rel err 1.38e-06): Same partial-period quadrature error as L0.
- `fourier_periodmismatch_penalty_L2` (basis, rel err 1.38e-06): Same partial-period quadrature error as L0.
- `monomial_eval_n3_d2` (basis, rel err inf): R's second-derivative coefficient is prod_r (e - 2r) rather than prod_r (e - r), so D^2 t^e comes out as e (e - 2) t^(e-2).
- `monomial_eval_n4_d2` (basis, rel err 1.00e+00): Same second-derivative coefficient defect as n3_d2.
- `monomial_eval_n5_d2` (basis, rel err 5.00e-01): Same second-derivative coefficient defect as n3_d2.
- `monomial_eval_n6_d2` (basis, rel err 3.33e-01): Same second-derivative coefficient defect as n3_d2.
- `monomial_customexp_eval_d2` (basis, rel err 3.33e-01): Same second-derivative coefficient defect as n3_d2.
- `power_cfg1_eval_d1` (basis, rel err 5.00e-01): R zeroes the derivative of every power except e = 2, and returns 4t there instead of 2t.
- `power_cfg2_eval_d1` (basis, rel err inf): Same first-derivative defect as cfg1.
- `power_cfg2_penalty_L0` (basis, rel err 0): Where e_i + e_j = -1 the antiderivative is a logarithm; R evaluates 0/0 and returns NaN. The exact entry is log(3 / 0.5).
- `polygonal_n5_eval_d1` (basis, rel err 5.00e+00): R's eval.basis ignores nderiv for a polygonal basis and returns the values themselves; its own eval.penalty uses the true derivatives.
- `polygonal_n11_eval_d1` (basis, rel err 1.10e+01): Same ignored-nderiv defect as n5.
- `polygonal_n20_eval_d1` (basis, rel err 2.00e+01): Same ignored-nderiv defect as n5.
- `inprod_fourier_L0_0` (core, rel err 4.26e-05): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. Fabel uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_fourier_L1_1` (core, rel err 1.17e-04): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. Fabel uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_fourier_L2_2` (core, rel err 1.99e-04): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. Fabel uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_basis_bspline_x_fourier` (core, rel err 5.73e-05): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. Fabel uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `deriv_fd_bspline_order6_L2` (core, rel err 1.27e-02): deriv.fd re-expands D^2 x in the *original* order-6 basis, which holds only C^4 functions while D^2 x is C^2; the result is 1.3% off (max 59.9 on a curve of size 4655). Fabel returns the exact order-4 spline.
- `fd_mul_bspline` (core, rel err 1.31e-01): times.fd keeps the factors' break points, so the order-7 product space is C^5 where the true product is only C^2 and cannot be represented. R's answer is 12.8% away from the exact product (max 0.254 on a curve of size 1.99), and is not even the L2 projection onto its own basis. Fabel raises the interior knot multiplicities and is exact.
- `fd_power2` (core, rel err 5.14e-05): ^.fd projects onto an arbitrary uniform refinement (80 intervals, order 7) whose knots miss the curve's own breaks at k/7, leaving a 5.7e-5 relative error. Fabel squares the curve exactly.
- `smooth_basis_weather_fourier65_harmonic_lambda0.01-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda100-coefs` (smoothing, rel err 2.01e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-gcv` (smoothing, rel err 9.86e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda100-y2cMap` (smoothing, rel err 2.69e-07): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-fitted` (smoothing, rel err 1.59e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-coefs` (smoothing, rel err 3.68e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-y2cMap` (smoothing, rel err 4.93e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-fitted` (smoothing, rel err 2.91e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda1e+06-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda1e+06-y2cMap` (smoothing, rel err 4.96e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-gcv` (smoothing, rel err 4.86e-09): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-y2cMap` (smoothing, rel err 1.11e-09): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-coefs` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-df` (smoothing, rel err 1.33e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-gcv` (smoothing, rel err 4.93e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-SSE` (smoothing, rel err 5.09e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-y2cMap` (smoothing, rel err 8.52e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-fitted` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-coefs` (smoothing, rel err 4.38e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-df` (smoothing, rel err 1.61e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-gcv` (smoothing, rel err 2.34e-05): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-SSE` (smoothing, rel err 1.08e-05): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-y2cMap` (smoothing, rel err 6.47e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-fitted` (smoothing, rel err 4.39e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `lambda2gcv_weather_log10lambda_2-gcv` (smoothing, rel err 9.86e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `df2lambda_weather_df10-lambda` (smoothing, rel err 7.70e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so Fabel bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `df2lambda_weather_df20-lambda` (smoothing, rel err 9.83e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so Fabel bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `df2lambda_weather_df30-lambda` (smoothing, rel err 8.56e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so Fabel bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `lambda2gcv_weather_grid_mean-mean_gcv` (smoothing, rel err 1.46e-08): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basisPar_growth_hgtm-coefs` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-df` (smoothing, rel err 1.33e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-gcv` (smoothing, rel err 4.93e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-SSE` (smoothing, rel err 5.09e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-y2cMap` (smoothing, rel err 8.52e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-fitted` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `data2fd_weather_fourier-coefs` (smoothing, rel err 2.01e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form Fabel is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `data2fd_growth_hgtf-coefs` (smoothing, rel err 5.80e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) Fabel is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the Fabel-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_monotone_growth_hgtf_3girls-Wfdobj_coefs` (smoothing, rel err 8.92e-03): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_growth_hgtf_3girls-beta_slope` (smoothing, rel err 7.16e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_growth_hgtf_3girls-deriv1` (smoothing, rel err 6.77e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-Wfdobj_coefs` (smoothing, rel err 2.33e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-beta_slope` (smoothing, rel err 2.00e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-deriv1` (smoothing, rel err 1.96e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. Fabel and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_pos_synthetic_2curves-Wfdobj_coefs` (smoothing, rel err 6.49e-05): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for Fabel, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum Fabel reports.
- `smooth_pos_synthetic_2curves-fitted` (smoothing, rel err 1.75e-05): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for Fabel, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum Fabel reports.
- `smooth_pos_weather_precip_pr_rupert-Wfdobj_coefs` (smoothing, rel err 5.03e-08): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for Fabel, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum Fabel reports.

## Every check

| Module | Case | Symbol | rtol | Max abs err | Max rel err | Status |
|---|---|---|---:|---:|---:|---|
| basis | `bspline_eval_k1_n4_dom0_1_d0` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k1_n7_dom0_1_d0` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k1_n10_dom0_1_d0` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k1_n15_dom0_1_d0` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k1_n23_dom0_1_d0` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k2_n4_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k2_n4_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k2_n4_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.39e-16 | 6.25e-16 | pass |
| basis | `bspline_eval_k2_n7_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k2_n7_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k2_n7_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.25e-16 | 1.12e-15 | pass |
| basis | `bspline_eval_k2_n10_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k2_n10_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k2_n10_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.53e-16 | 2.06e-15 | pass |
| basis | `bspline_eval_k2_n15_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k2_n15_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k2_n15_dom0_1_L0` | `BSpline` | 1.00e-08 | 6.25e-17 | 1.31e-15 | pass |
| basis | `bspline_eval_k2_n23_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k2_n23_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k2_n23_dom0_1_L0` | `BSpline` | 1.00e-08 | 7.29e-17 | 2.40e-15 | pass |
| basis | `bspline_eval_k3_n4_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k3_n4_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k3_n4_dom0_1_d2` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k3_n4_dom0_1_L0` | `BSpline` | 1.00e-08 | 4.86e-17 | 2.91e-16 | pass |
| basis | `bspline_penalty_k3_n4_dom0_1_L1` | `BSpline` | 1.00e-08 | 2.22e-16 | 8.33e-17 | pass |
| basis | `bspline_eval_k3_n7_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k3_n7_dom0_1_d1` | `BSpline` | 1.00e-08 | 1.78e-15 | 1.78e-16 | pass |
| basis | `bspline_eval_k3_n7_dom0_1_d2` | `BSpline` | 1.00e-08 | 7.11e-15 | 9.47e-17 | pass |
| basis | `bspline_penalty_k3_n7_dom0_1_L0` | `BSpline` | 1.00e-08 | 3.47e-17 | 3.15e-16 | pass |
| basis | `bspline_penalty_k3_n7_dom0_1_L1` | `BSpline` | 1.00e-08 | 7.11e-15 | 1.07e-15 | pass |
| basis | `bspline_eval_k3_n10_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k3_n10_dom0_1_d1` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_eval_k3_n10_dom0_1_d2` | `BSpline` | 1.00e-08 | 0 | 0 | pass |
| basis | `bspline_penalty_k3_n10_dom0_1_L0` | `BSpline` | 1.00e-08 | 4.51e-17 | 6.56e-16 | pass |
| basis | `bspline_penalty_k3_n10_dom0_1_L1` | `BSpline` | 1.00e-08 | 1.60e-14 | 1.50e-15 | pass |
| basis | `bspline_eval_k3_n15_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k3_n15_dom0_1_d1` | `BSpline` | 1.00e-08 | 1.78e-15 | 6.83e-17 | pass |
| basis | `bspline_eval_k3_n15_dom0_1_d2` | `BSpline` | 1.00e-08 | 5.68e-14 | 1.12e-16 | pass |
| basis | `bspline_penalty_k3_n15_dom0_1_L0` | `BSpline` | 1.00e-08 | 4.16e-17 | 9.84e-16 | pass |
| basis | `bspline_penalty_k3_n15_dom0_1_L1` | `BSpline` | 1.00e-08 | 3.20e-14 | 1.84e-15 | pass |
| basis | `bspline_eval_k3_n23_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k3_n23_dom0_1_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 8.46e-17 | pass |
| basis | `bspline_eval_k3_n23_dom0_1_d2` | `BSpline` | 1.00e-08 | 2.27e-13 | 1.72e-16 | pass |
| basis | `bspline_penalty_k3_n23_dom0_1_L0` | `BSpline` | 1.00e-08 | 2.95e-17 | 1.13e-15 | pass |
| basis | `bspline_penalty_k3_n23_dom0_1_L1` | `BSpline` | 1.00e-08 | 7.28e-14 | 2.60e-15 | pass |
| basis | `bspline_eval_k4_n4_dom0_1_d0` | `BSpline` | 1.00e-08 | 5.55e-17 | 5.55e-17 | pass |
| basis | `bspline_eval_k4_n4_dom0_1_d1` | `BSpline` | 1.00e-08 | 4.44e-16 | 1.48e-16 | pass |
| basis | `bspline_eval_k4_n4_dom0_1_d2` | `BSpline` | 1.00e-08 | 1.78e-15 | 1.48e-16 | pass |
| basis | `bspline_penalty_k4_n4_dom0_1_L0` | `BSpline` | 1.00e-08 | 8.57e-01 | 8.57e-01 | xfail |
| basis | `bspline_penalty_k4_n4_dom0_1_L1` | `BSpline` | 1.00e-08 | 2.40e+00 | 1.33e+00 | xfail |
| basis | `bspline_penalty_k4_n4_dom0_1_L2` | `BSpline` | 1.00e-08 | 3.60e+01 | 3.00e+00 | xfail |
| basis | `bspline_eval_k4_n7_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k4_n7_dom0_1_d1` | `BSpline` | 1.00e-08 | 1.78e-15 | 1.48e-16 | pass |
| basis | `bspline_eval_k4_n7_dom0_1_d2` | `BSpline` | 1.00e-08 | 1.42e-14 | 9.87e-17 | pass |
| basis | `bspline_penalty_k4_n7_dom0_1_L0` | `BSpline` | 1.00e-08 | 8.67e-17 | 7.24e-16 | pass |
| basis | `bspline_penalty_k4_n7_dom0_1_L1` | `BSpline` | 1.00e-08 | 7.11e-15 | 9.87e-16 | pass |
| basis | `bspline_penalty_k4_n7_dom0_1_L2` | `BSpline` | 1.00e-08 | 4.55e-13 | 2.96e-16 | pass |
| basis | `bspline_eval_k4_n10_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k4_n10_dom0_1_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 1.69e-16 | pass |
| basis | `bspline_eval_k4_n10_dom0_1_d2` | `BSpline` | 1.00e-08 | 5.68e-14 | 1.29e-16 | pass |
| basis | `bspline_penalty_k4_n10_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.51e-16 | 2.20e-15 | pass |
| basis | `bspline_penalty_k4_n10_dom0_1_L1` | `BSpline` | 1.00e-08 | 1.60e-14 | 1.27e-15 | pass |
| basis | `bspline_penalty_k4_n10_dom0_1_L2` | `BSpline` | 1.00e-08 | 5.46e-12 | 6.63e-16 | pass |
| basis | `bspline_eval_k4_n15_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k4_n15_dom0_1_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 9.87e-17 | pass |
| basis | `bspline_eval_k4_n15_dom0_1_d2` | `BSpline` | 1.00e-08 | 1.14e-13 | 8.77e-17 | pass |
| basis | `bspline_penalty_k4_n15_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.11e-16 | 2.78e-15 | pass |
| basis | `bspline_penalty_k4_n15_dom0_1_L1` | `BSpline` | 1.00e-08 | 1.39e-13 | 6.41e-15 | pass |
| basis | `bspline_penalty_k4_n15_dom0_1_L2` | `BSpline` | 1.00e-08 | 1.38e-10 | 3.33e-15 | pass |
| basis | `bspline_eval_k4_n23_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k4_n23_dom0_1_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 5.92e-17 | pass |
| basis | `bspline_eval_k4_n23_dom0_1_d2` | `BSpline` | 1.00e-08 | 4.55e-13 | 1.26e-16 | pass |
| basis | `bspline_penalty_k4_n23_dom0_1_L0` | `BSpline` | 1.00e-08 | 6.42e-17 | 2.68e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom0_1_L1` | `BSpline` | 1.00e-08 | 2.13e-13 | 5.92e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom0_1_L2` | `BSpline` | 1.00e-08 | 6.40e-10 | 3.33e-15 | pass |
| basis | `bspline_eval_k6_n7_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_eval_k6_n7_dom0_1_d1` | `BSpline` | 1.00e-08 | 1.78e-15 | 1.78e-16 | pass |
| basis | `bspline_eval_k6_n7_dom0_1_d2` | `BSpline` | 1.00e-08 | 1.42e-14 | 1.18e-16 | pass |
| basis | `bspline_penalty_k6_n7_dom0_1_L0` | `BSpline` | 1.00e-08 | 2.62e-15 | 3.94e-14 | pass |
| basis | `bspline_penalty_k6_n7_dom0_1_L1` | `BSpline` | 1.00e-08 | 5.51e-14 | 9.91e-15 | pass |
| basis | `bspline_penalty_k6_n7_dom0_1_L2` | `BSpline` | 1.00e-08 | 7.96e-13 | 9.22e-16 | pass |
| basis | `bspline_eval_k6_n10_dom0_1_d0` | `BSpline` | 1.00e-08 | 1.94e-16 | 1.94e-16 | pass |
| basis | `bspline_eval_k6_n10_dom0_1_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 1.42e-16 | pass |
| basis | `bspline_eval_k6_n10_dom0_1_d2` | `BSpline` | 1.00e-08 | 1.14e-13 | 1.52e-16 | pass |
| basis | `bspline_penalty_k6_n10_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.71e-15 | 2.83e-14 | pass |
| basis | `bspline_penalty_k6_n10_dom0_1_L1` | `BSpline` | 1.00e-08 | 8.17e-13 | 5.88e-14 | pass |
| basis | `bspline_penalty_k6_n10_dom0_1_L2` | `BSpline` | 1.00e-08 | 8.73e-11 | 6.48e-15 | pass |
| basis | `bspline_eval_k6_n15_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k6_n15_dom0_1_d1` | `BSpline` | 1.00e-08 | 9.33e-15 | 1.87e-16 | pass |
| basis | `bspline_eval_k6_n15_dom0_1_d2` | `BSpline` | 1.00e-08 | 4.55e-13 | 1.52e-16 | pass |
| basis | `bspline_penalty_k6_n15_dom0_1_L0` | `BSpline` | 1.00e-08 | 8.55e-16 | 2.17e-14 | pass |
| basis | `bspline_penalty_k6_n15_dom0_1_L1` | `BSpline` | 1.00e-08 | 1.63e-12 | 5.88e-14 | pass |
| basis | `bspline_penalty_k6_n15_dom0_1_L2` | `BSpline` | 1.00e-08 | 6.98e-10 | 6.48e-15 | pass |
| basis | `bspline_eval_k6_n23_dom0_1_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k6_n23_dom0_1_d1` | `BSpline` | 1.00e-08 | 7.11e-15 | 7.89e-17 | pass |
| basis | `bspline_eval_k6_n23_dom0_1_d2` | `BSpline` | 1.00e-08 | 9.09e-13 | 9.36e-17 | pass |
| basis | `bspline_penalty_k6_n23_dom0_1_L0` | `BSpline` | 1.00e-08 | 1.18e-15 | 5.39e-14 | pass |
| basis | `bspline_penalty_k6_n23_dom0_1_L1` | `BSpline` | 1.00e-08 | 2.96e-12 | 5.93e-14 | pass |
| basis | `bspline_penalty_k6_n23_dom0_1_L2` | `BSpline` | 1.00e-08 | 6.40e-09 | 1.02e-14 | pass |
| basis | `bspline_eval_k4_n23_dom0_365_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k4_n23_dom0_365_d1` | `BSpline` | 1.00e-08 | 6.94e-18 | 4.22e-17 | pass |
| basis | `bspline_eval_k4_n23_dom0_365_d2` | `BSpline` | 1.00e-08 | 3.47e-18 | 1.28e-16 | pass |
| basis | `bspline_penalty_k4_n23_dom0_365_L0` | `BSpline` | 1.00e-08 | 2.53e-14 | 2.89e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom0_365_L1` | `BSpline` | 1.00e-08 | 7.91e-16 | 8.02e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom0_365_L2` | `BSpline` | 1.00e-08 | 1.47e-17 | 3.73e-15 | pass |
| basis | `bspline_eval_k4_n23_dom-1_2_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_eval_k4_n23_dom-1_2_d1` | `BSpline` | 1.00e-08 | 1.78e-15 | 8.88e-17 | pass |
| basis | `bspline_eval_k4_n23_dom-1_2_d2` | `BSpline` | 1.00e-08 | 5.68e-14 | 1.42e-16 | pass |
| basis | `bspline_penalty_k4_n23_dom-1_2_L0` | `BSpline` | 1.00e-08 | 2.05e-16 | 2.85e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom-1_2_L1` | `BSpline` | 1.00e-08 | 4.62e-14 | 3.85e-15 | pass |
| basis | `bspline_penalty_k4_n23_dom-1_2_L2` | `BSpline` | 1.00e-08 | 1.36e-11 | 1.92e-15 | pass |
| basis | `bspline_breaks1_eval_d0` | `BSpline` | 1.00e-08 | 4.44e-16 | 4.44e-16 | pass |
| basis | `bspline_breaks1_eval_d1` | `BSpline` | 1.00e-08 | 3.55e-15 | 1.18e-16 | pass |
| basis | `bspline_breaks1_eval_d2` | `BSpline` | 1.00e-08 | 1.14e-13 | 1.14e-16 | pass |
| basis | `bspline_breaks1_penalty_L0` | `BSpline` | 1.00e-08 | 3.05e-16 | 2.66e-15 | pass |
| basis | `bspline_breaks1_penalty_L1` | `BSpline` | 1.00e-08 | 3.20e-14 | 1.78e-15 | pass |
| basis | `bspline_breaks1_penalty_L2` | `BSpline` | 1.00e-08 | 5.46e-12 | 1.95e-16 | pass |
| basis | `bspline_breaks2_eval_d0` | `BSpline` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `bspline_breaks2_eval_d1` | `BSpline` | 1.00e-08 | 6.94e-18 | 1.16e-16 | pass |
| basis | `bspline_breaks2_eval_d2` | `BSpline` | 1.00e-08 | 4.34e-19 | 1.28e-16 | pass |
| basis | `bspline_breaks2_penalty_L0` | `BSpline` | 1.00e-08 | 2.66e-14 | 5.88e-16 | pass |
| basis | `bspline_breaks2_penalty_L1` | `BSpline` | 1.00e-08 | 9.02e-17 | 2.51e-15 | pass |
| basis | `bspline_breaks2_penalty_L2` | `BSpline` | 1.00e-08 | 1.36e-19 | 7.70e-16 | pass |
| basis | `bspline_breaks3_eval_d0` | `BSpline` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `bspline_breaks3_eval_d1` | `BSpline` | 1.00e-08 | 4.44e-16 | 7.40e-17 | pass |
| basis | `bspline_breaks3_eval_d2` | `BSpline` | 1.00e-08 | 1.78e-15 | 4.93e-17 | pass |
| basis | `bspline_breaks3_penalty_L0` | `BSpline` | 1.00e-08 | 2.22e-16 | 7.12e-16 | pass |
| basis | `bspline_breaks3_penalty_L1` | `BSpline` | 1.00e-08 | 4.61e-15 | 1.28e-15 | pass |
| basis | `bspline_breaks3_penalty_L2` | `BSpline` | 1.00e-08 | 5.68e-14 | 2.96e-16 | pass |
| basis | `fourier_eval_n3_d0` | `Fourier` | 1.00e-08 | 2.22e-16 | 1.57e-16 | pass |
| basis | `fourier_eval_n3_d1` | `Fourier` | 1.00e-08 | 1.78e-15 | 2.00e-16 | pass |
| basis | `fourier_eval_n3_d2` | `Fourier` | 1.00e-08 | 1.42e-14 | 2.55e-16 | pass |
| basis | `fourier_eval_n3_d3` | `Fourier` | 1.00e-08 | 5.68e-14 | 1.62e-16 | pass |
| basis | `fourier_penalty_n3_L0` | `Fourier` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `fourier_penalty_n3_L1` | `Fourier` | 1.00e-08 | 0 | 0 | pass |
| basis | `fourier_penalty_n3_L2` | `Fourier` | 1.00e-08 | 2.27e-13 | 1.46e-16 | pass |
| basis | `fourier_eval_n5_d0` | `Fourier` | 1.00e-08 | 2.22e-16 | 1.57e-16 | pass |
| basis | `fourier_eval_n5_d1` | `Fourier` | 1.00e-08 | 3.55e-15 | 2.00e-16 | pass |
| basis | `fourier_eval_n5_d2` | `Fourier` | 1.00e-08 | 5.68e-14 | 2.55e-16 | pass |
| basis | `fourier_eval_n5_d3` | `Fourier` | 1.00e-08 | 4.55e-13 | 1.62e-16 | pass |
| basis | `fourier_penalty_n5_L0` | `Fourier` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `fourier_penalty_n5_L1` | `Fourier` | 1.00e-08 | 0 | 0 | pass |
| basis | `fourier_penalty_n5_L2` | `Fourier` | 1.00e-08 | 3.64e-12 | 1.46e-16 | pass |
| basis | `fourier_eval_n9_d0` | `Fourier` | 1.00e-08 | 4.55e-15 | 3.22e-15 | pass |
| basis | `fourier_eval_n9_d1` | `Fourier` | 1.00e-08 | 8.53e-14 | 2.40e-15 | pass |
| basis | `fourier_eval_n9_d2` | `Fourier` | 1.00e-08 | 1.62e-12 | 1.81e-15 | pass |
| basis | `fourier_eval_n9_d3` | `Fourier` | 1.00e-08 | 3.09e-11 | 1.38e-15 | pass |
| basis | `fourier_penalty_n9_L0` | `Fourier` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `fourier_penalty_n9_L1` | `Fourier` | 1.00e-08 | 0 | 0 | pass |
| basis | `fourier_penalty_n9_L2` | `Fourier` | 1.00e-08 | 5.82e-11 | 1.46e-16 | pass |
| basis | `fourier_eval_n65_d0` | `Fourier` | 1.00e-08 | 4.02e-14 | 2.84e-14 | pass |
| basis | `fourier_eval_n65_d1` | `Fourier` | 1.00e-08 | 7.74e-12 | 2.72e-14 | pass |
| basis | `fourier_eval_n65_d2` | `Fourier` | 1.00e-08 | 1.51e-09 | 2.63e-14 | pass |
| basis | `fourier_eval_n65_d3` | `Fourier` | 1.00e-08 | 2.94e-07 | 2.55e-14 | pass |
| basis | `fourier_penalty_n65_L0` | `Fourier` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `fourier_penalty_n65_L1` | `Fourier` | 1.00e-08 | 1.46e-11 | 3.60e-16 | pass |
| basis | `fourier_penalty_n65_L2` | `Fourier` | 1.00e-08 | 4.77e-07 | 2.92e-16 | pass |
| basis | `fourier_periodmismatch_eval_d0` | `Fourier` | 1.00e-08 | 1.72e-15 | 1.72e-15 | pass |
| basis | `fourier_periodmismatch_eval_d1` | `Fourier` | 1.00e-08 | 1.64e-14 | 1.31e-15 | pass |
| basis | `fourier_periodmismatch_eval_d2` | `Fourier` | 1.00e-08 | 1.53e-13 | 9.67e-16 | pass |
| basis | `fourier_periodmismatch_eval_d3` | `Fourier` | 1.00e-08 | 1.45e-12 | 7.30e-16 | pass |
| basis | `fourier_periodmismatch_penalty_L0` | `Fourier` | 1.00e-08 | 6.92e-07 | 1.38e-06 | xfail |
| basis | `fourier_periodmismatch_penalty_L1` | `Fourier` | 1.00e-08 | 1.09e-04 | 1.38e-06 | xfail |
| basis | `fourier_periodmismatch_penalty_L2` | `Fourier` | 1.00e-08 | 1.72e-02 | 1.38e-06 | xfail |
| basis | `fourier_even_nbasis_request_eval_d0` | `Fourier` | 1.00e-08 | 2.22e-16 | 1.57e-16 | pass |
| basis | `fourier_eval_dom365_n65_d0` | `Fourier` | 1.00e-08 | 2.10e-15 | 2.84e-14 | pass |
| basis | `fourier_eval_dom365_n65_d1` | `Fourier` | 1.00e-08 | 1.11e-15 | 2.72e-14 | pass |
| basis | `fourier_eval_dom365_n65_d2` | `Fourier` | 1.00e-08 | 5.92e-16 | 2.64e-14 | pass |
| basis | `fourier_eval_dom365_n65_d3` | `Fourier` | 1.00e-08 | 3.16e-16 | 2.55e-14 | pass |
| basis | `fourier_penalty_dom365_n65_L0` | `Fourier` | 1.00e-08 | 2.22e-16 | 2.22e-16 | pass |
| basis | `fourier_penalty_dom365_n65_L1` | `Fourier` | 1.00e-08 | 5.55e-17 | 1.83e-16 | pass |
| basis | `fourier_penalty_dom365_n65_L2` | `Fourier` | 1.00e-08 | 2.78e-17 | 3.01e-16 | pass |
| basis | `monomial_eval_n1_d0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n1_d1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n1_d2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n1_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n1_L1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n1_L2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n2_d0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n2_d1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n2_d2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n2_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n2_L1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n2_L2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n3_d0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n3_d1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n3_d2` | `Monomial` | 1.00e-08 | 2.00e+00 | inf | xfail |
| basis | `monomial_penalty_n3_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n3_L1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n3_L2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n4_d0` | `Monomial` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `monomial_eval_n4_d1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n4_d2` | `Monomial` | 1.00e-08 | 3.00e+00 | 1.00e+00 | xfail |
| basis | `monomial_penalty_n4_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n4_L1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n4_L2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n5_d0` | `Monomial` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `monomial_eval_n5_d1` | `Monomial` | 1.00e-08 | 4.44e-16 | 1.11e-16 | pass |
| basis | `monomial_eval_n5_d2` | `Monomial` | 1.00e-08 | 4.00e+00 | 5.00e-01 | xfail |
| basis | `monomial_penalty_n5_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n5_L1` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n5_L2` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_eval_n6_d0` | `Monomial` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `monomial_eval_n6_d1` | `Monomial` | 1.00e-08 | 8.88e-16 | 1.78e-16 | pass |
| basis | `monomial_eval_n6_d2` | `Monomial` | 1.00e-08 | 5.00e+00 | 3.33e-01 | xfail |
| basis | `monomial_penalty_n6_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_penalty_n6_L1` | `Monomial` | 1.00e-08 | 2.22e-16 | 7.99e-17 | pass |
| basis | `monomial_penalty_n6_L2` | `Monomial` | 1.00e-08 | 7.11e-15 | 1.24e-16 | pass |
| basis | `monomial_customexp_eval_d0` | `Monomial` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `monomial_customexp_eval_d1` | `Monomial` | 1.00e-08 | 8.88e-16 | 1.78e-16 | pass |
| basis | `monomial_customexp_eval_d2` | `Monomial` | 1.00e-08 | 5.00e+00 | 3.33e-01 | xfail |
| basis | `monomial_customexp_penalty_L0` | `Monomial` | 1.00e-08 | 0 | 0 | pass |
| basis | `monomial_customexp_penalty_L1` | `Monomial` | 1.00e-08 | 2.22e-16 | 7.99e-17 | pass |
| basis | `monomial_customexp_penalty_L2` | `Monomial` | 1.00e-08 | 7.11e-15 | 1.24e-16 | pass |
| basis | `exponential_cfg1_eval_d0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg1_eval_d1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg1_eval_d2` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg1_penalty_L0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg1_penalty_L1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg1_penalty_L2` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_eval_d0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_eval_d1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_eval_d2` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_penalty_L0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_penalty_L1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg2_penalty_L2` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_eval_d0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_eval_d1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_eval_d2` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_penalty_L0` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_penalty_L1` | `Exponential` | 1.00e-08 | 0 | 0 | pass |
| basis | `exponential_cfg3_penalty_L2` | `Exponential` | 1.00e-08 | 9.09e-13 | 1.67e-16 | pass |
| basis | `power_cfg1_eval_d0` | `Power` | 1.00e-08 | 0 | 0 | pass |
| basis | `power_cfg1_eval_d1` | `Power` | 1.00e-08 | 4.00e+00 | 5.00e-01 | xfail |
| basis | `power_cfg1_penalty_L0` | `Power` | 1.00e-08 | 0 | 0 | pass |
| basis | `power_cfg1_penalty_L1` | `Power` | 1.00e-08 | 0 | 0 | pass |
| basis | `power_cfg2_eval_d0` | `Power` | 1.00e-08 | 0 | 0 | pass |
| basis | `power_cfg2_eval_d1` | `Power` | 1.00e-08 | 4.00e+00 | inf | xfail |
| basis | `power_cfg2_penalty_L0` | `Power` | 1.00e-08 | 0 | 0 | xfail |
| basis | `power_cfg2_penalty_L1` | `Power` | 1.00e-08 | 0 | 0 | pass |
| basis | `constant_eval_dom0_1_d0` | `Constant` | 1.00e-08 | 0 | 0 | pass |
| basis | `constant_penalty_dom0_1_L0` | `Constant` | 1.00e-08 | 0 | 0 | pass |
| basis | `constant_eval_dom0_365_d0` | `Constant` | 1.00e-08 | 0 | 0 | pass |
| basis | `constant_penalty_dom0_365_L0` | `Constant` | 1.00e-08 | 0 | 0 | pass |
| basis | `polygonal_n5_eval_d0` | `Polygonal` | 1.00e-08 | 0 | 0 | pass |
| basis | `polygonal_n5_eval_d1` | `Polygonal` | 1.00e-08 | 5.00e+00 | 5.00e+00 | xfail |
| basis | `polygonal_n5_penalty_L0` | `Polygonal` | 1.00e-08 | 8.33e-17 | 5.00e-16 | pass |
| basis | `polygonal_n5_penalty_L1` | `Polygonal` | 1.00e-08 | 0 | 0 | pass |
| basis | `polygonal_n11_eval_d0` | `Polygonal` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `polygonal_n11_eval_d1` | `Polygonal` | 1.00e-08 | 1.10e+01 | 1.10e+01 | xfail |
| basis | `polygonal_n11_penalty_L0` | `Polygonal` | 1.00e-08 | 1.53e-16 | 2.29e-15 | pass |
| basis | `polygonal_n11_penalty_L1` | `Polygonal` | 1.00e-08 | 3.55e-15 | 1.78e-16 | pass |
| basis | `polygonal_n20_eval_d0` | `Polygonal` | 1.00e-08 | 1.11e-16 | 1.11e-16 | pass |
| basis | `polygonal_n20_eval_d1` | `Polygonal` | 1.00e-08 | 2.00e+01 | 2.00e+01 | xfail |
| basis | `polygonal_n20_penalty_L0` | `Polygonal` | 1.00e-08 | 9.71e-17 | 2.77e-15 | pass |
| basis | `polygonal_n20_penalty_L1` | `Polygonal` | 1.00e-08 | 7.11e-15 | 1.87e-16 | pass |
| core | `evalfd_bspline_n3curves_d0` | `FData` | 1.00e-06 | 3.33e-16 | 2.23e-16 | pass |
| core | `evalfd_fourier_n3curves_d0` | `FData` | 1.00e-06 | 9.77e-15 | 1.55e-15 | pass |
| core | `evalfd_bspline_n3curves_d1` | `FData` | 1.00e-06 | 7.11e-15 | 1.62e-16 | pass |
| core | `evalfd_fourier_n3curves_d1` | `FData` | 1.00e-06 | 9.95e-14 | 1.00e-15 | pass |
| core | `evalfd_bspline_n3curves_d2` | `FData` | 1.00e-06 | 2.27e-13 | 2.52e-16 | pass |
| core | `evalfd_fourier_n3curves_d2` | `FData` | 1.00e-06 | 3.41e-12 | 1.77e-15 | pass |
| core | `inprod_bspline_L0_0` | `inprod` | 1.00e-06 | 9.99e-16 | 1.63e-15 | pass |
| core | `inprod_bspline_L1_1` | `inprod` | 1.00e-06 | 9.95e-14 | 1.61e-15 | pass |
| core | `inprod_bspline_L2_0` | `inprod` | 1.00e-06 | 2.49e-14 | 5.00e-16 | pass |
| core | `inprod_bspline_L0_2` | `inprod` | 1.00e-06 | 5.68e-14 | 1.14e-15 | pass |
| core | `inprod_fourier_L0_0` | `inprod` | 1.00e-06 | 4.21e-04 | 4.26e-05 | xfail |
| core | `inprod_fourier_L1_1` | `inprod` | 1.00e-06 | 2.65e-01 | 1.17e-04 | xfail |
| core | `inprod_fourier_L2_2` | `inprod` | 1.00e-06 | 1.67e+02 | 1.99e-04 | xfail |
| core | `inprod_basis_bspline_x_fourier` | `inprod` | 1.00e-06 | 9.85e-06 | 5.73e-05 | xfail |
| core | `mean_fd_bspline_n10curves` | `FData` | 1.00e-06 | 1.11e-16 | 2.06e-16 | pass |
| core | `sd_fd_bspline_n10curves` | `FData` | 1.00e-06 | 8.88e-16 | 9.08e-16 | pass |
| core | `center_fd_bspline_n10curves` | `FData` | 1.00e-06 | 2.22e-16 | 1.16e-16 | pass |
| core | `var_fd_bspline_n10curves` | `FData` | 1.00e-06 | 2.22e-16 | 2.06e-16 | pass |
| core | `deriv_fd_bspline_order6_L2` | `FData` | 1.00e-06 | 5.99e+01 | 1.27e-02 | xfail |
| core | `deriv_fd_fourier_L1` | `FData` | 1.00e-06 | 6.31e-14 | 6.27e-16 | pass |
| core | `fd_add_bspline` | `FData` | 1.00e-06 | 0 | 0 | pass |
| core | `fd_mul_bspline` | `FData` | 1.00e-06 | 2.54e-01 | 1.31e-01 | xfail |
| core | `fd_scalar_mul` | `FData` | 1.00e-06 | 0 | 0 | pass |
| core | `fd_power2` | `FData` | 1.00e-06 | 1.16e-04 | 5.14e-05 | xfail |
| core | `evalfd_lfd_int2lfd2_bspline_order6` | `LDO` | 1.00e-06 | 4.55e-13 | 9.77e-17 | pass |
| core | `evalfd_lfd_vec2lfd_bspline` | `LDO` | 1.00e-06 | 1.14e-13 | 1.22e-16 | pass |
| core | `evalfd_lfd_harmonic_accelerator_fourier` | `LDO` | 1.00e-06 | 2.33e-18 | 5.31e-15 | pass |
| core | `bifd_construct_4d_coef` | `BiFData` | 1.00e-06 | 0 | 0 | pass |
| datasets | `test_canadian_weather` | `datasets.load_*` | 1.00e-12 | 3.21e-09 | 5.94e-14 | pass |
| datasets | `test_growth` | `datasets.load_*` | 1.00e-12 | 2.62e-10 | 1.56e-15 | pass |
| datasets | `test_gait` | `datasets.load_*` | 1.00e-12 | 0 | 0 | pass |
| datasets | `test_handwrit` | `datasets.load_*` | 1.00e-12 | 1.17e-14 | 2.32e-14 | pass |
| datasets | `test_handwrit_time` | `datasets.load_*` | 1.00e-12 | 2.33e-10 | 1.45e-16 | pass |
| datasets | `test_pinch` | `datasets.load_*` | 1.00e-12 | 1.18e-11 | 2.29e-15 | pass |
| datasets | `test_melanoma` | `datasets.load_*` | 1.00e-12 | 2.91e-11 | 3.98e-16 | pass |
| datasets | `test_refinery` | `datasets.load_*` | 1.00e-12 | 5.68e-14 | 3.38e-16 | pass |
| datasets | `test_seabird` | `datasets.load_*` | 1.00e-12 | 9.09e-12 | 7.37e-16 | pass |
| datasets | `test_regina_precip` | `datasets.load_*` | 1.00e-12 | 4.55e-13 | 2.15e-16 | pass |
| datasets | `test_catalog_scope` | `datasets.load_*` | — | — | — | pass |
| io | `test_fd_coefs_and_eval_match_r` | `read_rds` | 1.00e-08 | 2.22e-16 | 9.98e-17 | pass |
| io | `test_bifd_coefs_and_eval_match_r` | `read_rds` | 1.00e-08 | 5.55e-17 | 2.43e-17 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-coefs` | `smooth` | 1.00e-08 | 9.09e-10 | 2.88e-12 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-df` | `smooth` | 1.00e-08 | 1.42e-14 | 2.19e-16 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-gcv` | `smooth` | 1.00e-08 | 3.22e-15 | 5.25e-15 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-SSE` | `smooth` | 1.00e-08 | 1.77e-11 | 5.22e-15 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-penmat` | `smooth` | 1.00e-08 | 3.86e-08 | 1.38e-06 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-y2cMap` | `smooth` | 1.00e-08 | 2.85e-11 | 3.86e-10 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda0.01-fitted` | `smooth` | 1.00e-08 | 7.60e-11 | 2.27e-12 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-coefs` | `smooth` | 1.00e-08 | 6.33e-07 | 2.01e-09 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-df` | `smooth` | 1.00e-08 | 5.26e-13 | 9.96e-15 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-gcv` | `smooth` | 1.00e-08 | 5.97e-09 | 9.86e-09 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-SSE` | `smooth` | 1.00e-08 | 3.96e-06 | 1.08e-09 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-penmat` | `smooth` | 1.00e-08 | 3.86e-08 | 1.38e-06 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-y2cMap` | `smooth` | 1.00e-08 | 1.99e-08 | 2.69e-07 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda100-fitted` | `smooth` | 1.00e-08 | 5.29e-08 | 1.59e-09 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-coefs` | `smooth` | 1.00e-08 | 1.16e-08 | 3.68e-11 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-df` | `smooth` | 1.00e-08 | 3.20e-14 | 1.22e-15 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-gcv` | `smooth` | 1.00e-08 | 1.26e-10 | 1.26e-10 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-SSE` | `smooth` | 1.00e-08 | 9.83e-08 | 1.70e-11 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-penmat` | `smooth` | 1.00e-08 | 3.86e-08 | 1.38e-06 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-y2cMap` | `smooth` | 1.00e-08 | 3.65e-10 | 4.93e-09 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda10000-fitted` | `smooth` | 1.00e-08 | 9.70e-10 | 2.91e-11 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-coefs` | `smooth` | 1.00e-08 | 1.17e-10 | 3.70e-13 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-df` | `smooth` | 1.00e-08 | 3.55e-15 | 2.89e-16 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-gcv` | `smooth` | 1.00e-08 | 1.17e-12 | 9.43e-13 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-SSE` | `smooth` | 1.00e-08 | 9.98e-10 | 1.32e-13 | pass |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-penmat` | `smooth` | 1.00e-08 | 3.86e-08 | 1.38e-06 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-y2cMap` | `smooth` | 1.00e-08 | 3.67e-12 | 4.96e-11 | xfail |
| smoothing | `smooth_basis_weather_fourier65_harmonic_lambda1e+06-fitted` | `smooth` | 1.00e-08 | 9.78e-12 | 2.94e-13 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-coefs` | `smooth` | 1.00e-08 | 1.42e-07 | 7.27e-10 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-df` | `smooth` | 1.00e-08 | 1.37e-09 | 1.08e-10 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-gcv` | `smooth` | 1.00e-08 | 2.95e-09 | 4.86e-09 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-SSE` | `smooth` | 1.00e-08 | 3.12e-07 | 3.54e-09 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-penmat` | `smooth` | 1.00e-08 | 5.96e-08 | 2.23e-16 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-y2cMap` | `smooth` | 1.00e-08 | 1.03e-09 | 1.11e-09 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-fitted` | `smooth` | 1.00e-08 | 1.42e-07 | 7.27e-10 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-coefs` | `smooth` | 1.00e-08 | 1.09e-05 | 5.63e-08 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-df` | `smooth` | 1.00e-08 | 1.07e-07 | 1.33e-08 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-gcv` | `smooth` | 1.00e-08 | 6.25e-07 | 4.93e-07 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-SSE` | `smooth` | 1.00e-08 | 2.03e-04 | 5.09e-07 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-penmat` | `smooth` | 1.00e-08 | 5.96e-08 | 2.23e-16 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-y2cMap` | `smooth` | 1.00e-08 | 6.47e-08 | 8.52e-08 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-fitted` | `smooth` | 1.00e-08 | 1.09e-05 | 5.63e-08 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-coefs` | `smooth` | 1.00e-08 | 8.58e-04 | 4.38e-06 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-df` | `smooth` | 1.00e-08 | 8.74e-06 | 1.61e-06 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-gcv` | `smooth` | 1.00e-08 | 8.80e-05 | 2.34e-05 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-SSE` | `smooth` | 1.00e-08 | 1.59e-02 | 1.08e-05 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-penmat` | `smooth` | 1.00e-08 | 5.96e-08 | 2.23e-16 | pass |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-y2cMap` | `smooth` | 1.00e-08 | 3.57e-06 | 6.47e-06 | xfail |
| smoothing | `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-fitted` | `smooth` | 1.00e-08 | 8.58e-04 | 4.39e-06 | xfail |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-coefs` | `smooth` | 1.00e-08 | 4.44e-15 | 1.45e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-gcv` | `smooth` | 1.00e-08 | 9.71e-17 | 8.74e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-SSE` | `smooth` | 1.00e-08 | 2.22e-15 | 6.92e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-y2cMap` | `smooth` | 1.00e-08 | 4.44e-16 | 5.17e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0-fitted` | `smooth` | 1.00e-08 | 2.00e-15 | 1.05e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-coefs` | `smooth` | 1.00e-08 | 3.30e-15 | 1.09e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-gcv` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-SSE` | `smooth` | 1.00e-08 | 4.44e-16 | 1.38e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-y2cMap` | `smooth` | 1.00e-08 | 6.66e-16 | 7.94e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda0.001-fitted` | `smooth` | 1.00e-08 | 1.22e-15 | 6.41e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-coefs` | `smooth` | 1.00e-08 | 2.66e-15 | 2.08e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-df` | `smooth` | 1.00e-08 | 4.44e-15 | 7.27e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-gcv` | `smooth` | 1.00e-08 | 2.22e-16 | 4.41e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-SSE` | `smooth` | 1.00e-08 | 3.55e-15 | 1.83e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-y2cMap` | `smooth` | 1.00e-08 | 5.00e-16 | 1.44e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve1_lambda1-fitted` | `smooth` | 1.00e-08 | 2.66e-15 | 2.19e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-coefs` | `smooth` | 1.00e-08 | 4.44e-15 | 1.31e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-gcv` | `smooth` | 1.00e-08 | 2.78e-17 | 1.82e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-SSE` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-y2cMap` | `smooth` | 1.00e-08 | 4.44e-16 | 5.17e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0-fitted` | `smooth` | 1.00e-08 | 2.22e-15 | 1.07e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-coefs` | `smooth` | 1.00e-08 | 2.89e-15 | 8.53e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-gcv` | `smooth` | 1.00e-08 | 2.78e-17 | 1.83e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-SSE` | `smooth` | 1.00e-08 | 1.78e-15 | 1.92e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-y2cMap` | `smooth` | 1.00e-08 | 6.66e-16 | 7.94e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda0.001-fitted` | `smooth` | 1.00e-08 | 1.67e-15 | 8.04e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-coefs` | `smooth` | 1.00e-08 | 2.22e-15 | 1.52e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-df` | `smooth` | 1.00e-08 | 4.44e-15 | 7.27e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-gcv` | `smooth` | 1.00e-08 | 5.55e-17 | 9.31e-17 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-SSE` | `smooth` | 1.00e-08 | 2.13e-14 | 3.69e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-y2cMap` | `smooth` | 1.00e-08 | 5.00e-16 | 1.44e-15 | pass |
| smoothing | `smooth_basis_synthetic_bspline4_n12_ncurve3_lambda1-fitted` | `smooth` | 1.00e-08 | 2.22e-15 | 1.84e-15 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-coefs` | `smooth` | 1.00e-08 | 1.55e-15 | 8.26e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-df` | `smooth` | 1.00e-08 | 1.78e-15 | 1.65e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-gcv` | `smooth` | 1.00e-08 | 6.94e-18 | 1.25e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-SSE` | `smooth` | 1.00e-08 | 1.11e-16 | 1.61e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-y2cMap` | `smooth` | 1.00e-08 | 7.77e-16 | 8.73e-16 | pass |
| smoothing | `smooth_basis_synthetic_wtvec-fitted` | `smooth` | 1.00e-08 | 1.55e-15 | 8.74e-16 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-coefs` | `smooth` | 1.00e-08 | 1.55e-15 | 1.18e-15 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-gcv` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-SSE` | `smooth` | 1.00e-08 | 3.33e-16 | 5.50e-16 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-penmat` | `smooth` | 1.00e-08 | 1.44e-15 | 2.11e-15 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-y2cMap` | `smooth` | 1.00e-08 | 7.77e-16 | 8.38e-16 | pass |
| smoothing | `smooth_basis_synthetic_lfd_integer0-fitted` | `smooth` | 1.00e-08 | 7.77e-16 | 7.18e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-coefs` | `smooth` | 1.00e-08 | 6.66e-16 | 3.83e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-df` | `smooth` | 1.00e-08 | 8.88e-16 | 1.64e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-gcv` | `smooth` | 1.00e-08 | 2.78e-17 | 2.86e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-SSE` | `smooth` | 1.00e-08 | 1.78e-15 | 5.25e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-penmat` | `smooth` | 1.00e-08 | 2.27e-12 | 7.58e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-y2cMap` | `smooth` | 1.00e-08 | 2.22e-16 | 5.85e-16 | pass |
| smoothing | `smooth_basis_synthetic_y_bare_vector-fitted` | `smooth` | 1.00e-08 | 4.23e-16 | 4.22e-16 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-coefs` | `smooth` | 1.00e-08 | 2.89e-15 | 2.90e-14 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-gcv` | `smooth` | 1.00e-08 | 2.20e-20 | 1.78e-15 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-SSE` | `smooth` | 1.00e-08 | 5.55e-17 | 1.75e-16 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-penmat` | `smooth` | 1.00e-08 | 5.83e-19 | 1.60e-14 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-y2cMap` | `smooth` | 1.00e-08 | 7.91e-16 | 4.02e-15 | pass |
| smoothing | `smooth_basis_handwrit_3d_y-fitted` | `smooth` | 1.00e-08 | 4.61e-16 | 1.17e-14 | pass |
| smoothing | `lambda2gcv_weather_log10lambda_m2-gcv` | `gcv_curve` | 1.00e-08 | 3.30e-15 | 5.38e-15 | pass |
| smoothing | `lambda2gcv_weather_log10lambda_0-gcv` | `gcv_curve` | 1.00e-08 | 3.23e-11 | 5.28e-11 | pass |
| smoothing | `lambda2gcv_weather_log10lambda_2-gcv` | `gcv_curve` | 1.00e-08 | 5.97e-09 | 9.86e-09 | xfail |
| smoothing | `lambda2gcv_weather_log10lambda_4-gcv` | `gcv_curve` | 1.00e-08 | 1.26e-10 | 1.26e-10 | pass |
| smoothing | `lambda2gcv_weather_log10lambda_6-gcv` | `gcv_curve` | 1.00e-08 | 4.98e-10 | 4.01e-10 | pass |
| smoothing | `lambda2df_weather_lambda0.01-df` | `lambda_to_df` | 1.00e-08 | 1.42e-14 | 2.19e-16 | pass |
| smoothing | `lambda2df_weather_lambda1-df` | `lambda_to_df` | 1.00e-08 | 5.68e-14 | 8.78e-16 | pass |
| smoothing | `lambda2df_weather_lambda100-df` | `lambda_to_df` | 1.00e-08 | 6.32e-13 | 1.20e-14 | pass |
| smoothing | `df2lambda_weather_df10-lambda` | `df_to_lambda` | 1.00e-08 | 2.73e+03 | 7.70e-04 | xfail |
| smoothing | `df2lambda_weather_df20-lambda` | `df_to_lambda` | 1.00e-08 | 5.08e+01 | 9.83e-04 | xfail |
| smoothing | `df2lambda_weather_df30-lambda` | `df_to_lambda` | 1.00e-08 | 3.77e+00 | 8.56e-04 | xfail |
| smoothing | `lambda2gcv_weather_grid_mean-mean_gcv` | `gcv_curve` | 1.00e-08 | 1.20e-08 | 1.46e-08 | xfail |
| smoothing | `lambda2gcv_weather_grid_mean-argmin_log10lambda` | `gcv_curve` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basisPar_weather_fourier-coefs` | `smooth` | 1.00e-08 | 3.98e-13 | 1.26e-15 | pass |
| smoothing | `smooth_basisPar_weather_fourier-df` | `smooth` | 1.00e-08 | 2.13e-14 | 5.51e-16 | pass |
| smoothing | `smooth_basisPar_weather_fourier-gcv` | `smooth` | 1.00e-08 | 2.11e-15 | 2.77e-15 | pass |
| smoothing | `smooth_basisPar_weather_fourier-SSE` | `smooth` | 1.00e-08 | 4.55e-12 | 1.03e-15 | pass |
| smoothing | `smooth_basisPar_weather_fourier-penmat` | `smooth` | 1.00e-08 | 2.78e-17 | 3.01e-16 | pass |
| smoothing | `smooth_basisPar_weather_fourier-y2cMap` | `smooth` | 1.00e-08 | 9.66e-16 | 1.31e-14 | pass |
| smoothing | `smooth_basisPar_weather_fourier-fitted` | `smooth` | 1.00e-08 | 1.24e-13 | 3.73e-15 | pass |
| smoothing | `smooth_basisPar_growth_hgtm-coefs` | `smooth` | 1.00e-08 | 1.09e-05 | 5.63e-08 | xfail |
| smoothing | `smooth_basisPar_growth_hgtm-df` | `smooth` | 1.00e-08 | 1.07e-07 | 1.33e-08 | xfail |
| smoothing | `smooth_basisPar_growth_hgtm-gcv` | `smooth` | 1.00e-08 | 6.25e-07 | 4.93e-07 | xfail |
| smoothing | `smooth_basisPar_growth_hgtm-SSE` | `smooth` | 1.00e-08 | 2.03e-04 | 5.09e-07 | xfail |
| smoothing | `smooth_basisPar_growth_hgtm-penmat` | `smooth` | 1.00e-08 | 5.96e-08 | 2.23e-16 | pass |
| smoothing | `smooth_basisPar_growth_hgtm-y2cMap` | `smooth` | 1.00e-08 | 6.47e-08 | 8.52e-08 | xfail |
| smoothing | `smooth_basisPar_growth_hgtm-fitted` | `smooth` | 1.00e-08 | 1.09e-05 | 5.63e-08 | xfail |
| smoothing | `smooth_basisPar_synthetic-coefs` | `smooth` | 1.00e-08 | 1.55e-15 | 1.38e-15 | pass |
| smoothing | `smooth_basisPar_synthetic-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basisPar_synthetic-gcv` | `smooth` | 1.00e-08 | 1.39e-17 | 1.66e-16 | pass |
| smoothing | `smooth_basisPar_synthetic-SSE` | `smooth` | 1.00e-08 | 2.22e-16 | 1.25e-16 | pass |
| smoothing | `smooth_basisPar_synthetic-penmat` | `smooth` | 1.00e-08 | 7.11e-15 | 4.06e-16 | pass |
| smoothing | `smooth_basisPar_synthetic-y2cMap` | `smooth` | 1.00e-08 | 6.66e-16 | 8.55e-16 | pass |
| smoothing | `smooth_basisPar_synthetic-fitted` | `smooth` | 1.00e-08 | 7.77e-16 | 7.98e-16 | pass |
| smoothing | `data2fd_synthetic_bspline-coefs` | `smooth` | 1.00e-08 | 8.88e-16 | 7.95e-16 | pass |
| smoothing | `data2fd_weather_fourier-coefs` | `smooth` | 1.00e-08 | 6.33e-07 | 2.01e-09 | xfail |
| smoothing | `data2fd_growth_hgtf-coefs` | `smooth` | 1.00e-08 | 1.06e-05 | 5.80e-08 | xfail |
| smoothing | `smooth_monotone_growth_hgtf_3girls-Wfdobj_coefs` | `smooth` | 1.00e-05 | 8.72e-02 | 8.92e-03 | xfail |
| smoothing | `smooth_monotone_growth_hgtf_3girls-fitted` | `smooth` | 1.00e-05 | 5.00e-05 | 3.01e-07 | pass |
| smoothing | `smooth_monotone_growth_hgtf_3girls-beta_intercept` | `smooth` | 1.00e-05 | 5.00e-05 | 6.34e-07 | pass |
| smoothing | `smooth_monotone_growth_hgtf_3girls-beta_slope` | `smooth` | 1.00e-05 | 3.98e-01 | 7.16e-02 | xfail |
| smoothing | `smooth_monotone_growth_hgtf_3girls-deriv1` | `smooth` | 1.00e-05 | 2.67e-01 | 6.77e-02 | xfail |
| smoothing | `smooth_monotone_synthetic_sigmoid-Wfdobj_coefs` | `smooth` | 1.00e-05 | 1.99e-02 | 2.33e-02 | xfail |
| smoothing | `smooth_monotone_synthetic_sigmoid-fitted` | `smooth` | 1.00e-05 | 5.76e-06 | 1.13e-06 | pass |
| smoothing | `smooth_monotone_synthetic_sigmoid-beta_intercept` | `smooth` | 1.00e-05 | 5.76e-06 | 3.02e-06 | pass |
| smoothing | `smooth_monotone_synthetic_sigmoid-beta_slope` | `smooth` | 1.00e-05 | 6.34e-02 | 2.00e-02 | xfail |
| smoothing | `smooth_monotone_synthetic_sigmoid-deriv1` | `smooth` | 1.00e-05 | 3.08e-02 | 1.96e-02 | xfail |
| smoothing | `smooth_pos_synthetic_2curves-Wfdobj_coefs` | `smooth` | 1.00e-05 | 3.10e-05 | 6.49e-05 | xfail |
| smoothing | `smooth_pos_synthetic_2curves-fitted` | `smooth` | 1.00e-05 | 2.82e-05 | 1.75e-05 | xfail |
| smoothing | `smooth_pos_weather_precip_pr_rupert-Wfdobj_coefs` | `smooth` | 1.00e-05 | 1.82e-06 | 5.03e-08 | xfail |
| smoothing | `smooth_pos_weather_precip_pr_rupert-fitted` | `smooth` | 1.00e-05 | 4.85e-06 | 3.78e-07 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-coefs` | `smooth` | 1.00e-08 | 4.27e-15 | 1.07e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-df` | `smooth` | 1.00e-08 | 8.88e-16 | 1.17e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-gcv` | `smooth` | 1.00e-08 | 1.73e-17 | 1.37e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-SSE` | `smooth` | 1.00e-08 | 1.25e-16 | 1.28e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-penmat` | `smooth` | 1.00e-08 | 1.33e-15 | 4.44e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-y2cMap` | `smooth` | 1.00e-08 | 7.77e-16 | 5.68e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case1_n20-fitted` | `smooth` | 1.00e-08 | 3.11e-15 | 7.80e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-coefs` | `smooth` | 1.00e-08 | 4.00e-15 | 1.12e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-df` | `smooth` | 1.00e-08 | 1.78e-15 | 2.68e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-gcv` | `smooth` | 1.00e-08 | 2.26e-17 | 2.17e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-SSE` | `smooth` | 1.00e-08 | 1.39e-16 | 2.85e-15 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-penmat` | `smooth` | 1.00e-08 | 1.33e-15 | 4.44e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-y2cMap` | `smooth` | 1.00e-08 | 2.00e-15 | 7.88e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case2_n15-fitted` | `smooth` | 1.00e-08 | 1.33e-15 | 4.29e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-coefs` | `smooth` | 1.00e-08 | 1.78e-15 | 4.45e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-df` | `smooth` | 1.00e-08 | 0 | 0 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-gcv` | `smooth` | 1.00e-08 | 1.39e-17 | 7.78e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-SSE` | `smooth` | 1.00e-08 | 8.33e-17 | 3.89e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-penmat` | `smooth` | 1.00e-08 | 1.33e-15 | 4.44e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-y2cMap` | `smooth` | 1.00e-08 | 7.77e-16 | 6.24e-16 | pass |
| smoothing | `smooth_basis_irregular_argvals_case3_n25-fitted` | `smooth` | 1.00e-08 | 6.66e-16 | 1.69e-16 | pass |
