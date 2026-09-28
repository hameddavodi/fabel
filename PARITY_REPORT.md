# Parity report

fdatools 1.0.0 against R `fda` 6.3.0 (R version 4.6.1 (2026-06-24)).
Generated 2026-09-28 by `tools/parity_report.py` from a live run of `tests/parity`.

Errors are measured on every `assert_allclose` a parity test makes:
absolute error is `max |fdatools - R|`, relative error is the normwise
`max |fdatools - R| / max |R|` (entrywise ratios are meaningless on matrix entries
that are zero up to rounding). The summary columns show
the worst value over **passing** cases only; strict xfails (cases where R is
demonstrably the less accurate side) are listed with their measured reason below.

**Totals:** 935 parity checks -- 724 pass, 211 strict xfail (R defect), 0 fail, 0 skipped.

## By public symbol

| Symbol | R counterpart | Golden cases | Checks | Max abs err | Max rel err | Status |
|---|---|---:|---:|---:|---:|---|
| `FData` | fd, eval.fd, mean.fd, sd.fd, center.fd, deriv.fd, +.fd | 17 | 22 | 3.41e-12 | 1.77e-15 | pass, 3 xfail (R defect) |
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
| `SmoothResult` | smooth.basis return list | 2 | 19 | 6.77e-12 | 3.66e-13 | pass, 4 xfail (R defect) |
| `Smoother` | (new: sklearn estimator) | 0 | 0 | — | — | no golden cases |
| `gcv_curve` | lambda2gcv | 6 | 7 | 4.98e-10 | 4.01e-10 | pass, 2 xfail (R defect) |
| `lambda_to_df` | lambda2df | 3 | 3 | 6.32e-13 | 1.20e-14 | pass |
| `df_to_lambda` | df2lambda | 3 | 3 | — | — | all 3 checks xfail (R defect) |
| `FPCA` | pca.fd, varmx.pca.fd | 10 | 55 | 5.36e-08 | 1.23e-09 | pass, 32 xfail (R defect) |
| `FCCA` | cca.fd | 2 | 10 | 3.93e-12 | 6.32e-13 | pass |
| `fregress` | fRegress, predict.fRegress, fRegress.stderr, fRegress.CV | 9 | 60 | 8.86e-09 | 6.74e-09 | pass, 34 xfail (R defect) |
| `FRegress` | fRegress (estimator form) | 0 | 0 | — | — | no golden cases |
| `linmod` | linmod | 4 | 25 | 4.49e-11 | 4.42e-12 | pass, 6 xfail (R defect) |
| `register` | register.fd, AmpPhaseDecomp | 6 | 56 | 1.22e-04 | 1.93e-05 | pass, 31 xfail (R defect) |
| `landmark_register` | landmarkreg | 2 | 4 | 4.92e-07 | 5.89e-07 | pass, 3 xfail (R defect) |
| `Registrator` | register.fd (estimator form) | 0 | 0 | — | — | no golden cases |
| `PDA` | pda.fd, pda.overlay | 9 | 52 | 1.24e-09 | 1.32e-06 | pass |
| `PDAStability` | eigen.pda | 2 | 10 | 9.35e-12 | 1.23e-12 | pass, 4 xfail (R defect) |
| `phase_plane` | phaseplanePlot | 0 | 0 | — | — | no golden cases |
| `sparse_mean` | smooth.sparse.mean | 3 | 3 | 1.17e-12 | 6.77e-15 | pass |
| `PACE` | pcaPACE, covPACE, scoresPACE | 4 | 27 | 3.35e-10 | 4.01e-11 | pass, 11 xfail (R defect) |
| `fit_density` | density.fd (dropped in fda 6.3.0; checked via intensity.fd / n) | 5 | 36 | 1.33e-04 | 4.27e-06 | pass, 4 xfail (R defect) |
| `fit_intensity` | intensity.fd | 4 | 20 | 1.64e-04 | 2.79e-06 | pass, 1 xfail (R defect) |
| `profile_ode` | CSTR2in, CSTR2, CSTRfitLS, CSTRfn, CSTRres, CSTRsse, quadset | 21 | 21 | 2.04e-06 | 1.39e-06 | pass |
| `stats.cov` | var.fd | 1 | 1 | 4.41e-13 | 4.70e-15 | pass |
| `stats.cor` | cor.fd | 1 | 1 | 3.22e-15 | 4.03e-15 | pass |
| `stats.depth` | fdepth | 1 | 6 | 7.11e-15 | 4.01e-16 | pass |
| `stats.boxplot` | fbplot, boxplot.fd | 1 | 2 | 7.82e-14 | 1.38e-14 | pass |
| `stats.f_test` | Fperm.fd | 1 | 5 | 0 | 0 | pass, 3 xfail (R defect) |
| `stats.t_test` | tperm.fd | 1 | 5 | 3.38e-14 | 1.44e-14 | pass |
| `stats.confidence_band` | smooth.basis y2cMap variance, fRegress.stderr | 4 | 24 | 1.48e-12 | 1.22e-13 | pass, 6 xfail (R defect) |
| `datasets.load_*` | data(package='fda') | 11 | 12 | 3.21e-09 | 5.94e-14 | pass |
| `nn.BasisLayer` | (new: PyTorch layer) | 0 | 0 | — | — | no golden cases |
| `nn.FDataDataset` | (new: PyTorch dataset) | 0 | 0 | — | — | no golden cases |
| `from_pandas` | (new) | 0 | 0 | — | — | no golden cases |
| `to_pandas` | (new) | 0 | 0 | — | — | no golden cases |
| `to_xarray` | (new) | 0 | 0 | — | — | no golden cases |
| `read_rds` | readRDS on fd / bifd / basisfd objects | 0 | 2 | 2.22e-16 | 9.98e-17 | pass |

## Strict xfails (R fda is the less accurate side)

- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-beta1]` (bands, rel err 2.15e-04): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-df]` (bands, rel err 3.11e-08): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-sigma2]` (bands, rel err 1.80e-05): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-bvar]` (bands, rel err 9.20e-05): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-stderr0]` (bands, rel err 3.49e-06): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-stderr1]` (bands, rel err 4.93e-05): R integrates int temp_i(t) theta_k(t) dt approximately; fdatools' integrals are exact. R's own dense trapezoid rule (365001 points) gives int temp_1 theta_4 = 943.854155886, fdatools' value to 1e-13, while R's inprod() gives 947.596583972 (4.0e-3 off); fRegress's integrals carry the same kind of error (see tests/parity/test_regression.py). Downstream, measured relative to the largest entry: beta1 2.2e-4, df 3.1e-8, sigma2 1.8e-5, bvar 9.2e-5, the intercept standard error 3.5e-6 and the pointwise standard error of the temperature coefficient 4.9e-5.
- `bspline_penalty_k4_n4_dom0_1_L0` (basis, rel err 8.57e-01): With no interior knots R returns the monomial Gram (the Hilbert matrix) instead of the Bernstein-basis Gram; fdatools' [0, 0] entry is 1/7.
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
- `inprod_fourier_L0_0` (core, rel err 4.26e-05): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. fdatools uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_fourier_L1_1` (core, rel err 1.17e-04): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. fdatools uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_fourier_L2_2` (core, rel err 1.99e-04): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. fdatools uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `inprod_basis_bspline_x_fourier` (core, rel err 5.73e-05): R's inprod integrates by Romberg iteration with a convergence tolerance of 1e-4 and returns values ~1.3e-4 away from the exact integral. Its public API exposes no tolerance argument and no exact alternative outside inprod.bspline, which does not accept a Fourier basis. fdatools uses Gauss-Legendre panels that are exact for these integrands -- its Gram matrices match R's own eval.penalty to 1e-8.
- `deriv_fd_bspline_order6_L2` (core, rel err 1.27e-02): deriv.fd re-expands D^2 x in the *original* order-6 basis, which holds only C^4 functions while D^2 x is C^2; the result is 1.3% off (max 59.9 on a curve of size 4655). fdatools returns the exact order-4 spline.
- `fd_mul_bspline` (core, rel err 1.31e-01): times.fd keeps the factors' break points, so the order-7 product space is C^5 where the true product is only C^2 and cannot be represented. R's answer is 12.8% away from the exact product (max 0.254 on a curve of size 1.99), and is not even the L2 projection onto its own basis. fdatools raises the interior knot multiplicities and is exact.
- `fd_power2` (core, rel err 5.14e-05): ^.fd projects onto an arbitrary uniform refinement (80 intervals, order 7) whose knots miss the curve's own breaks at k/7, leaving a 5.7e-5 relative error. fdatools squares the curve exactly.
- `pca_fd_weather_nharm2_harmlambda0-harmonics` (decomposition, rel err 3.52e-09): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_weather_nharm2_harmlambda0-scores` (decomposition, rel err 1.91e-03): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `pca_fd_weather_nharm2_harmlambda10000-harmonics` (decomposition, rel err 2.51e-11): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_weather_nharm2_harmlambda10000-scores` (decomposition, rel err 7.69e-03): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `pca_fd_weather_nharm4_harmlambda0-harmonics` (decomposition, rel err 5.51e-09): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_weather_nharm4_harmlambda0-scores` (decomposition, rel err 2.23e-03): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `pca_fd_weather_nharm4_harmlambda10000-harmonics` (decomposition, rel err 3.76e-11): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_weather_nharm4_harmlambda10000-scores` (decomposition, rel err 1.74e-03): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-values` (decomposition, rel err 4.41e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-harmonics` (decomposition, rel err 1.11e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-scores` (decomposition, rel err 3.91e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-varprop` (decomposition, rel err 4.10e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-rotmat` (decomposition, rel err 1.10e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-values` (decomposition, rel err 2.03e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-harmonics` (decomposition, rel err 4.85e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-scores` (decomposition, rel err 2.24e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-varprop` (decomposition, rel err 2.20e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-rotmat` (decomposition, rel err 5.41e-03): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `pca_fd_growth_hgtm_nharm3-harmonics` (decomposition, rel err 4.57e-05): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_growth_hgtm_nharm3-scores` (decomposition, rel err 4.51e-05): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `varmx_pca_fd_growth_hgtm_nharm3-values` (decomposition, rel err 7.21e-05): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_growth_hgtm_nharm3-harmonics` (decomposition, rel err 1.21e-04): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_growth_hgtm_nharm3-scores` (decomposition, rel err 6.90e-05): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_growth_hgtm_nharm3-varprop` (decomposition, rel err 6.13e-05): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `varmx_pca_fd_growth_hgtm_nharm3-rotmat` (decomposition, rel err 1.23e-04): R's varmx() stops before the varimax criterion is stationary and returns a rotation that is not orthogonal. On R's own 501-point harmonic values the gradient of the criterion on the rotation group (the skew part of L' dV/dL) is 4.0e-5 / 2.2e-4 / 6.7e-3 at R's rotation (weather lambda 0 / weather lambda 1e4 / growth) against 8.2e-14 / 7.7e-14 / 3.0e-11 at fdatools', and max\|T'T - I\| is 1.8e-7 / 1.3e-7 / 9.6e-8 for R's T. On the weather cases R's criterion is also lower (0.04072609184 against 0.04072610195; 0.04092609258 against 0.04092669536). Growth inherits the inprod() Gram error of its source pca.fd case as well.
- `pca_fd_synthetic_bspline_n10curves-values` (decomposition, rel err 5.30e-05): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_synthetic_bspline_n10curves-harmonics` (decomposition, rel err 4.92e-04): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_synthetic_bspline_n10curves-scores` (decomposition, rel err 1.90e-04): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `pca_fd_synthetic_bspline_n10curves-varprop` (decomposition, rel err 3.17e-05): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_gait_multivariate_nharm3-values` (decomposition, rel err 7.56e-10): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_gait_multivariate_nharm3-harmonics` (decomposition, rel err 6.87e-08): R's pca.fd builds the covariance operator W V W with the Gram matrix W from inprod(), whose Romberg quadrature is only good to 4-5 digits, while it normalises with the exact Gram matrix. inprod(basis, basis) differs from R's own exact eval.penalty(basis, 0) by 9.1e-5 (absolute) on the synthetic order-4 spline, 2.6e-4 (relative) on the growth order-6 spline and 1.38e-6 on both Fourier bases (weather, gait), where the error sits on the highest-frequency pair. Substituting R's inprod() Gram into fdatools' eigenproblem reproduces R's values and harmonics to 1e-15 (lambda = 0) and 1.9e-11 (weather lambda = 1e4, which also carries R's harmonic accelerator penalty error, see test_smoothing.py). With the exact Gram the gaps are 5.3e-5 (synthetic values), 1.3e-6 (the small gait eigenvalues), 4.9e-4 / 4.6e-5 (synthetic / growth harmonics), 1.38e-6 (weather harmonics, rows 63-64) and 6.9e-8 of the largest coefficient (gait harmonics, spread over the coupled hip/knee blocks). fdatools uses the exact Gram matrix.
- `pca_fd_gait_multivariate_nharm3-scores` (decomposition, rel err 6.51e-04): R's pca.fd scores are inprod(centred curves, harmonics), a Romberg quadrature good to 4-5 digits, not the exact C' W h. Measured inside R: inprod() differs from t(C) %*% eval.penalty(basis, 0) %*% h by 0.660 on weather scores of size ~150 and by 2.2e-3 on growth scores of size 48.7. fdatools' scores are the exact inner products; the relative gaps asserted here are 1.7e-3 to 7.7e-3 (weather), 4.5e-5 (growth), 1.9e-4 (synthetic) and 6.5e-4 (gait).
- `intensity_sine_bspline11_order3_L1_lam0-coefs` (density, rel err 2.17e-06): R's intensity.fd integrates exp(W) with a quadrature that is off by up to 1.6e-4 in its reported criterion f; its coefficients are therefore not stationary for the exact criterion, and fdatools' (exact Gauss-Legendre integral, gradient norm < 2e-13) give a lower exact criterion in every case (test_fdatools_improves_on_r_criterion). Exact gradient norm at R's coefficients 4.07e-5. Worst entry: coefficient 7, R 0.11139673524, fdatools 0.11140196492 (abs 5.2e-6, rel 4.7e-5).
- `density_regina_bspline13_L2_lam0.1-coefs` (density, rel err 2.64e-06): R's intensity.fd integrates exp(W) with a quadrature that is off by up to 1.6e-4 in its reported criterion f; its coefficients are therefore not stationary for the exact criterion, and fdatools' (exact Gauss-Legendre integral, gradient norm < 2e-13) give a lower exact criterion in every case (test_fdatools_improves_on_r_criterion). Exact gradient norm at R's coefficients 4.18e-5. Worst entry: coefficient 6, R -0.037409071057, fdatools -0.037405675073 (abs 3.4e-6, rel 9.1e-5).
- `density_regina_bspline13_L2_lam0.1-log_intensity` (density, rel err 1.24e-06): R's intensity.fd integrates exp(W) with a quadrature that is off by up to 1.6e-4 in its reported criterion f; its coefficients are therefore not stationary for the exact criterion, and fdatools' (exact Gauss-Legendre integral, gradient norm < 2e-13) give a lower exact criterion in every case (test_fdatools_improves_on_r_criterion). Worst entry: W(42.85), R 0.033795023457, fdatools 0.033793807749 (abs 1.2e-6, rel 3.6e-5); the intensity itself agrees to 4.3e-6.
- `density_gamma_bspline9_L2_lam1-log_intensity` (density, rel err 1.38e-07): R's intensity.fd integrates exp(W) with a quadrature that is off by up to 1.6e-4 in its reported criterion f; its coefficients are therefore not stationary for the exact criterion, and fdatools' (exact Gauss-Legendre integral, gradient norm < 2e-13) give a lower exact criterion in every case (test_fdatools_improves_on_r_criterion). Exact gradient norm at R's coefficients 1.60e-5. Worst entry: W(8.7), R 0.010768158864, fdatools 0.010768044513 (abs 1.1e-7, rel 1.1e-5).
- `density_truncnormal_monomial3_L1_lam0-log_intensity` (density, rel err 3.07e-07): R's intensity.fd integrates exp(W) with a quadrature that is off by up to 1.6e-4 in its reported criterion f; its coefficients are therefore not stationary for the exact criterion, and fdatools' (exact Gauss-Legendre integral, gradient norm < 2e-13) give a lower exact criterion in every case (test_fdatools_improves_on_r_criterion). Exact gradient norm at R's coefficients 5.16e-5. Worst entry: W(-2.1), R -0.033503115309, fdatools -0.033503752936 (abs 6.4e-7, rel 1.9e-5).
- `test_golden_field[linmod_weather_logprecip_on_temp-beta0estfd_coefs]` (linmod, rel err 4.68e-07): R's integrals over the Fourier bases are approximate. Its harmonic-accelerator penalty of the Fourier(11) coefficient basis (eval.penalty) is 1.2e-4 (relative) off the closed form k^2 w^6 (k^2 - 1)^2, and its Fourier Gram inprod(fb11, fb11) 1.4e-6 off the identity. fdatools' exact fit is 4.7e-7 (intercept), 1.9e-6 (surface) and 3.7e-7 (fitted curves) off R, relative to the largest entry.
- `test_golden_field[linmod_weather_logprecip_on_temp-beta1estbifd_coefs]` (linmod, rel err 1.88e-06): R's integrals over the Fourier bases are approximate. Its harmonic-accelerator penalty of the Fourier(11) coefficient basis (eval.penalty) is 1.2e-4 (relative) off the closed form k^2 w^6 (k^2 - 1)^2, and its Fourier Gram inprod(fb11, fb11) 1.4e-6 off the identity. fdatools' exact fit is 4.7e-7 (intercept), 1.9e-6 (surface) and 3.7e-7 (fitted curves) off R, relative to the largest entry.
- `test_golden_field[linmod_weather_logprecip_on_temp-yhatfdobj_coefs]` (linmod, rel err 3.72e-07): R's integrals over the Fourier bases are approximate. Its harmonic-accelerator penalty of the Fourier(11) coefficient basis (eval.penalty) is 1.2e-4 (relative) off the closed form k^2 w^6 (k^2 - 1)^2, and its Fourier Gram inprod(fb11, fb11) 1.4e-6 off the identity. fdatools' exact fit is 4.7e-7 (intercept), 1.9e-6 (surface) and 3.7e-7 (fitted curves) off R, relative to the largest entry.
- `test_golden_field[linmod_synthetic_bspline-beta0estfd_coefs]` (linmod, rel err 5.91e-05): R's inprod of cubic B-spline bases is 2.2e-5 (Gram) and 1.1e-4 (the covariate integrals int x_i theta_s) off the exact values; fdatools' exact intercept and surface are 5.9e-5 and 9.6e-5 off R. With R's own matrices (r_integrals) fdatools reproduces R's coefficients to 1e-15 (test_normal_equations_with_r_integrals). R's fitted curves are further a 201-point least-squares fit, not the L2 projection: 1.4e-3 off fdatools (test_r_fitted_curves_are_a_grid_least_squares_fit).
- `test_golden_field[linmod_synthetic_bspline-beta1estbifd_coefs]` (linmod, rel err 9.55e-05): R's inprod of cubic B-spline bases is 2.2e-5 (Gram) and 1.1e-4 (the covariate integrals int x_i theta_s) off the exact values; fdatools' exact intercept and surface are 5.9e-5 and 9.6e-5 off R. With R's own matrices (r_integrals) fdatools reproduces R's coefficients to 1e-15 (test_normal_equations_with_r_integrals). R's fitted curves are further a 201-point least-squares fit, not the L2 projection: 1.4e-3 off fdatools (test_r_fitted_curves_are_a_grid_least_squares_fit).
- `test_golden_field[linmod_synthetic_bspline-yhatfdobj_coefs]` (linmod, rel err 1.38e-03): R's inprod of cubic B-spline bases is 2.2e-5 (Gram) and 1.1e-4 (the covariate integrals int x_i theta_s) off the exact values; fdatools' exact intercept and surface are 5.9e-5 and 9.6e-5 off R. With R's own matrices (r_integrals) fdatools reproduces R's coefficients to 1e-15 (test_normal_equations_with_r_integrals). R's fitted curves are further a 201-point least-squares fit, not the L2 projection: 1.4e-3 off fdatools (test_r_fitted_curves_are_a_grid_least_squares_fit).
- `test_monotone_derivatives_match_eval_monfd[0-monfd_growth_girl1]` (monotone_derivs, rel err 2.44e-06): R's eval.monfd(t, W, 0) integrates exp W numerically: measured against scipy.integrate.quad (tol 1e-14) R is off by 2.4e-6 (growth) / 1.8e-6 (fixed W) of the range, 9e-5 relative pointwise; fdatools by at most 4.3e-14
- `test_monotone_derivatives_match_eval_monfd[0-monfd_posfd_fixed_w_two_curves]` (monotone_derivs, rel err 1.81e-06): R's eval.monfd(t, W, 0) integrates exp W numerically: measured against scipy.integrate.quad (tol 1e-14) R is off by 2.4e-6 (growth) / 1.8e-6 (fixed W) of the range, 9e-5 relative pointwise; fdatools by at most 4.3e-14
- `test_positive_derivatives_match_eval_posfd[2]` (monotone_derivs, rel err 8.66e-01): R's eval.posfd(t, W, 2) returns exp(W) D2W (to 1.5e-16), not D2 exp W = exp(W) (D2W + (DW)^2); measured 0.87 relative difference
- `test_scaled_monotone_derivatives_match_predict[0]` (monotone_derivs, rel err 1.26e-06): R's eval.monfd(t, W, 0) integrates exp W numerically: measured against scipy.integrate.quad (tol 1e-14) R is off by 2.4e-6 (growth) / 1.8e-6 (fixed W) of the range, 9e-5 relative pointwise; fdatools by at most 4.3e-14
- `pace_growth_bspline6_covlambda0-values` (pace, rel err 7.21e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_growth_bspline6_covlambda0-varprop` (pace, rel err 6.78e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_growth_bspline6_covlambda0-harmonics` (pace, rel err 9.03e-05): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_growth_bspline6_covlambda10-cov_coefs` (pace, rel err 1.16e-05): R's covPACE penalty is lambda (W (x) P + P (x) W) with W = inprod(basis, basis) (Romberg, 4.0e-5 off the exact Gram for this order-4 spline on [1, 18]) and the exact P = eval.penalty(basis, 2). With R's inprod() W (recorded in the golden file) fdatools reproduces R's surface to 1.9e-14 (test_same_math_with_r_inprod_gram); with the exact W the relative gap is 1.2e-5. Unpenalised surfaces (cov lambda 0) match to 6e-14.
- `pace_growth_bspline6_covlambda10-values` (pace, rel err 4.48e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_growth_bspline6_covlambda10-varprop` (pace, rel err 4.05e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_growth_bspline6_covlambda10-harmonics` (pace, rel err 5.86e-05): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_weather_fourier5_bspline7-values` (pace, rel err 8.95e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_weather_fourier5_bspline7-varprop` (pace, rel err 7.03e-06): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_weather_fourier5_bspline7-harmonics` (pace, rel err 1.59e-04): R's pcaPACE builds the covariance operator J C J' with J = inprod(harmonic basis, covariance basis), a Romberg quadrature good to 4-5 digits, and normalises with the exact eval.penalty Gram. Measured max\|inprod(b, b) - eval.penalty(b, 0)\|: 4.0e-5 for the order-4 spline with 6 functions on [1, 18], 2.9e-4 for the one with 7 functions on [1, 365]; it is exactly 0 for monomial bases. With R's own inprod() J (recorded in the golden file) fdatools' eigenproblem reproduces R's values to 1.5e-15 and harmonics to 1.7e-14 (test_same_math_with_r_inprod_gram). With the exact J the relative gaps are 7.2e-6 / 4.5e-6 / 9.0e-6 (values) and 9.0e-5 / 5.9e-5 / 1.6e-4 (harmonics) for growth cov-lambda 0 / growth cov-lambda 10 / weather. fdatools uses the exact J.
- `pace_weather_fourier5_bspline7-scores` (pace, rel err 4.52e+00): R's scoresPACE does not compute the conditional expectation Lambda Xi_i' Sigma_i^-1 (y_i - mu_i). Measured as a black box: every finite row of R's score matrix is exactly parallel (to 1e-12) to lambda_k * xi_k(tau) for ONE integer day tau (test_r_scores_use_one_point), so each curve's score vector has rank one in the harmonics, and 33 of the 35 rows are NA. On a dense regular probe (20 points, all curves at t = 1..20) curve i used tau = i, and a unit change of the first observation moved the scores about 365 times more than the same change of the fifth; changing another curve's data moved them too. The same function errors ('evalarg contains 1 NA') for non-integer times such as growth ages. fdatools implements the published PACE estimator (Yao, Mueller & Wang 2005).
- `test_limits_match_r_limvals[pda_forced_refinery_constant]` (pda_forcing, rel err 2.00e+00): R 6.3.0 eigen.pda limvals is not the equilibrium -A(t)^-1 f(t): measured on refinery (one first-order equation) it is exactly -z* (sign flipped); for the second-order equation it is (0, -a/b0) where z* = (a/b0, 0); for the forced two-equation system it is (-0.5000000338, -1.9e-8) where z* = (0.4999999929, -0.3861262521) -- the second equation's forcing is lost. fdatools returns z*, checked against the ODE solution in tests/unit/test_dynamics.py.
- `test_limits_match_r_limvals[pda_forced_refinery_two_forcings]` (pda_forcing, rel err 2.00e+00): R 6.3.0 eigen.pda limvals is not the equilibrium -A(t)^-1 f(t): measured on refinery (one first-order equation) it is exactly -z* (sign flipped); for the second-order equation it is (0, -a/b0) where z* = (a/b0, 0); for the forced two-equation system it is (-0.5000000338, -1.9e-8) where z* = (0.4999999929, -0.3861262521) -- the second equation's forcing is lost. fdatools returns z*, checked against the ODE solution in tests/unit/test_dynamics.py.
- `test_limits_match_r_limvals[pda_forced_order2_constant]` (pda_forcing, rel err 1.00e+00): R 6.3.0 eigen.pda limvals is not the equilibrium -A(t)^-1 f(t): measured on refinery (one first-order equation) it is exactly -z* (sign flipped); for the second-order equation it is (0, -a/b0) where z* = (a/b0, 0); for the forced two-equation system it is (-0.5000000338, -1.9e-8) where z* = (0.4999999929, -0.3861262521) -- the second equation's forcing is lost. fdatools returns z*, checked against the ODE solution in tests/unit/test_dynamics.py.
- `test_limits_match_r_limvals[pda_forced_system_order1]` (pda_forcing, rel err 2.00e+00): R 6.3.0 eigen.pda limvals is not the equilibrium -A(t)^-1 f(t): measured on refinery (one first-order equation) it is exactly -z* (sign flipped); for the second-order equation it is (0, -a/b0) where z* = (a/b0, 0); for the forced two-equation system it is (-0.5000000338, -1.9e-8) where z* = (0.4999999929, -0.3861262521) -- the second equation's forcing is lost. fdatools returns z*, checked against the ODE solution in tests/unit/test_dynamics.py.
- `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-regfd_coefs]` (registration, rel err 9.77e-01): landmarkreg's regfd is not x(h(t)): for girl 3 it has a coefficient of 3360.37 and a value of 421.27 cm at age 1.2 (the curve spans 67.6-183.2 cm), where x(h(1.2)) = 80.25 cm. It is sampled through R's faulty inverse warp (above), whose steep start leaves the first years almost unsampled. fdatools projects x(h(t)) on the registration grid.
- `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-warpfd_coefs]` (registration, rel err 4.23e-06): landmarkreg's warpfd is not the exact warp of its own Wfd. Integrating exp(W) for R's recorded Wfd to rounding error gives warps that differ from R's warpfd by up to 3.5e-5 years; the order-6 spline coefficients amplify that to 7.6e-5, so 6 of 350 coefficients miss 1e-5 relative (worst 6.9e-5). No 1025-point trapezoid or grid choice reproduces R's values, so the gap is R's quadrature, not a different warp: fdatools' own Wfd agrees with R's to 5.9e-7.
- `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-warpinvfd_coefs]` (registration, rel err 9.89e-02): landmarkreg's warpinvfd is not the inverse of its warpfd: warpinvfd(warpfd(t)) - t reaches 1.73 years (girl 3, t = 2.48), and warpinvfd rises from 1 to 3.53 over the first year for girl 1, where the true inverse reaches 2.17. fdatools inverts its warps by Newton's method; its inverse composes with the warp to the identity within 2e-8.
- `test_registration_matches_r[register_fd_growth_hgtf_to_mean-regfd_coefs]` (registration, rel err 1.60e-03): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's own discretisation (grid mean over 351 points, crit=2, lambda=1), R's Wfd scores 0.621731 / 0.304731 / 0.912346 / 0.114607 / 0.310316 / 0.070639 / 0.815999 / 0.194697 / 0.282611 / 0.646418 for the ten girls, with gradients of max-norm 0.021 to 0.371 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.619192 / 0.304554 / 0.911651 / 0.114439 / 0.310197 / 0.070245 / 0.808155 / 0.194462 / 0.282337 / 0.639359 -- lower on every curve. R's own dbglev=2 trace agrees: it stops girls 1-4 with gradient lengths 0.0516 / 0.1017 / 0.028 / 0.0081 after 'Criterion increased, terminating iterations' and 'Reset twice, terminating'. The line search fails because R's criterion is not a function of the warp it returns: at a constant W (h(t) = t exactly) with x0 = x = 1000 t on [1, 18] it reports 1.1564 instead of 0 and a gradient of 3026; the moments of its criterion put its internal warp at about (1 - 1.0e-4) t there. The optima differ by up to 0.089 in W.
- `test_registration_matches_r[register_fd_growth_hgtf_to_mean-warpfd_coefs]` (registration, rel err 7.49e-03): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's own discretisation (grid mean over 351 points, crit=2, lambda=1), R's Wfd scores 0.621731 / 0.304731 / 0.912346 / 0.114607 / 0.310316 / 0.070639 / 0.815999 / 0.194697 / 0.282611 / 0.646418 for the ten girls, with gradients of max-norm 0.021 to 0.371 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.619192 / 0.304554 / 0.911651 / 0.114439 / 0.310197 / 0.070245 / 0.808155 / 0.194462 / 0.282337 / 0.639359 -- lower on every curve. R's own dbglev=2 trace agrees: it stops girls 1-4 with gradient lengths 0.0516 / 0.1017 / 0.028 / 0.0081 after 'Criterion increased, terminating iterations' and 'Reset twice, terminating'. The line search fails because R's criterion is not a function of the warp it returns: at a constant W (h(t) = t exactly) with x0 = x = 1000 t on [1, 18] it reports 1.1564 instead of 0 and a gradient of 3026; the moments of its criterion put its internal warp at about (1 - 1.0e-4) t there. The optima differ by up to 0.089 in W.
- `test_registration_matches_r[register_fd_growth_hgtf_to_mean-Wfd_coefs]` (registration, rel err 5.91e-02): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's own discretisation (grid mean over 351 points, crit=2, lambda=1), R's Wfd scores 0.621731 / 0.304731 / 0.912346 / 0.114607 / 0.310316 / 0.070639 / 0.815999 / 0.194697 / 0.282611 / 0.646418 for the ten girls, with gradients of max-norm 0.021 to 0.371 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.619192 / 0.304554 / 0.911651 / 0.114439 / 0.310197 / 0.070245 / 0.808155 / 0.194462 / 0.282337 / 0.639359 -- lower on every curve. R's own dbglev=2 trace agrees: it stops girls 1-4 with gradient lengths 0.0516 / 0.1017 / 0.028 / 0.0081 after 'Criterion increased, terminating iterations' and 'Reset twice, terminating'. The line search fails because R's criterion is not a function of the warp it returns: at a constant W (h(t) = t exactly) with x0 = x = 1000 t on [1, 18] it reports 1.1564 instead of 0 and a gradient of 3026; the moments of its criterion put its internal warp at about (1 - 1.0e-4) t there. The optima differ by up to 0.089 in W.
- `test_registration_matches_r[register_fd_weather_periodic_crit1-regfd_coefs]` (registration, rel err 5.13e-02): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's discretisation (grid mean over 651 points, crit=1, lambda=1, periodic shift), R's (Wfd, shift) scores 0.634340 / 0.286287 / 0.101681 / 1.122814 / 0.806865 for the five stations, with gradients of max-norm 0.0088 to 0.205 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.631612 / 0.285765 / 0.101659 / 1.120003 / 0.798741 -- lower on every curve. Shifts and warps trade off along a shallow valley, so the stationary point is far from R's: shifts 6.269 / -4.475 / 2.127 / 8.781 / -16.230 days against R's 6.820 / -4.428 / 2.104 / 9.895 / -18.618, and W up to 4.1 apart. R's internal warp is also not the one it returns (see the growth case).
- `test_registration_matches_r[register_fd_weather_periodic_crit1-warpfd_coefs]` (registration, rel err 1.01e-01): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's discretisation (grid mean over 651 points, crit=1, lambda=1, periodic shift), R's (Wfd, shift) scores 0.634340 / 0.286287 / 0.101681 / 1.122814 / 0.806865 for the five stations, with gradients of max-norm 0.0088 to 0.205 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.631612 / 0.285765 / 0.101659 / 1.120003 / 0.798741 -- lower on every curve. Shifts and warps trade off along a shallow valley, so the stationary point is far from R's: shifts 6.269 / -4.475 / 2.127 / 8.781 / -16.230 days against R's 6.820 / -4.428 / 2.104 / 9.895 / -18.618, and W up to 4.1 apart. R's internal warp is also not the one it returns (see the growth case).
- `test_registration_matches_r[register_fd_weather_periodic_crit1-Wfd_coefs]` (registration, rel err 2.27e+00): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's discretisation (grid mean over 651 points, crit=1, lambda=1, periodic shift), R's (Wfd, shift) scores 0.634340 / 0.286287 / 0.101681 / 1.122814 / 0.806865 for the five stations, with gradients of max-norm 0.0088 to 0.205 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.631612 / 0.285765 / 0.101659 / 1.120003 / 0.798741 -- lower on every curve. Shifts and warps trade off along a shallow valley, so the stationary point is far from R's: shifts 6.269 / -4.475 / 2.127 / 8.781 / -16.230 days against R's 6.820 / -4.428 / 2.104 / 9.895 / -18.618, and W up to 4.1 apart. R's internal warp is also not the one it returns (see the growth case).
- `test_registration_matches_r[register_fd_weather_periodic_crit1-shift]` (registration, rel err 1.28e-01): register.fd returns a point that is not a minimum of its own criterion. Evaluated with R's discretisation (grid mean over 651 points, crit=1, lambda=1, periodic shift), R's (Wfd, shift) scores 0.634340 / 0.286287 / 0.101681 / 1.122814 / 0.806865 for the five stations, with gradients of max-norm 0.0088 to 0.205 there; fdatools' Newton iterate is stationary (gradient < 1e-10) at 0.631612 / 0.285765 / 0.101659 / 1.120003 / 0.798741 -- lower on every curve. Shifts and warps trade off along a shallow valley, so the stationary point is far from R's: shifts 6.269 / -4.475 / 2.127 / 8.781 / -16.230 days against R's 6.820 / -4.428 / 2.104 / 9.895 / -18.618, and W up to 4.1 apart. R's internal warp is also not the one it returns (see the growth case).
- `test_warps_from_r_latent_match_r[register_fd_growth_hgtf_to_mean-regfd_coefs]` (registration, rel err 3.58e-04): R's regfd is not the registered curve of the warp R returns. Given R's own Wfd, x(h(t)) projected on the curve basis differs from R's regfd by 3.6e-4 (growth) / 3.2e-4 (weather) relative. With x(t) = t on a basis that also holds the warp, R's regfd differs from R's warpfd by 9e-5 although both are the same function x(h(t)) = h(t); and R's criterion at a constant W (h(t) = t exactly) is not zero -- its internal warp is about (1 - 1.0e-4) t there, not t.
- `test_warps_from_r_latent_match_r[register_fd_weather_periodic_crit1-warpfd_coefs]` (registration, rel err 2.16e-06): R integrates exp(W) for warpfd with the trapezoidal rule on 1025 points (monfn) -- reproduced to 2e-12. For the weather warps that rule is off the exact integral by up to 4.7e-4 days, which moves the least-squares warp coefficients by up to 8.1e-4; the smallest coefficient (-0.1143) is then 2.5e-3 out relative. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_warps_from_r_latent_match_r[register_fd_weather_periodic_crit1-regfd_coefs]` (registration, rel err 3.19e-04): R's regfd is not the registered curve of the warp R returns. Given R's own Wfd, x(h(t)) projected on the curve basis differs from R's regfd by 3.6e-4 (growth) / 3.2e-4 (weather) relative. With x(t) = t on a basis that also holds the warp, R's regfd differs from R's warpfd by 9e-5 although both are the same function x(h(t)) = h(t); and R's criterion at a constant W (h(t) = t exactly) is not zero -- its internal warp is about (1 - 1.0e-4) t there, not t.
- `test_register_matches_r[register_fd_gait_hip_knee_crit2-regfd_coefs]` (registration_multivariate, rel err 1.80e-02): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=0.01, first variable), R's Wfd scores 4.005844 / 1.376484 / 1.724921 / 3.736952 / 7.200670 for the five boys, with gradients of max-norm 0.19 to 2.09 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 3.993467 / 1.369739 / 1.716418 / 3.722099 / 6.791805 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 2.01 in W.
- `test_register_matches_r[register_fd_gait_hip_knee_crit2-warpfd_coefs]` (registration_multivariate, rel err 4.17e-02): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=0.01, first variable), R's Wfd scores 4.005844 / 1.376484 / 1.724921 / 3.736952 / 7.200670 for the five boys, with gradients of max-norm 0.19 to 2.09 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 3.993467 / 1.369739 / 1.716418 / 3.722099 / 6.791805 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 2.01 in W.
- `test_register_matches_r[register_fd_gait_hip_knee_crit2-Wfd_coefs]` (registration_multivariate, rel err 4.03e-01): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=0.01, first variable), R's Wfd scores 4.005844 / 1.376484 / 1.724921 / 3.736952 / 7.200670 for the five boys, with gradients of max-norm 0.19 to 2.09 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 3.993467 / 1.369739 / 1.716418 / 3.722099 / 6.791805 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 2.01 in W.
- `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-regfd_coefs]` (registration_multivariate, rel err 1.17e-01): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=1, lambda=0.01, periodic shift, first variable), R's (Wfd, shift) scores 4.772789 / 4.817748 / 13.194021 / 18.376605 for the four boys, with gradients of max-norm 1.09 to 6.20 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 4.578179 / 3.473667 / 12.702899 / 13.294800 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). Shifts differ by up to 2.32 and W by up to 7.6.
- `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-warpfd_coefs]` (registration_multivariate, rel err 1.79e-01): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=1, lambda=0.01, periodic shift, first variable), R's (Wfd, shift) scores 4.772789 / 4.817748 / 13.194021 / 18.376605 for the four boys, with gradients of max-norm 1.09 to 6.20 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 4.578179 / 3.473667 / 12.702899 / 13.294800 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). Shifts differ by up to 2.32 and W by up to 7.6.
- `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-Wfd_coefs]` (registration_multivariate, rel err 1.22e+00): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=1, lambda=0.01, periodic shift, first variable), R's (Wfd, shift) scores 4.772789 / 4.817748 / 13.194021 / 18.376605 for the four boys, with gradients of max-norm 1.09 to 6.20 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 4.578179 / 3.473667 / 12.702899 / 13.294800 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). Shifts differ by up to 2.32 and W by up to 7.6.
- `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-shift]` (registration_multivariate, rel err 9.17e-01): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=1, lambda=0.01, periodic shift, first variable), R's (Wfd, shift) scores 4.772789 / 4.817748 / 13.194021 / 18.376605 for the four boys, with gradients of max-norm 1.09 to 6.20 there; fdatools' Newton iterate is stationary (gradient < 1e-12) at 4.578179 / 3.473667 / 12.702899 / 13.294800 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). Shifts differ by up to 2.32 and W by up to 7.6.
- `test_register_matches_r[register_fd_handwriting_xy_crit2-regfd_coefs]` (registration_multivariate, rel err 2.68e-03): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=1, first variable), R's Wfd scores 0.495565 / 0.367009 / 1.008572 / 0.924378 / 1.371460 for the five samples, with gradients of max-norm 0.064 to 0.26 there; fdatools' Newton iterate is stationary (gradient < 1e-11) at 0.495236 / 0.366892 / 1.008538 / 0.924099 / 1.371380 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 0.0073 in W.
- `test_register_matches_r[register_fd_handwriting_xy_crit2-warpfd_coefs]` (registration_multivariate, rel err 1.58e-04): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=1, first variable), R's Wfd scores 0.495565 / 0.367009 / 1.008572 / 0.924378 / 1.371460 for the five samples, with gradients of max-norm 0.064 to 0.26 there; fdatools' Newton iterate is stationary (gradient < 1e-11) at 0.495236 / 0.366892 / 1.008538 / 0.924099 / 1.371380 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 0.0073 in W.
- `test_register_matches_r[register_fd_handwriting_xy_crit2-Wfd_coefs]` (registration_multivariate, rel err 9.86e-02): register.fd returns a point that is not a minimum of its own criterion. With R's discretisation (grid mean over 211 points, crit=2, lambda=1, first variable), R's Wfd scores 0.495565 / 0.367009 / 1.008572 / 0.924378 / 1.371460 for the five samples, with gradients of max-norm 0.064 to 0.26 there; fdatools' Newton iterate is stationary (gradient < 1e-11) at 0.495236 / 0.366892 / 1.008538 / 0.924099 / 1.371380 -- lower on every curve (checked by test_fdatools_optimum_is_below_r_point). The optima differ by up to 0.0073 in W.
- `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_crit2-warpfd_coefs]` (registration_multivariate, rel err 7.96e-06): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_crit2-regfd_coefs]` (registration_multivariate, rel err 4.66e-03): R's regfd is not x(h(t)) for the warp R returns: it differs from R's own register.newfd(yfd, Wfd) -- x at R's trapezoid warp, reproduced to 1e-14 -- by 4.66e-3 (gait) / 1.02e-4 (handwriting) normwise relative, and from x(h(t) + shift) at R's trapezoid warp by 3.99e-3 (periodic gait). The univariate cases show the same internal warp defect (test_registration.py).
- `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_periodic_crit1-warpfd_coefs]` (registration_multivariate, rel err 2.30e-06): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_periodic_crit1-regfd_coefs]` (registration_multivariate, rel err 3.99e-03): R's regfd is not x(h(t)) for the warp R returns: it differs from R's own register.newfd(yfd, Wfd) -- x at R's trapezoid warp, reproduced to 1e-14 -- by 4.66e-3 (gait) / 1.02e-4 (handwriting) normwise relative, and from x(h(t) + shift) at R's trapezoid warp by 3.99e-3 (periodic gait). The univariate cases show the same internal warp defect (test_registration.py).
- `test_warps_from_r_latent_match_r[register_fd_handwriting_xy_crit2-warpfd_coefs]` (registration_multivariate, rel err 4.95e-08): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_warps_from_r_latent_match_r[register_fd_handwriting_xy_crit2-regfd_coefs]` (registration_multivariate, rel err 1.02e-04): R's regfd is not x(h(t)) for the warp R returns: it differs from R's own register.newfd(yfd, Wfd) -- x at R's trapezoid warp, reproduced to 1e-14 -- by 4.66e-3 (gait) / 1.02e-4 (handwriting) normwise relative, and from x(h(t) + shift) at R's trapezoid warp by 3.99e-3 (periodic gait). The univariate cases show the same internal warp defect (test_registration.py).
- `test_apply_matches_register_newfd[register_fd_gait_hip_knee_crit2]` (registration_multivariate, rel err 1.47e-05): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_apply_matches_register_newfd[register_fd_gait_hip_knee_periodic_crit1]` (registration_multivariate, rel err 5.02e-06): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_apply_matches_register_newfd[register_fd_handwriting_xy_crit2]` (registration_multivariate, rel err 5.53e-07): R computes the warp h(t) from Wfd with a 1025-point trapezoid rule and linear interpolation; that rule reproduces R's warpfd and register.newfd output to 1e-14 (test_r_warps_use_a_trapezoid_rule). It is off the exact integral by up to 9.4e-5 (gait, of 20) / 5.4e-5 (periodic gait) / 9.7e-5 (handwriting, of 2300), which moves R's warpfd coefficients by up to 3.9e-5 / 3.4e-5 / 6.3e-5 relative and R's register.newfd coefficients by 1.5e-5 / 5.0e-6 / 5.5e-7 normwise, up to 1.5e-2 / 7.4e-3 / 1.8e-5 relative on the smallest coefficients. fdatools integrates exp(W) by Gauss-Legendre to rounding error.
- `test_landmark_latent_matches_r[Wfd_coefs_hip]` (registration_multivariate, rel err 1.93e-05): landmarkreg's Wfd does not meet the landmarks it is fitted to: its warps miss the target landmark by 3.9e-6 / 3.9e-7 / 1.3e-7 / 1.6e-6 / 3.3e-7 / 1.1e-5 for the six boys (smooth.morph objective 1.5e-11 to 1.3e-10), while fdatools' meet it to 1e-10 (objective < 1e-23, checked by test_fdatools_landmark_warps_meet_the_landmarks). The Wfd of the boy whose landmark is 0.01 from the target is 2.8e-3 in size, so its 1e-6 gap is 3.9e-4 relative; the largest gap is 3.4e-6 in W.
- `test_landmark_latent_matches_r[Wfd_coefs_knee]` (registration_multivariate, rel err 1.93e-05): landmarkreg's Wfd does not meet the landmarks it is fitted to: its warps miss the target landmark by 3.9e-6 / 3.9e-7 / 1.3e-7 / 1.6e-6 / 3.3e-7 / 1.1e-5 for the six boys (smooth.morph objective 1.5e-11 to 1.3e-10), while fdatools' meet it to 1e-10 (objective < 1e-23, checked by test_fdatools_landmark_warps_meet_the_landmarks). The Wfd of the boy whose landmark is 0.01 from the target is 2.8e-3 in size, so its 1e-6 gap is 3.9e-4 relative; the largest gap is 3.4e-6 in W.
- `test_golden_field[fregress_scalar_precip_on_temp-OCV]` (regression, rel err 1.90e-08): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_scalar_precip_on_temp-Cmat]` (regression, rel err 1.19e-09): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_scalar_precip_on_temp-Dmat]` (regression, rel err 1.28e-09): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_scalar_precip_on_temp-beta1]` (regression, rel err 1.73e-08): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_cv_scalar_precip_on_temp-SSE.CV]` (regression, rel err 1.84e-08): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_cv_scalar_precip_on_temp-errfd.cv]` (regression, rel err 5.77e-08): R integrates int x_i(t) theta_k(t) dt approximately. The response basis is an orthonormal Fourier basis, so the exact integral is simply the k-th coefficient of x_i; against that R's Dmat is up to 3.4e-8 (relative) out and its Cmat up to 8.2e-8. fdatools' Cmat and Dmat are the exact values.
- `test_golden_field[fregress_functional_temp_on_region-yhatfdobj_coefs]` (regression, rel err 1.13e-06): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-Cmat]` (regression, rel err 9.50e-07): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-Dmat]` (regression, rel err 1.38e-06): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-beta0]` (regression, rel err 1.20e-06): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-beta1]` (regression, rel err 1.16e-06): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-beta2]` (regression, rel err 1.03e-06): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_functional_temp_on_region-beta3]` (regression, rel err 8.31e-07): R's Gram matrix of the Fourier(65) coefficient basis is not the identity: the same approximate integration gives inprod(fbasis5, fbasis5)[4, 4] = 0.9999986 for an orthonormal basis. With it R's Dmat is 1.38e-6 and its betas 1.13e-6 out (relative to the largest entry); entry by entry the near-zero coefficients differ by up to 1.1e-2 relative.
- `test_golden_field[fregress_concurrent_logprecip_on_temp-yhatfdobj_coefs]` (regression, rel err 6.22e-03): R's concurrent-model integrals are 0.1-0.3% off. R's own dense periodic trapezoid rule (20000 points, exact for trigonometric polynomials) gives sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) and Dmat[6] 4172.91107 (3.4e-3 low). fdatools reproduces the dense values (5861.19665, 4187.25630). The betas inherit a 1.2e-2 relative error.
- `test_golden_field[fregress_concurrent_logprecip_on_temp-Cmat]` (regression, rel err 1.07e-03): R's concurrent-model integrals are 0.1-0.3% off. R's own dense periodic trapezoid rule (20000 points, exact for trigonometric polynomials) gives sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) and Dmat[6] 4172.91107 (3.4e-3 low). fdatools reproduces the dense values (5861.19665, 4187.25630). The betas inherit a 1.2e-2 relative error.
- `test_golden_field[fregress_concurrent_logprecip_on_temp-Dmat]` (regression, rel err 3.44e-03): R's concurrent-model integrals are 0.1-0.3% off. R's own dense periodic trapezoid rule (20000 points, exact for trigonometric polynomials) gives sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) and Dmat[6] 4172.91107 (3.4e-3 low). fdatools reproduces the dense values (5861.19665, 4187.25630). The betas inherit a 1.2e-2 relative error.
- `test_golden_field[fregress_concurrent_logprecip_on_temp-beta0]` (regression, rel err 1.19e-02): R's concurrent-model integrals are 0.1-0.3% off. R's own dense periodic trapezoid rule (20000 points, exact for trigonometric polynomials) gives sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) and Dmat[6] 4172.91107 (3.4e-3 low). fdatools reproduces the dense values (5861.19665, 4187.25630). The betas inherit a 1.2e-2 relative error.
- `test_golden_field[fregress_concurrent_logprecip_on_temp-beta1]` (regression, rel err 6.61e-03): R's concurrent-model integrals are 0.1-0.3% off. R's own dense periodic trapezoid rule (20000 points, exact for trigonometric polynomials) gives sum_i int temp_i^2 theta_1^2 = 5861.19665 and int theta_1 sum_i temp_i logprecip_i = 4187.25630; fRegress's Cmat[6, 6] is 5855.72278 (9.3e-4 low) and Dmat[6] 4172.91107 (3.4e-3 low). fdatools reproduces the dense values (5861.19665, 4187.25630). The betas inherit a 1.2e-2 relative error.
- `test_golden_field[fregress_synthetic_2scalar_1functional-df]` (regression, rel err 1.51e-06): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-yhatfdobj]` (regression, rel err 3.06e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-gcv]` (regression, rel err 1.96e-06): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-OCV]` (regression, rel err 3.02e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-Cmat]` (regression, rel err 1.16e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-Dmat]` (regression, rel err 1.01e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-beta0]` (regression, rel err 1.07e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-beta1]` (regression, rel err 1.41e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-beta2]` (regression, rel err 3.86e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional-beta3]` (regression, rel err 7.03e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-yhatfdobj]` (regression, rel err 7.00e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-gcv]` (regression, rel err 1.20e-06): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-OCV]` (regression, rel err 3.63e-04): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-Cmat]` (regression, rel err 1.16e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-Dmat]` (regression, rel err 1.01e-05): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-beta3]` (regression, rel err 2.08e-03): R's int x_i(t) theta_k(t) dt for the B-spline covariate is off by up to 1.7e-3 (relative): fRegress's Dmat[9] is -0.2107866 while R's own exact inprod.bspline and a 200001-point trapezoid rule both give -0.2104203, which is fdatools' value. Every quantity downstream of the design inherits the error (betas 1.7e-4 / 4.3e-2 at lambda 1e-2 / 0, fitted values 1.6e-4).
- `smooth_basis_weather_fourier65_harmonic_lambda0.01-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda100-coefs` (smoothing, rel err 2.01e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-gcv` (smoothing, rel err 9.86e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda100-y2cMap` (smoothing, rel err 2.69e-07): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda100-fitted` (smoothing, rel err 1.59e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-coefs` (smoothing, rel err 3.68e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-y2cMap` (smoothing, rel err 4.93e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda10000-fitted` (smoothing, rel err 2.91e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_weather_fourier65_harmonic_lambda1e+06-penmat` (smoothing, rel err 1.38e-06): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02.
- `smooth_basis_weather_fourier65_harmonic_lambda1e+06-y2cMap` (smoothing, rel err 4.96e-11): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-gcv` (smoothing, rel err 4.86e-09): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda0.01-y2cMap` (smoothing, rel err 1.11e-09): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-coefs` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-df` (smoothing, rel err 1.33e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-gcv` (smoothing, rel err 4.93e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-SSE` (smoothing, rel err 5.09e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-y2cMap` (smoothing, rel err 8.52e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda1-fitted` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-coefs` (smoothing, rel err 4.38e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-df` (smoothing, rel err 1.61e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-gcv` (smoothing, rel err 2.34e-05): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-SSE` (smoothing, rel err 1.08e-05): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-y2cMap` (smoothing, rel err 6.47e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basis_growth_hgtm_bspline6_lfd4_lambda100-fitted` (smoothing, rel err 4.39e-06): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `lambda2gcv_weather_log10lambda_2-gcv` (smoothing, rel err 9.86e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `df2lambda_weather_df10-lambda` (smoothing, rel err 7.70e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so fdatools bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `df2lambda_weather_df20-lambda` (smoothing, rel err 9.83e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so fdatools bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `df2lambda_weather_df30-lambda` (smoothing, rel err 8.56e-04): R's df2lambda stops its search short of solving df(lambda) = df. Feeding R's own answer back through df(lambda) gives 9.998754682 / 20.003244908 / 30.004195052 for the targets 10 / 20 / 30. df(lambda) is smooth and strictly decreasing, so fdatools bisects it to 1e-12 and lands on lambda 3547717.97 / 51683.3869 / 4407.61711, which reproduce the targets to twelve decimals; the resulting lambdas differ from R's by 7.7e-4 / 9.8e-4 / 8.6e-4 (relative).
- `lambda2gcv_weather_grid_mean-mean_gcv` (smoothing, rel err 1.46e-08): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `smooth_basisPar_growth_hgtm-coefs` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-df` (smoothing, rel err 1.33e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-gcv` (smoothing, rel err 4.93e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-SSE` (smoothing, rel err 5.09e-07): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-y2cMap` (smoothing, rel err 8.52e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_basisPar_growth_hgtm-fitted` (smoothing, rel err 5.63e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `data2fd_weather_fourier-coefs` (smoothing, rel err 2.01e-09): R's Fourier penalty under the harmonic accelerator L = w^2 D + D^3 (w = 2*pi/365) is wrong in its 6th significant digit. L annihilates the constant and the fundamental and maps cos(k w t) to k w^3 (k^2 - 1) sin(k w t), so on the normalised Fourier basis the penalty is exactly diagonal with entries k^2 w^6 (k^2 - 1)^2. Against that closed form fdatools is 2.1e-15 (relative) out and R is 1.383e-6: R's worst entry is (64, 64), exact 2.788516685330926e-02 against R's 2.788512828061242e-02. At lambda >= 100 the wrong penalty moves the solution of (Phi'Phi + lambda R) c = Phi'y, so the coefficients (1.0e-6 / 1.4e-6 relative at lambda = 1e2 / 1e4), the hat map, the fitted values and the GCV score inherit the same error.
- `data2fd_growth_hgtf-coefs` (smoothing, rel err 5.80e-08): Not reachable in double precision by either implementation: 31 ages on 35 order-6 B-splines give cond(Phi'Phi + lambda R) = 1.0e7 / 8.3e8 / 7.8e10 at lambda = 0.01 / 1 / 100. Both penalty matrices are correct -- against a 50-digit mpmath integration of the D^4 penalty (exact, since D^4 of an order-6 spline is piecewise linear) fdatools is 2.2e-16 out and R 1.7e-16, i.e. both sit on the rounding floor. But solving the two systems in 60-digit arithmetic gives coefficients that differ by 7.3e-10 / 5.4e-8 / 3.6e-6 -- the same size as the fdatools-vs-R gap being asserted here. The 1e-8 golden tolerance is below the noise floor of the problem, not a sign that either solver is wrong.
- `smooth_monotone_growth_hgtf_3girls-Wfdobj_coefs` (smoothing, rel err 8.92e-03): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_growth_hgtf_3girls-beta_slope` (smoothing, rel err 7.16e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_growth_hgtf_3girls-deriv1` (smoothing, rel err 6.77e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-Wfdobj_coefs` (smoothing, rel err 2.33e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-beta_slope` (smoothing, rel err 2.00e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_monotone_synthetic_sigmoid-deriv1` (smoothing, rel err 1.96e-02): smooth.monotone has an exactly flat direction, so the coefficients are not identified: x(t) = beta0 + beta1 * integral exp(W) is unchanged by W -> W + s together with beta1 -> beta1 * exp(-s), and the roughness penalty is unchanged too because D^m annihilates constants. fdatools and R settle on different members of that one-parameter family: the measured shift is constant across every coefficient of a curve (-0.08723 / -0.06885 / -0.05558 for the three girls, -0.01985 for the sigmoid), and beta1 differs by exactly the matching factor (4.76127 against 4.36367 = 4.76127 * exp(-0.08723)). What the model does identify agrees: the fitted values to 6.6e-7 / 3.0e-6 and beta0 to the same.
- `smooth_pos_synthetic_2curves-Wfdobj_coefs` (smoothing, rel err 6.49e-05): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for fdatools, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum fdatools reports.
- `smooth_pos_synthetic_2curves-fitted` (smoothing, rel err 1.75e-05): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for fdatools, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum fdatools reports.
- `smooth_pos_weather_precip_pr_rupert-Wfdobj_coefs` (smoothing, rel err 5.03e-08): R's smooth.pos returns a point that is not stationary for its own criterion mean((y - exp(Phi c))^2) + lambda c'Rc. The gradient there has norm 2.5e-4 (synthetic) and 1.2e-6 (Prince Rupert) against 2.6e-7 and 9.6e-13 for fdatools, and R's criterion value is the higher of the two (1.77593716848 against 1.77593715997; 3.19972667966087 against 3.19972667965875). One Gauss-Newton step from R's coefficients moves them by 3.43e-6 -- the whole 3.48e-6 gap being asserted -- so R is literally one un-taken iteration short of the minimum fdatools reports.
- `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-Fobs]` (stats, rel err 9.85e-05): R's fRegress assembles its normal equations with inprod(), whose Romberg quadrature is inexact on a Fourier basis: R's inprod(basis, basis) for the 25-function Fourier basis on [0, 365] departs from the identity -- the exact Gram matrix of that orthonormal basis -- by up to 1.383e-6, and the right-hand side Dmat carries a 7.28e-6 relative error. Feeding R's own Cmat and Dmat (read from fRegress()$Cmat / $Dmat) through fdatools' solver and F statistic reproduces R's betas to 1.5e-15 and R's Fobs to 3.0e-15, so the whole gap is R's quadrature; with the exact Gram matrices the maximal F statistic is 0.44257145697174 against R's 0.44261503693042 (9.85e-5 relative). Every refit under a permutation inherits the same error: Fnull is off by 3.9e-5 (median) and 2.5e-4 (worst) relative, and its 5% quantile qval by 1.1e-4. The permutations themselves are R's exact draws (see RRandom), and pval, which only compares Fobs with Fnull, matches.
- `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-Fnull]` (stats, rel err 4.65e-05): R's fRegress assembles its normal equations with inprod(), whose Romberg quadrature is inexact on a Fourier basis: R's inprod(basis, basis) for the 25-function Fourier basis on [0, 365] departs from the identity -- the exact Gram matrix of that orthonormal basis -- by up to 1.383e-6, and the right-hand side Dmat carries a 7.28e-6 relative error. Feeding R's own Cmat and Dmat (read from fRegress()$Cmat / $Dmat) through fdatools' solver and F statistic reproduces R's betas to 1.5e-15 and R's Fobs to 3.0e-15, so the whole gap is R's quadrature; with the exact Gram matrices the maximal F statistic is 0.44257145697174 against R's 0.44261503693042 (9.85e-5 relative). Every refit under a permutation inherits the same error: Fnull is off by 3.9e-5 (median) and 2.5e-4 (worst) relative, and its 5% quantile qval by 1.1e-4. The permutations themselves are R's exact draws (see RRandom), and pval, which only compares Fobs with Fnull, matches.
- `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-qval]` (stats, rel err 1.11e-04): R's fRegress assembles its normal equations with inprod(), whose Romberg quadrature is inexact on a Fourier basis: R's inprod(basis, basis) for the 25-function Fourier basis on [0, 365] departs from the identity -- the exact Gram matrix of that orthonormal basis -- by up to 1.383e-6, and the right-hand side Dmat carries a 7.28e-6 relative error. Feeding R's own Cmat and Dmat (read from fRegress()$Cmat / $Dmat) through fdatools' solver and F statistic reproduces R's betas to 1.5e-15 and R's Fobs to 3.0e-15, so the whole gap is R's quadrature; with the exact Gram matrices the maximal F statistic is 0.44257145697174 against R's 0.44261503693042 (9.85e-5 relative). Every refit under a permutation inherits the same error: Fnull is off by 3.9e-5 (median) and 2.5e-4 (worst) relative, and its 5% quantile qval by 1.1e-4. The permutations themselves are R's exact draws (see RRandom), and pval, which only compares Fobs with Fnull, matches.

## Every check

| Module | Case | Symbol | rtol | Max abs err | Max rel err | Status |
|---|---|---|---:|---:|---:|---|
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-df]` | `stats.confidence_band` | 1.00e-08 | 7.28e-14 | 5.27e-15 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-SSE]` | `stats.confidence_band` | 1.00e-08 | 1.48e-12 | 3.08e-15 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-sigma2]` | `stats.confidence_band` | 1.00e-08 | 1.25e-16 | 3.20e-15 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-stderr]` | `stats.confidence_band` | 1.00e-08 | 7.08e-15 | 1.10e-13 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-stderr_deriv1]` | `stats.confidence_band` | 1.00e-08 | 2.88e-16 | 7.77e-14 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-z]` | `stats.confidence_band` | 1.00e-08 | 0 | 0 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-lower_curve1]` | `stats.confidence_band` | 1.00e-08 | 6.71e-14 | 1.10e-13 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_scalar-upper_curve1]` | `stats.confidence_band` | 1.00e-08 | 9.48e-14 | 1.22e-13 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_pointwise-varvec]` | `stats.confidence_band` | 1.00e-08 | 9.95e-15 | 1.01e-13 | pass |
| bands | `test_golden_field[band_smooth_logprec_bspline53_sigma_pointwise-stderr]` | `stats.confidence_band` | 1.00e-08 | 5.82e-15 | 1.01e-13 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-beta]` | `stats.confidence_band` | 1.00e-08 | 1.85e-15 | 4.51e-16 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-sigma2]` | `stats.confidence_band` | 1.00e-08 | 1.04e-16 | 2.97e-15 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-bvar]` | `stats.confidence_band` | 1.00e-08 | 8.33e-16 | 3.65e-15 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-stderr]` | `stats.confidence_band` | 1.00e-08 | 8.88e-16 | 1.86e-15 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-lower]` | `stats.confidence_band` | 1.00e-08 | 1.82e-15 | 5.75e-16 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_meantemp_latitude-upper]` | `stats.confidence_band` | 1.00e-08 | 3.55e-15 | 7.04e-16 | pass |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-beta1]` | `stats.confidence_band` | 1.00e-08 | 9.48e-07 | 2.15e-04 | xfail |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-df]` | `stats.confidence_band` | 1.00e-08 | 2.49e-07 | 3.11e-08 | xfail |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-sigma2]` | `stats.confidence_band` | 1.00e-08 | 3.43e-07 | 1.80e-05 | xfail |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-bvar]` | `stats.confidence_band` | 1.00e-08 | 1.56e-09 | 9.20e-05 | xfail |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-stderr0]` | `stats.confidence_band` | 1.00e-08 | 2.76e-09 | 3.49e-06 | xfail |
| bands | `test_golden_field[band_fregress_logannualprec_on_temp_bspline7-stderr1]` | `stats.confidence_band` | 1.00e-08 | 2.03e-07 | 4.93e-05 | xfail |
| bands | `test_every_defect_names_a_real_field` | `stats.confidence_band` | — | — | — | pass |
| bands | `test_every_case_has_an_accessor` | `stats.confidence_band` | — | — | — | pass |
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
| datasets | `test_missing_release_dir_skips` | `datasets.load_*` | — | — | — | pass |
| decomposition | `test_every_golden_case_is_replayed` | `FPCA` | — | — | — | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda0-values` | `FPCA` | 1.00e-08 | 5.36e-08 | 3.53e-12 | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda0-harmonics` | `FPCA` | 1.00e-08 | 3.25e-09 | 3.52e-09 | xfail |
| decomposition | `pca_fd_weather_nharm2_harmlambda0-scores` | `FPCA` | 1.00e-08 | 6.60e-01 | 1.91e-03 | xfail |
| decomposition | `pca_fd_weather_nharm2_harmlambda0-varprop` | `FPCA` | 1.00e-08 | 1.53e-12 | 1.73e-12 | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda0-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda10000-values` | `FPCA` | 1.00e-08 | 2.62e-10 | 1.73e-14 | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda10000-harmonics` | `FPCA` | 1.00e-08 | 2.32e-11 | 2.51e-11 | xfail |
| decomposition | `pca_fd_weather_nharm2_harmlambda10000-scores` | `FPCA` | 1.00e-08 | 2.65e+00 | 7.69e-03 | xfail |
| decomposition | `pca_fd_weather_nharm2_harmlambda10000-varprop` | `FPCA` | 1.00e-08 | 7.36e-15 | 8.29e-15 | pass |
| decomposition | `pca_fd_weather_nharm2_harmlambda10000-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda0-values` | `FPCA` | 1.00e-08 | 5.36e-08 | 3.53e-12 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda0-harmonics` | `FPCA` | 1.00e-08 | 5.09e-09 | 5.51e-09 | xfail |
| decomposition | `pca_fd_weather_nharm4_harmlambda0-scores` | `FPCA` | 1.00e-08 | 7.70e-01 | 2.23e-03 | xfail |
| decomposition | `pca_fd_weather_nharm4_harmlambda0-varprop` | `FPCA` | 1.00e-08 | 1.53e-12 | 1.73e-12 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda0-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda10000-values` | `FPCA` | 1.00e-08 | 2.62e-10 | 1.73e-14 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda10000-harmonics` | `FPCA` | 1.00e-08 | 3.47e-11 | 3.76e-11 | xfail |
| decomposition | `pca_fd_weather_nharm4_harmlambda10000-scores` | `FPCA` | 1.00e-08 | 6.01e-01 | 1.74e-03 | xfail |
| decomposition | `pca_fd_weather_nharm4_harmlambda10000-varprop` | `FPCA` | 1.00e-08 | 7.36e-15 | 8.29e-15 | pass |
| decomposition | `pca_fd_weather_nharm4_harmlambda10000-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-values` | `FPCA` | 1.00e-08 | 2.88e+01 | 4.41e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-harmonics` | `FPCA` | 1.00e-08 | 7.77e-04 | 1.11e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-scores` | `FPCA` | 1.00e-08 | 9.18e-01 | 3.91e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-varprop` | `FPCA` | 1.00e-08 | 1.56e-03 | 4.10e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda0-rotmat` | `FPCA` | 1.00e-08 | 8.07e-04 | 1.10e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-values` | `FPCA` | 1.00e-08 | 1.39e+01 | 2.03e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-harmonics` | `FPCA` | 1.00e-08 | 3.35e-03 | 4.85e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-scores` | `FPCA` | 1.00e-08 | 5.15e-01 | 2.24e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-varprop` | `FPCA` | 1.00e-08 | 8.80e-04 | 2.20e-03 | xfail |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-meanfd` | `FPCA` | 1.00e-08 | 2.84e-14 | 1.47e-16 | pass |
| decomposition | `varmx_pca_fd_pca_fd_weather_nharm4_harmlambda10000-rotmat` | `FPCA` | 1.00e-08 | 3.90e-03 | 5.41e-03 | xfail |
| decomposition | `pca_fd_growth_hgtm_nharm3-values` | `FPCA` | 1.00e-08 | 4.62e-10 | 9.22e-13 | pass |
| decomposition | `pca_fd_growth_hgtm_nharm3-harmonics` | `FPCA` | 1.00e-08 | 3.19e-05 | 4.57e-05 | xfail |
| decomposition | `pca_fd_growth_hgtm_nharm3-scores` | `FPCA` | 1.00e-08 | 2.20e-03 | 4.51e-05 | xfail |
| decomposition | `pca_fd_growth_hgtm_nharm3-varprop` | `FPCA` | 1.00e-08 | 1.25e-12 | 1.38e-12 | pass |
| decomposition | `pca_fd_growth_hgtm_nharm3-meanfd` | `FPCA` | 1.00e-08 | 5.68e-14 | 3.16e-16 | pass |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-values` | `FPCA` | 1.00e-08 | 1.98e-02 | 7.21e-05 | xfail |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-harmonics` | `FPCA` | 1.00e-08 | 8.91e-05 | 1.21e-04 | xfail |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-scores` | `FPCA` | 1.00e-08 | 2.37e-03 | 6.90e-05 | xfail |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-varprop` | `FPCA` | 1.00e-08 | 3.03e-05 | 6.13e-05 | xfail |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-meanfd` | `FPCA` | 1.00e-08 | 5.68e-14 | 3.16e-16 | pass |
| decomposition | `varmx_pca_fd_growth_hgtm_nharm3-rotmat` | `FPCA` | 1.00e-08 | 1.10e-04 | 1.23e-04 | xfail |
| decomposition | `pca_fd_synthetic_bspline_n10curves-values` | `FPCA` | 1.00e-08 | 1.04e-04 | 5.30e-05 | xfail |
| decomposition | `pca_fd_synthetic_bspline_n10curves-harmonics` | `FPCA` | 1.00e-08 | 8.51e-04 | 4.92e-04 | xfail |
| decomposition | `pca_fd_synthetic_bspline_n10curves-scores` | `FPCA` | 1.00e-08 | 5.60e-04 | 1.90e-04 | xfail |
| decomposition | `pca_fd_synthetic_bspline_n10curves-varprop` | `FPCA` | 1.00e-08 | 1.50e-05 | 3.17e-05 | xfail |
| decomposition | `pca_fd_synthetic_bspline_n10curves-meanfd` | `FPCA` | 1.00e-08 | 1.11e-16 | 1.42e-16 | pass |
| decomposition | `pca_fd_gait_multivariate_nharm3-values` | `FPCA` | 1.00e-08 | 5.39e-07 | 7.56e-10 | xfail |
| decomposition | `pca_fd_gait_multivariate_nharm3-harmonics` | `FPCA` | 1.00e-08 | 6.05e-08 | 6.87e-08 | xfail |
| decomposition | `pca_fd_gait_multivariate_nharm3-scores` | `FPCA` | 1.00e-08 | 3.72e-02 | 6.51e-04 | xfail |
| decomposition | `pca_fd_gait_multivariate_nharm3-varprop` | `FPCA` | 1.00e-08 | 5.52e-10 | 1.23e-09 | pass |
| decomposition | `pca_fd_gait_multivariate_nharm3-meanfd` | `FPCA` | 1.00e-08 | 5.68e-14 | 4.24e-16 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda100-ccacorr` | `FCCA` | 1.00e-08 | 1.68e-15 | 1.68e-15 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda100-ccawtfd1` | `FCCA` | 1.00e-08 | 2.99e-13 | 4.57e-13 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda100-ccawtfd2` | `FCCA` | 1.00e-08 | 4.12e-13 | 6.32e-13 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda100-ccavar1` | `FCCA` | 1.00e-08 | 3.93e-12 | 3.85e-13 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda100-ccavar2` | `FCCA` | 1.00e-08 | 2.07e-12 | 3.22e-13 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda1e+06-ccacorr` | `FCCA` | 1.00e-08 | 1.11e-15 | 1.17e-15 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda1e+06-ccawtfd1` | `FCCA` | 1.00e-08 | 4.50e-15 | 5.02e-15 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda1e+06-ccawtfd2` | `FCCA` | 1.00e-08 | 1.55e-15 | 1.58e-15 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda1e+06-ccavar1` | `FCCA` | 1.00e-08 | 8.17e-14 | 1.59e-15 | pass |
| decomposition | `cca_fd_weather_temp_vs_logprecip_lambda1e+06-ccavar2` | `FCCA` | 1.00e-08 | 3.91e-14 | 3.44e-15 | pass |
| decomposition | `test_pca_harmonics_follow_r_sign_rule` | `FPCA` | — | — | — | pass |
| density | `intensity_homogeneous_bspline13_L1_lam10-coefs` | `fit_intensity` | 1.00e-05 | 1.35e-06 | 1.16e-06 | pass |
| density | `intensity_homogeneous_bspline13_L1_lam10-log_intensity` | `fit_intensity` | 1.00e-05 | 4.29e-07 | 3.71e-07 | pass |
| density | `intensity_homogeneous_bspline13_L1_lam10-intensity` | `fit_intensity` | 1.00e-05 | 9.94e-07 | 3.12e-07 | pass |
| density | `intensity_homogeneous_bspline13_L1_lam10-f` | `fit_intensity` | 1.00e-05 | 2.34e-06 | 7.22e-08 | pass |
| density | `intensity_sine_bspline23_L2_lam1-coefs` | `fit_intensity` | 1.00e-05 | 7.65e-07 | 3.29e-07 | pass |
| density | `intensity_sine_bspline23_L2_lam1-log_intensity` | `fit_intensity` | 1.00e-05 | 5.92e-07 | 2.58e-07 | pass |
| density | `intensity_sine_bspline23_L2_lam1-intensity` | `fit_intensity` | 1.00e-05 | 5.88e-06 | 5.92e-07 | pass |
| density | `intensity_sine_bspline23_L2_lam1-f` | `fit_intensity` | 1.00e-05 | 1.54e-04 | 1.81e-06 | pass |
| density | `intensity_sine_bspline11_order3_L1_lam0-coefs` | `fit_intensity` | 1.00e-05 | 7.30e-06 | 2.17e-06 | xfail |
| density | `intensity_sine_bspline11_order3_L1_lam0-log_intensity` | `fit_intensity` | 1.00e-05 | 4.17e-06 | 1.62e-06 | pass |
| density | `intensity_sine_bspline11_order3_L1_lam0-intensity` | `fit_intensity` | 1.00e-05 | 3.65e-05 | 2.79e-06 | pass |
| density | `intensity_sine_bspline11_order3_L1_lam0-f` | `fit_intensity` | 1.00e-05 | 1.64e-04 | 1.87e-06 | pass |
| density | `intensity_sine_fourier7_L2_lam0.1-coefs` | `fit_intensity` | 1.00e-05 | 4.49e-07 | 6.49e-08 | pass |
| density | `intensity_sine_fourier7_L2_lam0.1-log_intensity` | `fit_intensity` | 1.00e-05 | 2.15e-07 | 9.03e-08 | pass |
| density | `intensity_sine_fourier7_L2_lam0.1-intensity` | `fit_intensity` | 1.00e-05 | 1.39e-06 | 1.29e-07 | pass |
| density | `intensity_sine_fourier7_L2_lam0.1-f` | `fit_intensity` | 1.00e-05 | 4.98e-07 | 5.81e-09 | pass |
| density | `density_regina_bspline13_L2_lam0.1-coefs` | `fit_density` | 1.00e-05 | 1.05e-05 | 2.64e-06 | xfail |
| density | `density_regina_bspline13_L2_lam0.1-log_intensity` | `fit_density` | 1.00e-05 | 4.27e-06 | 1.24e-06 | xfail |
| density | `density_regina_bspline13_L2_lam0.1-intensity` | `fit_density` | 1.00e-05 | 1.33e-04 | 4.27e-06 | pass |
| density | `density_regina_bspline13_L2_lam0.1-f` | `fit_density` | 1.00e-05 | 3.11e-05 | 1.01e-07 | pass |
| density | `density_regina_bspline13_L2_lam0.1-density` | `fit_density` | 1.00e-05 | 6.46e-07 | 4.27e-06 | pass |
| density | `density_regina_bspline13_L2_lam0.1-log_density` | `fit_density` | 1.00e-05 | 4.27e-06 | 5.12e-07 | pass |
| density | `density_regina_bspline13_L2_lam10-coefs` | `fit_density` | 1.00e-05 | 3.02e-06 | 8.68e-07 | pass |
| density | `density_regina_bspline13_L2_lam10-log_intensity` | `fit_density` | 1.00e-05 | 1.43e-06 | 4.11e-07 | pass |
| density | `density_regina_bspline13_L2_lam10-intensity` | `fit_density` | 1.00e-05 | 4.64e-05 | 1.43e-06 | pass |
| density | `density_regina_bspline13_L2_lam10-f` | `fit_density` | 1.00e-05 | 1.16e-05 | 3.80e-08 | pass |
| density | `density_regina_bspline13_L2_lam10-density` | `fit_density` | 1.00e-05 | 2.25e-07 | 1.43e-06 | pass |
| density | `density_regina_bspline13_L2_lam10-log_density` | `fit_density` | 1.00e-05 | 1.43e-06 | 2.15e-07 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-coefs` | `fit_density` | 1.00e-05 | 8.40e-07 | 1.87e-07 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-log_intensity` | `fit_density` | 1.00e-05 | 8.40e-07 | 1.92e-07 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-intensity` | `fit_density` | 1.00e-05 | 2.87e-06 | 3.60e-08 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-f` | `fit_density` | 1.00e-05 | 2.06e-07 | 3.52e-10 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-density` | `fit_density` | 1.00e-05 | 1.44e-08 | 3.60e-08 | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01-log_density` | `fit_density` | 1.00e-05 | 8.40e-07 | 8.69e-08 | pass |
| density | `density_gamma_bspline9_L2_lam1-coefs` | `fit_density` | 1.00e-05 | 1.17e-06 | 2.37e-07 | pass |
| density | `density_gamma_bspline9_L2_lam1-log_intensity` | `fit_density` | 1.00e-05 | 6.83e-07 | 1.38e-07 | xfail |
| density | `density_gamma_bspline9_L2_lam1-intensity` | `fit_density` | 1.00e-05 | 1.42e-05 | 3.09e-07 | pass |
| density | `density_gamma_bspline9_L2_lam1-f` | `fit_density` | 1.00e-05 | 2.48e-06 | 7.56e-09 | pass |
| density | `density_gamma_bspline9_L2_lam1-density` | `fit_density` | 1.00e-05 | 9.47e-08 | 3.09e-07 | pass |
| density | `density_gamma_bspline9_L2_lam1-log_density` | `fit_density` | 1.00e-05 | 6.83e-07 | 6.87e-08 | pass |
| density | `density_truncnormal_monomial3_L1_lam0-coefs` | `fit_density` | 1.00e-05 | 1.05e-07 | 2.91e-08 | pass |
| density | `density_truncnormal_monomial3_L1_lam0-log_intensity` | `fit_density` | 1.00e-05 | 1.19e-06 | 3.07e-07 | xfail |
| density | `density_truncnormal_monomial3_L1_lam0-intensity` | `fit_density` | 1.00e-05 | 6.15e-06 | 1.29e-07 | pass |
| density | `density_truncnormal_monomial3_L1_lam0-f` | `fit_density` | 1.00e-05 | 1.18e-05 | 4.11e-08 | pass |
| density | `density_truncnormal_monomial3_L1_lam0-density` | `fit_density` | 1.00e-05 | 5.13e-08 | 1.29e-07 | pass |
| density | `density_truncnormal_monomial3_L1_lam0-log_density` | `fit_density` | 1.00e-05 | 1.19e-06 | 1.55e-07 | pass |
| density | `intensity_homogeneous_bspline13_L1_lam10` | `fit_intensity` | — | — | — | pass |
| density | `intensity_sine_bspline23_L2_lam1` | `fit_intensity` | — | — | — | pass |
| density | `intensity_sine_bspline11_order3_L1_lam0` | `fit_intensity` | — | — | — | pass |
| density | `intensity_sine_fourier7_L2_lam0.1` | `fit_intensity` | — | — | — | pass |
| density | `density_regina_bspline13_L2_lam0.1` | `fit_density` | — | — | — | pass |
| density | `density_regina_bspline13_L2_lam10` | `fit_density` | — | — | — | pass |
| density | `density_normal_bspline12_order5_L3_lam0.01` | `fit_density` | — | — | — | pass |
| density | `density_gamma_bspline9_L2_lam1` | `fit_density` | — | — | — | pass |
| density | `density_truncnormal_monomial3_L1_lam0` | `fit_density` | — | — | — | pass |
| density | `test_regina_sample_is_the_dataset` | `fit_density` | — | — | — | pass |
| dynamics | `test_every_golden_case_is_covered` | `PDA` | — | — | — | pass |
| dynamics | `test_input_curves_reproduce_the_golden_design` | `PDA` | 1.00e-12 | 7.22e-15 | 7.22e-15 | pass |
| dynamics | `test_weight_coefficients[pda_fd_order1_analytic_exp_decay-0-bwt_coefs-None]` | `PDA` | 1.00e-08 | 0 | 0 | pass |
| dynamics | `test_weight_coefficients[pda_fd_order2_lip_constant-0-bwt0_coefs-bwt0_vals]` | `PDA` | 1.00e-08 | 1.41e-12 | 2.59e-14 | pass |
| dynamics | `test_weight_coefficients[pda_fd_order2_lip_constant-1-bwt1_coefs-bwt1_vals]` | `PDA` | 1.00e-08 | 8.44e-14 | 3.97e-14 | pass |
| dynamics | `test_weight_coefficients[pda_fd_order2_lip_bspline11-0-bwt0_coefs-bwt0_vals]` | `PDA` | 1.00e-08 | 2.33e-12 | 1.10e-14 | pass |
| dynamics | `test_weight_coefficients[pda_fd_order2_lip_bspline11-1-bwt1_coefs-bwt1_vals]` | `PDA` | 1.00e-08 | 1.48e-12 | 1.03e-13 | pass |
| dynamics | `test_weight_values[pda_fd_order1_analytic_exp_decay-0-bwt_coefs-None]` | `PDA` | 1.00e-08 | 0 | 0 | pass |
| dynamics | `test_weight_values[pda_fd_order2_lip_constant-0-bwt0_coefs-bwt0_vals]` | `PDA` | 1.00e-08 | 1.41e-12 | 2.59e-14 | pass |
| dynamics | `test_weight_values[pda_fd_order2_lip_constant-1-bwt1_coefs-bwt1_vals]` | `PDA` | 1.00e-08 | 8.44e-14 | 3.97e-14 | pass |
| dynamics | `test_weight_values[pda_fd_order2_lip_bspline11-0-bwt0_coefs-bwt0_vals]` | `PDA` | 1.00e-08 | 1.88e-12 | 9.53e-15 | pass |
| dynamics | `test_weight_values[pda_fd_order2_lip_bspline11-1-bwt1_coefs-bwt1_vals]` | `PDA` | 1.00e-08 | 1.48e-12 | 1.03e-13 | pass |
| dynamics | `test_residual_coefficients[pda_fd_order1_analytic_exp_decay]` | `PDA` | 1.00e-08 | 8.92e-15 | 1.04e-09 | pass |
| dynamics | `test_residual_coefficients[pda_fd_order2_lip_bspline11]` | `PDA` | 1.00e-08 | 2.63e-10 | 4.05e-13 | pass |
| dynamics | `test_residual_coefficients[pda_fd_order2_lip_constant]` | `PDA` | 1.00e-08 | 5.40e-10 | 1.91e-13 | pass |
| dynamics | `test_transform_reproduces_the_residuals[pda_fd_order1_analytic_exp_decay]` | `PDA` | 1.00e-08 | 8.92e-15 | 1.04e-09 | pass |
| dynamics | `test_transform_reproduces_the_residuals[pda_fd_order2_lip_bspline11]` | `PDA` | 1.00e-08 | 2.63e-10 | 4.05e-13 | pass |
| dynamics | `test_transform_reproduces_the_residuals[pda_fd_order2_lip_constant]` | `PDA` | 1.00e-08 | 5.40e-10 | 1.91e-13 | pass |
| dynamics | `test_solution_matches_the_r_estimate` | `PDA` | 1.00e-08 | 6.20e-11 | 6.20e-11 | pass |
| dynamics | `test_solution_tracks_the_analytic_ground_truth` | `PDA` | 1.00e-15 | 0 | 0 | pass |
| io | `test_fd_coefs_and_eval_match_r` | `read_rds` | 1.00e-08 | 2.22e-16 | 9.98e-17 | pass |
| io | `test_bifd_coefs_and_eval_match_r` | `read_rds` | 1.00e-08 | 5.55e-17 | 2.43e-17 | pass |
| linmod | `test_golden_field[linmod_weather_logprecip_on_temp-beta0estfd_coefs]` | `linmod` | 1.00e-08 | 8.71e-06 | 4.68e-07 | xfail |
| linmod | `test_golden_field[linmod_weather_logprecip_on_temp-beta1estbifd_coefs]` | `linmod` | 1.00e-08 | 1.71e-06 | 1.88e-06 | xfail |
| linmod | `test_golden_field[linmod_weather_logprecip_on_temp-yhatfdobj_coefs]` | `linmod` | 1.00e-08 | 4.31e-06 | 3.72e-07 | xfail |
| linmod | `test_golden_field[linmod_synthetic_monomial-beta0estfd_coefs]` | `linmod` | 1.00e-08 | 3.21e-14 | 5.33e-14 | pass |
| linmod | `test_golden_field[linmod_synthetic_monomial-beta1estbifd_coefs]` | `linmod` | 1.00e-08 | 2.81e-13 | 1.27e-13 | pass |
| linmod | `test_golden_field[linmod_synthetic_monomial-yhatfdobj_coefs]` | `linmod` | 1.00e-08 | 9.79e-14 | 8.49e-14 | pass |
| linmod | `test_golden_field[linmod_synthetic_monomial_lambda0-beta0estfd_coefs]` | `linmod` | 1.00e-08 | 2.16e-13 | 5.32e-13 | pass |
| linmod | `test_golden_field[linmod_synthetic_monomial_lambda0-beta1estbifd_coefs]` | `linmod` | 1.00e-08 | 1.11e-11 | 1.09e-12 | pass |
| linmod | `test_golden_field[linmod_synthetic_monomial_lambda0-yhatfdobj_coefs]` | `linmod` | 1.00e-08 | 1.83e-12 | 1.13e-12 | pass |
| linmod | `test_golden_field[linmod_synthetic_bspline-beta0estfd_coefs]` | `linmod` | 1.00e-08 | 1.50e-05 | 5.91e-05 | xfail |
| linmod | `test_golden_field[linmod_synthetic_bspline-beta1estbifd_coefs]` | `linmod` | 1.00e-08 | 3.32e-04 | 9.55e-05 | xfail |
| linmod | `test_golden_field[linmod_synthetic_bspline-yhatfdobj_coefs]` | `linmod` | 1.00e-08 | 1.31e-03 | 1.38e-03 | xfail |
| linmod | `test_every_defect_names_a_real_field` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_synthetic_bspline-beta0estfd_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_synthetic_bspline-beta1estbifd_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_synthetic_bspline-yhatfdobj_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_weather_logprecip_on_temp-beta0estfd_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_weather_logprecip_on_temp-beta1estbifd_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_defect_is_as_measured[linmod_weather_logprecip_on_temp-yhatfdobj_coefs]` | `linmod` | — | — | — | pass |
| linmod | `test_normal_equations_with_r_integrals[linmod_synthetic_bspline]` | `linmod` | 1.00e-08 | 2.66e-15 | 7.67e-16 | pass |
| linmod | `test_normal_equations_with_r_integrals[linmod_synthetic_monomial]` | `linmod` | 1.00e-08 | 4.73e-13 | 2.14e-13 | pass |
| linmod | `test_normal_equations_with_r_integrals[linmod_synthetic_monomial_lambda0]` | `linmod` | 1.00e-08 | 4.49e-11 | 4.42e-12 | pass |
| linmod | `test_r_integrals_match_fdatools_where_r_is_exact[linmod_synthetic_monomial]` | `linmod` | 1.00e-12 | 1.78e-15 | 2.91e-16 | pass |
| linmod | `test_r_integrals_match_fdatools_where_r_is_exact[linmod_synthetic_monomial_lambda0]` | `linmod` | 1.00e-12 | 1.78e-15 | 2.78e-16 | pass |
| linmod | `test_r_fitted_curves_are_a_grid_least_squares_fit` | `linmod` | 1.00e-08 | 4.88e-15 | 5.13e-15 | pass |
| monotone_derivs | `test_every_golden_case_is_covered` | `SmoothResult` | — | — | — | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[0-monfd_growth_girl1]` | `SmoothResult` | 1.00e-08 | 4.52e-05 | 2.44e-06 | xfail |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[0-monfd_posfd_fixed_w_two_curves]` | `SmoothResult` | 1.00e-08 | 5.62e-06 | 1.81e-06 | xfail |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[1-monfd_growth_girl1]` | `SmoothResult` | 1.00e-08 | 4.44e-16 | 1.41e-16 | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[1-monfd_posfd_fixed_w_two_curves]` | `SmoothResult` | 1.00e-08 | 4.44e-16 | 1.93e-16 | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[2-monfd_growth_girl1]` | `SmoothResult` | 1.00e-08 | 2.66e-15 | 1.64e-15 | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[2-monfd_posfd_fixed_w_two_curves]` | `SmoothResult` | 1.00e-08 | 1.78e-15 | 2.06e-16 | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[3-monfd_growth_girl1]` | `SmoothResult` | 1.00e-08 | 6.42e-14 | 4.09e-14 | pass |
| monotone_derivs | `test_monotone_derivatives_match_eval_monfd[3-monfd_posfd_fixed_w_two_curves]` | `SmoothResult` | 1.00e-08 | 2.40e-14 | 1.61e-16 | pass |
| monotone_derivs | `test_positive_derivatives_match_eval_posfd[0]` | `SmoothResult` | 1.00e-08 | 4.44e-16 | 1.93e-16 | pass |
| monotone_derivs | `test_positive_derivatives_match_eval_posfd[1]` | `SmoothResult` | 1.00e-08 | 1.78e-15 | 2.06e-16 | pass |
| monotone_derivs | `test_positive_derivatives_match_eval_posfd[2]` | `SmoothResult` | 1.00e-08 | 8.14e+01 | 8.66e-01 | xfail |
| monotone_derivs | `test_scaled_monotone_derivatives_match_predict[0]` | `SmoothResult` | 1.00e-08 | 2.00e-04 | 1.26e-06 | xfail |
| monotone_derivs | `test_scaled_monotone_derivatives_match_predict[1]` | `SmoothResult` | 1.00e-08 | 1.78e-15 | 1.27e-16 | pass |
| monotone_derivs | `test_scaled_monotone_derivatives_match_predict[2]` | `SmoothResult` | 1.00e-08 | 1.15e-14 | 1.61e-15 | pass |
| monotone_derivs | `test_scaled_monotone_derivatives_match_predict[3]` | `SmoothResult` | 1.00e-08 | 2.83e-13 | 4.08e-14 | pass |
| monotone_derivs | `test_monotone_values_match_adaptive_quadrature[monfd_growth_girl1]` | `SmoothResult` | 1.00e-12 | 6.77e-12 | 3.66e-13 | pass |
| monotone_derivs | `test_monotone_values_match_adaptive_quadrature[monfd_posfd_fixed_w_two_curves]` | `SmoothResult` | 1.00e-12 | 4.07e-13 | 1.31e-13 | pass |
| monotone_derivs | `test_r_second_derivative_of_a_positive_function_drops_the_square_term` | `SmoothResult` | 1.00e-14 | 1.42e-14 | 1.51e-16 | pass |
| pace | `smooth_sparse_mean_growth_bspline_lambda0` | `sparse_mean` | 1.00e-08 | 9.81e-13 | 5.69e-15 | pass |
| pace | `smooth_sparse_mean_growth_bspline_lambda1` | `sparse_mean` | 1.00e-08 | 1.17e-12 | 6.77e-15 | pass |
| pace | `smooth_sparse_mean_weather_fourier5` | `sparse_mean` | 1.00e-08 | 1.31e-13 | 6.73e-16 | pass |
| pace | `test_sparse_mean_accepts_r_list_form` | `PACE` | 1.00e-08 | 9.81e-13 | 5.69e-15 | pass |
| pace | `pace_growth_bspline6_covlambda0-mean_coefs` | `PACE` | 1.00e-08 | 9.81e-13 | 5.69e-15 | pass |
| pace | `pace_growth_bspline6_covlambda0-cov_coefs` | `PACE` | 1.00e-08 | 6.62e-12 | 6.09e-14 | pass |
| pace | `pace_growth_bspline6_covlambda0-values` | `PACE` | 1.00e-08 | 3.78e-03 | 7.21e-06 | xfail |
| pace | `pace_growth_bspline6_covlambda0-varprop` | `PACE` | 1.00e-08 | 5.57e-06 | 6.78e-06 | xfail |
| pace | `pace_growth_bspline6_covlambda0-harmonics` | `PACE` | 1.00e-08 | 8.28e-05 | 9.03e-05 | xfail |
| pace | `pace_growth_bspline6_covlambda10-mean_coefs` | `PACE` | 1.00e-08 | 9.81e-13 | 5.69e-15 | pass |
| pace | `pace_growth_bspline6_covlambda10-cov_coefs` | `PACE` | 1.00e-08 | 1.10e-03 | 1.16e-05 | xfail |
| pace | `pace_growth_bspline6_covlambda10-values` | `PACE` | 1.00e-08 | 2.24e-03 | 4.48e-06 | xfail |
| pace | `pace_growth_bspline6_covlambda10-varprop` | `PACE` | 1.00e-08 | 3.33e-06 | 4.05e-06 | xfail |
| pace | `pace_growth_bspline6_covlambda10-harmonics` | `PACE` | 1.00e-08 | 5.42e-05 | 5.86e-05 | xfail |
| pace | `pace_growth_monomial3_penalised-mean_coefs` | `PACE` | 1.00e-08 | 1.79e-12 | 2.54e-14 | pass |
| pace | `pace_growth_monomial3_penalised-cov_coefs` | `PACE` | 1.00e-08 | 3.35e-10 | 4.01e-11 | pass |
| pace | `pace_growth_monomial3_penalised-values` | `PACE` | 1.00e-08 | 1.66e-10 | 3.50e-13 | pass |
| pace | `pace_growth_monomial3_penalised-varprop` | `PACE` | 1.00e-08 | 2.38e-13 | 2.86e-13 | pass |
| pace | `pace_growth_monomial3_penalised-harmonics` | `PACE` | 1.00e-08 | 4.15e-13 | 1.51e-12 | pass |
| pace | `pace_weather_fourier5_bspline7-mean_coefs` | `PACE` | 1.00e-08 | 1.31e-13 | 6.73e-16 | pass |
| pace | `pace_weather_fourier5_bspline7-cov_coefs` | `PACE` | 1.00e-08 | 9.80e-12 | 4.01e-14 | pass |
| pace | `pace_weather_fourier5_bspline7-values` | `PACE` | 1.00e-08 | 1.21e-01 | 8.95e-06 | xfail |
| pace | `pace_weather_fourier5_bspline7-varprop` | `PACE` | 1.00e-08 | 5.55e-06 | 7.03e-06 | xfail |
| pace | `pace_weather_fourier5_bspline7-harmonics` | `PACE` | 1.00e-08 | 2.38e-05 | 1.59e-04 | xfail |
| pace | `pace_weather_fourier5_bspline7-scores` | `PACE` | 1.00e-08 | 1.28e+02 | 4.52e+00 | xfail |
| pace | `pace_growth_bspline6_covlambda0` | `PACE` | 1.00e-08 | 6.62e-12 | 6.09e-14 | pass |
| pace | `pace_growth_bspline6_covlambda10` | `PACE` | 1.00e-08 | 1.78e-12 | 1.87e-14 | pass |
| pace | `pace_growth_monomial3_penalised` | `PACE` | 1.00e-08 | 3.35e-10 | 4.01e-11 | pass |
| pace | `pace_weather_fourier5_bspline7` | `PACE` | 1.00e-08 | 2.00e-11 | 4.01e-14 | pass |
| pace | `test_r_scores_use_one_point` | `PACE` | — | — | — | pass |
| pda_forcing | `test_weight_coefficients[pda_forced_order1_constant_exp]` | `PDA` | 1.00e-08 | 1.24e-14 | 3.11e-15 | pass |
| pda_forcing | `test_weight_coefficients[pda_forced_multicurve_bspline_weights]` | `PDA` | 1.00e-08 | 3.55e-15 | 1.57e-15 | pass |
| pda_forcing | `test_weight_coefficients[pda_forced_refinery_constant]` | `PDA` | 1.00e-08 | 6.94e-18 | 3.09e-16 | pass |
| pda_forcing | `test_weight_coefficients[pda_forced_refinery_two_forcings]` | `PDA` | 1.00e-08 | 1.75e-14 | 9.75e-13 | pass |
| pda_forcing | `test_weight_coefficients[pda_forced_order2_constant]` | `PDA` | 1.00e-08 | 3.55e-14 | 5.07e-13 | pass |
| pda_forcing | `test_forcing_weight_coefficients[pda_forced_order1_constant_exp]` | `PDA` | 1.00e-08 | 7.55e-15 | 3.77e-15 | pass |
| pda_forcing | `test_forcing_weight_coefficients[pda_forced_multicurve_bspline_weights]` | `PDA` | 1.00e-08 | 2.22e-15 | 1.79e-15 | pass |
| pda_forcing | `test_forcing_weight_coefficients[pda_forced_refinery_constant]` | `PDA` | 1.00e-08 | 5.55e-17 | 2.75e-16 | pass |
| pda_forcing | `test_forcing_weight_coefficients[pda_forced_refinery_two_forcings]` | `PDA` | 1.00e-08 | 1.45e-13 | 4.28e-13 | pass |
| pda_forcing | `test_forcing_weight_coefficients[pda_forced_order2_constant]` | `PDA` | 1.00e-08 | 3.55e-15 | 1.78e-16 | pass |
| pda_forcing | `test_residual_coefficients[pda_forced_order1_constant_exp]` | `PDA` | 1.00e-08 | 1.12e-14 | 2.61e-09 | pass |
| pda_forcing | `test_residual_coefficients[pda_forced_multicurve_bspline_weights]` | `PDA` | 1.00e-08 | 1.28e-14 | 8.57e-15 | pass |
| pda_forcing | `test_residual_coefficients[pda_forced_refinery_constant]` | `PDA` | 1.00e-08 | 1.25e-16 | 1.09e-15 | pass |
| pda_forcing | `test_residual_coefficients[pda_forced_refinery_two_forcings]` | `PDA` | 1.00e-08 | 2.81e-15 | 2.82e-14 | pass |
| pda_forcing | `test_residual_coefficients[pda_forced_order2_constant]` | `PDA` | 1.00e-08 | 1.85e-12 | 8.97e-13 | pass |
| pda_forcing | `test_transform_reproduces_the_residuals[pda_forced_order1_constant_exp]` | `PDA` | 1.00e-08 | 1.12e-14 | 2.61e-09 | pass |
| pda_forcing | `test_transform_reproduces_the_residuals[pda_forced_multicurve_bspline_weights]` | `PDA` | 1.00e-08 | 1.28e-14 | 8.57e-15 | pass |
| pda_forcing | `test_transform_reproduces_the_residuals[pda_forced_refinery_constant]` | `PDA` | 1.00e-08 | 1.25e-16 | 1.09e-15 | pass |
| pda_forcing | `test_transform_reproduces_the_residuals[pda_forced_refinery_two_forcings]` | `PDA` | 1.00e-08 | 2.81e-15 | 2.82e-14 | pass |
| pda_forcing | `test_transform_reproduces_the_residuals[pda_forced_order2_constant]` | `PDA` | 1.00e-08 | 1.85e-12 | 8.97e-13 | pass |
| pda_forcing | `test_weight_values_on_a_grid` | `PDA` | 1.00e-08 | 3.55e-15 | 1.79e-15 | pass |
| pda_forcing | `test_system_weights[pda_forced_system_order1]` | `PDA` | 1.00e-08 | 6.79e-13 | 1.32e-06 | pass |
| pda_forcing | `test_system_weights[eigen_pda_handwriting_system]` | `PDA` | 1.00e-08 | 1.67e-10 | 2.18e-12 | pass |
| pda_forcing | `test_system_forcing_weights` | `PDA` | 1.00e-08 | 3.58e-13 | 1.21e-13 | pass |
| pda_forcing | `test_system_residuals[pda_forced_system_order1]` | `PDA` | 1.00e-08 | 1.11e-13 | 1.71e-08 | pass |
| pda_forcing | `test_system_residuals[eigen_pda_handwriting_system]` | `PDA` | 1.00e-08 | 3.22e-13 | 3.57e-13 | pass |
| pda_forcing | `test_stability_eigenvalues[pda_forced_refinery_constant]` | `PDAStability` | 1.00e-08 | 6.94e-18 | 3.09e-16 | pass |
| pda_forcing | `test_stability_eigenvalues[pda_forced_refinery_two_forcings]` | `PDAStability` | 1.00e-08 | 1.75e-14 | 1.23e-12 | pass |
| pda_forcing | `test_stability_eigenvalues[pda_forced_order2_constant]` | `PDAStability` | 1.00e-08 | 2.66e-15 | 5.07e-13 | pass |
| pda_forcing | `test_stability_eigenvalues[pda_forced_system_order1]` | `PDAStability` | 1.00e-08 | 1.26e-12 | 3.16e-13 | pass |
| pda_forcing | `test_stability_eigenvalues[eigen_pda_lip_bspline11]` | `PDAStability` | 1.00e-08 | 7.41e-13 | 1.03e-13 | pass |
| pda_forcing | `test_stability_eigenvalues[eigen_pda_handwriting_system]` | `PDAStability` | 1.00e-08 | 9.35e-12 | 5.06e-13 | pass |
| pda_forcing | `test_unforced_limits_are_zero[eigen_pda_lip_bspline11]` | `PDA` | — | — | — | pass |
| pda_forcing | `test_unforced_limits_are_zero[eigen_pda_handwriting_system]` | `PDA` | — | — | — | pass |
| pda_forcing | `test_first_order_limits_are_r_limvals_with_the_sign_flipped[pda_forced_refinery_constant]` | `PDA` | 1.00e-08 | 1.78e-15 | 3.77e-16 | pass |
| pda_forcing | `test_first_order_limits_are_r_limvals_with_the_sign_flipped[pda_forced_refinery_two_forcings]` | `PDA` | 1.00e-08 | 1.24e-09 | 1.53e-11 | pass |
| pda_forcing | `test_limits_match_r_limvals[pda_forced_refinery_constant]` | `PDAStability` | 1.00e-08 | 9.42e+00 | 2.00e+00 | xfail |
| pda_forcing | `test_limits_match_r_limvals[pda_forced_refinery_two_forcings]` | `PDAStability` | 1.00e-08 | 1.62e+02 | 2.00e+00 | xfail |
| pda_forcing | `test_limits_match_r_limvals[pda_forced_order2_constant]` | `PDAStability` | 1.00e-08 | 5.01e-01 | 1.00e+00 | xfail |
| pda_forcing | `test_limits_match_r_limvals[pda_forced_system_order1]` | `PDAStability` | 1.00e-08 | 1.00e+00 | 2.00e+00 | xfail |
| pda_forcing | `test_every_golden_case_is_covered` | `PDA` | — | — | — | pass |
| pda_forcing | `test_lip_weights_and_residuals` | `PDA` | 1.00e-08 | 2.63e-10 | 4.05e-13 | pass |
| profiling | `test_data_quadrature_is_the_default_rule` | `profile_ode` | 1.00e-12 | 0 | 0 | pass |
| profiling | `quadset_n5_3breaks` | `profile_ode` | 1.00e-08 | 0 | 0 | pass |
| profiling | `quadset_n7_3breaks` | `profile_ode` | 1.00e-08 | 1.11e-16 | 1.78e-16 | pass |
| profiling | `quadset_n5_49breaks` | `profile_ode` | 1.00e-08 | 0 | 0 | pass |
| profiling | `cstr2in_all_cool_step` | `profile_ode` | — | — | — | pass |
| profiling | `cstr2in_all_hot_step` | `profile_ode` | — | — | — | pass |
| profiling | `cstr2in_Tc_hot_step` | `profile_ode` | — | — | — | pass |
| profiling | `cstr2in_Tc_cool_step` | `profile_ode` | — | — | — | pass |
| profiling | `cstr2_default` | `profile_ode` | 1.00e-08 | 4.55e-13 | 2.06e-15 | pass |
| profiling | `cstr2_hot_other_theta` | `profile_ode` | 1.00e-08 | 4.55e-13 | 1.10e-15 | pass |
| profiling | `cstr2_tc_step_constants` | `profile_ode` | 1.00e-08 | 1.14e-13 | 5.08e-16 | pass |
| profiling | `cstr_fitls_fit11_truth` | `profile_ode` | 1.00e-08 | 7.11e-13 | 1.16e-13 | pass |
| profiling | `cstr_fitls_fit01_other` | `profile_ode` | 1.00e-08 | 1.67e-13 | 7.19e-14 | pass |
| profiling | `cstr_fitls_fit10_other` | `profile_ode` | 1.00e-08 | 6.56e-14 | 3.41e-14 | pass |
| profiling | `test_cstr_fn_fails_in_r` | `profile_ode` | — | — | — | pass |
| profiling | `cstr_inner_fit11_truth` | `profile_ode` | 1.00e-08 | 5.25e-11 | 4.06e-14 | pass |
| profiling | `cstr_inner_fit01_other` | `profile_ode` | 1.00e-08 | 1.77e-11 | 2.55e-14 | pass |
| profiling | `cstr_inner_fit11_small_lambda` | `profile_ode` | 1.00e-08 | 5.63e-12 | 2.15e-14 | pass |
| profiling | `cstr_profile_fit11_kref_EoverR` | `profile_ode` | 1.00e-05 | 2.04e-06 | 1.39e-06 | pass |
| profiling | `cstr_profile_fit01_kref_EoverR` | `profile_ode` | 1.00e-05 | 1.47e-07 | 6.25e-08 | pass |
| profiling | `cstr_profile_fit11_kref_EoverR_a` | `profile_ode` | 1.00e-05 | 6.74e-07 | 2.31e-07 | pass |
| registration | `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-regfd_coefs]` | `landmark_register` | 1.00e-05 | 3.28e+03 | 9.77e-01 | xfail |
| registration | `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-warpfd_coefs]` | `landmark_register` | 1.00e-05 | 7.61e-05 | 4.23e-06 | xfail |
| registration | `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-warpinvfd_coefs]` | `landmark_register` | 1.00e-05 | 1.78e+00 | 9.89e-02 | xfail |
| registration | `test_registration_matches_r[landmarkreg_growth_hgtf_pubertal_spurt-Wfd_coefs]` | `landmark_register` | 1.00e-05 | 4.92e-07 | 5.89e-07 | pass |
| registration | `test_registration_matches_r[register_fd_growth_hgtf_to_mean-regfd_coefs]` | `register` | 1.00e-05 | 2.94e-01 | 1.60e-03 | xfail |
| registration | `test_registration_matches_r[register_fd_growth_hgtf_to_mean-warpfd_coefs]` | `register` | 1.00e-05 | 1.35e-01 | 7.49e-03 | xfail |
| registration | `test_registration_matches_r[register_fd_growth_hgtf_to_mean-Wfd_coefs]` | `register` | 1.00e-05 | 8.91e-02 | 5.91e-02 | xfail |
| registration | `test_registration_matches_r[register_fd_growth_hgtf_to_mean-shift]` | `register` | 1.00e-05 | 0 | 0 | pass |
| registration | `test_registration_matches_r[register_fd_weather_periodic_crit1-regfd_coefs]` | `register` | 1.00e-05 | 7.37e+00 | 5.13e-02 | xfail |
| registration | `test_registration_matches_r[register_fd_weather_periodic_crit1-warpfd_coefs]` | `register` | 1.00e-05 | 3.80e+01 | 1.01e-01 | xfail |
| registration | `test_registration_matches_r[register_fd_weather_periodic_crit1-Wfd_coefs]` | `register` | 1.00e-05 | 4.11e+00 | 2.27e+00 | xfail |
| registration | `test_registration_matches_r[register_fd_weather_periodic_crit1-shift]` | `register` | 1.00e-05 | 2.39e+00 | 1.28e-01 | xfail |
| registration | `test_warps_from_r_latent_match_r[register_fd_growth_hgtf_to_mean-warpfd_coefs]` | `register` | 1.00e-05 | 5.55e-05 | 3.07e-06 | pass |
| registration | `test_warps_from_r_latent_match_r[register_fd_growth_hgtf_to_mean-regfd_coefs]` | `register` | 1.00e-05 | 6.56e-02 | 3.58e-04 | xfail |
| registration | `test_warps_from_r_latent_match_r[register_fd_weather_periodic_crit1-warpfd_coefs]` | `register` | 1.00e-05 | 8.10e-04 | 2.16e-06 | xfail |
| registration | `test_warps_from_r_latent_match_r[register_fd_weather_periodic_crit1-regfd_coefs]` | `register` | 1.00e-05 | 4.59e-02 | 3.19e-04 | xfail |
| registration | `test_from_r_latent_keeps_r_latent_and_shift` | `register` | 1.00e-14 | 0 | 0 | pass |
| registration | `test_decompose_matches_ampphasedecomp[MS.amp]` | `register` | 1.00e-05 | 3.78e-05 | 8.55e-08 | pass |
| registration | `test_decompose_matches_ampphasedecomp[MS.pha]` | `register` | 1.00e-05 | 1.22e-04 | 1.94e-07 | pass |
| registration | `test_decompose_matches_ampphasedecomp[RSQR]` | `register` | 1.00e-05 | 6.76e-08 | 1.15e-07 | pass |
| registration | `test_decompose_matches_ampphasedecomp[C]` | `register` | 1.00e-05 | 2.34e-10 | 2.34e-10 | pass |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_crit2-regfd_coefs]` | `register` | 1.00e-05 | 2.70e+00 | 1.80e-02 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_crit2-warpfd_coefs]` | `register` | 1.00e-05 | 8.48e-01 | 4.17e-02 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_crit2-Wfd_coefs]` | `register` | 1.00e-05 | 2.01e+00 | 4.03e-01 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_crit2-shift]` | `register` | 1.00e-05 | 0 | 0 | pass |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-regfd_coefs]` | `register` | 1.00e-05 | 1.86e+01 | 1.17e-01 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-warpfd_coefs]` | `register` | 1.00e-05 | 3.95e+00 | 1.79e-01 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-Wfd_coefs]` | `register` | 1.00e-05 | 7.65e+00 | 1.22e+00 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_gait_hip_knee_periodic_crit1-shift]` | `register` | 1.00e-05 | 2.32e+00 | 9.17e-01 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_handwriting_xy_crit2-regfd_coefs]` | `register` | 1.00e-05 | 2.61e-01 | 2.68e-03 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_handwriting_xy_crit2-warpfd_coefs]` | `register` | 1.00e-05 | 3.64e-01 | 1.58e-04 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_handwriting_xy_crit2-Wfd_coefs]` | `register` | 1.00e-05 | 7.26e-03 | 9.86e-02 | xfail |
| registration_multivariate | `test_register_matches_r[register_fd_handwriting_xy_crit2-shift]` | `register` | 1.00e-05 | 0 | 0 | pass |
| registration_multivariate | `test_r_fits_multivariate_warps_to_the_first_variable[register_fd_gait_hip_knee_crit2]` | `register` | 1.00e-08 | 5.36e-10 | 1.07e-10 | pass |
| registration_multivariate | `test_r_fits_multivariate_warps_to_the_first_variable[register_fd_gait_hip_knee_periodic_crit1]` | `register` | 1.00e-08 | 2.60e-12 | 4.16e-13 | pass |
| registration_multivariate | `test_r_fits_multivariate_warps_to_the_first_variable[register_fd_handwriting_xy_crit2]` | `register` | 1.00e-08 | 2.31e-13 | 3.14e-12 | pass |
| registration_multivariate | `test_first_variable_weights_are_the_univariate_registration[register_fd_gait_hip_knee_crit2]` | `register` | 0 | 7.90e-14 | 1.61e-15 | pass |
| registration_multivariate | `test_first_variable_weights_are_the_univariate_registration[register_fd_gait_hip_knee_periodic_crit1]` | `register` | 0 | 1.85e-13 | 4.03e-15 | pass |
| registration_multivariate | `test_first_variable_weights_are_the_univariate_registration[register_fd_handwriting_xy_crit2]` | `register` | 0 | 4.50e-13 | 4.31e-14 | pass |
| registration_multivariate | `test_fdatools_optimum_is_below_r_point[register_fd_gait_hip_knee_crit2]` | `register` | — | — | — | pass |
| registration_multivariate | `test_fdatools_optimum_is_below_r_point[register_fd_gait_hip_knee_periodic_crit1]` | `register` | — | — | — | pass |
| registration_multivariate | `test_fdatools_optimum_is_below_r_point[register_fd_handwriting_xy_crit2]` | `register` | — | — | — | pass |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_crit2-warpfd_coefs]` | `register` | 1.00e-05 | 1.62e-04 | 7.96e-06 | xfail |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_crit2-regfd_coefs]` | `register` | 1.00e-05 | 6.98e-01 | 4.66e-03 | xfail |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_periodic_crit1-warpfd_coefs]` | `register` | 1.00e-05 | 5.07e-05 | 2.30e-06 | xfail |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_gait_hip_knee_periodic_crit1-regfd_coefs]` | `register` | 1.00e-05 | 6.35e-01 | 3.99e-03 | xfail |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_handwriting_xy_crit2-warpfd_coefs]` | `register` | 1.00e-05 | 1.14e-04 | 4.95e-08 | xfail |
| registration_multivariate | `test_warps_from_r_latent_match_r[register_fd_handwriting_xy_crit2-regfd_coefs]` | `register` | 1.00e-05 | 9.94e-03 | 1.02e-04 | xfail |
| registration_multivariate | `test_from_r_latent_keeps_r_latent_and_shift` | `register` | 1.00e-14 | 0 | 0 | pass |
| registration_multivariate | `test_apply_matches_register_newfd[register_fd_gait_hip_knee_crit2]` | `register` | 1.00e-05 | 2.20e-03 | 1.47e-05 | xfail |
| registration_multivariate | `test_apply_matches_register_newfd[register_fd_gait_hip_knee_periodic_crit1]` | `register` | 1.00e-05 | 7.02e-04 | 5.02e-06 | xfail |
| registration_multivariate | `test_apply_matches_register_newfd[register_fd_handwriting_xy_crit2]` | `register` | 1.00e-05 | 5.39e-05 | 5.53e-07 | xfail |
| registration_multivariate | `test_r_warps_use_a_trapezoid_rule[register_fd_gait_hip_knee_crit2]` | `register` | 1.00e-12 | 4.55e-13 | 4.81e-15 | pass |
| registration_multivariate | `test_r_warps_use_a_trapezoid_rule[register_fd_gait_hip_knee_periodic_crit1]` | `register` | 1.00e-12 | 5.97e-13 | 8.22e-15 | pass |
| registration_multivariate | `test_r_warps_use_a_trapezoid_rule[register_fd_handwriting_xy_crit2]` | `register` | 1.00e-12 | 9.32e-12 | 2.37e-14 | pass |
| registration_multivariate | `test_landmark_latent_matches_r[Wfd_coefs_hip]` | `register` | 1.00e-05 | 3.35e-06 | 1.93e-05 | xfail |
| registration_multivariate | `test_landmark_latent_matches_r[Wfd_coefs_knee]` | `register` | 1.00e-05 | 3.35e-06 | 1.93e-05 | xfail |
| registration_multivariate | `test_r_landmark_warps_do_not_depend_on_the_variable` | `register` | 1.00e-14 | 0 | 0 | pass |
| registration_multivariate | `test_fdatools_landmark_warps_meet_the_landmarks` | `register` | 0 | 3.35e-06 | 1.93e-05 | pass |
| registration_multivariate | `test_multivariate_landmarks_register_every_variable` | `register` | 0 | 5.68e-14 | 3.65e-16 | pass |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-df]` | `fregress` | 1.00e-08 | 2.66e-15 | 4.44e-16 | pass |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-yhatfdobj]` | `fregress` | 1.00e-08 | 8.86e-09 | 2.79e-09 | pass |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-gcv]` | `fregress` | 1.00e-08 | 2.93e-18 | 4.39e-15 | pass |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-OCV]` | `fregress` | 1.00e-08 | 2.41e-08 | 1.90e-08 | xfail |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-Cmat]` | `fregress` | 1.00e-08 | 5.56e-03 | 1.19e-09 | xfail |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-Dmat]` | `fregress` | 1.00e-08 | 4.60e-05 | 1.28e-09 | xfail |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-beta0]` | `fregress` | 1.00e-08 | 6.40e-11 | 6.74e-09 | pass |
| regression | `test_golden_field[fregress_scalar_precip_on_temp-beta1]` | `fregress` | 1.00e-08 | 2.78e-10 | 1.73e-08 | xfail |
| regression | `test_golden_field[predict_fregress_scalar_precip_on_temp-predicted]` | `fregress` | 1.00e-08 | 8.86e-09 | 2.79e-09 | pass |
| regression | `test_golden_field[fregress_cv_scalar_precip_on_temp-SSE.CV]` | `fregress` | 1.00e-08 | 2.33e-08 | 1.84e-08 | xfail |
| regression | `test_golden_field[fregress_cv_scalar_precip_on_temp-errfd.cv]` | `fregress` | 1.00e-08 | 2.94e-08 | 5.77e-08 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-yhatfdobj_coefs]` | `fregress` | 1.00e-08 | 2.75e-04 | 1.13e-06 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-Cmat]` | `fregress` | 1.00e-08 | 1.69e-05 | 9.50e-07 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-Dmat]` | `fregress` | 1.00e-08 | 3.93e-03 | 1.38e-06 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-beta0]` | `fregress` | 1.00e-08 | 2.10e-04 | 1.20e-06 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-beta1]` | `fregress` | 1.00e-08 | 2.75e-04 | 1.16e-06 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-beta2]` | `fregress` | 1.00e-08 | 1.55e-04 | 1.03e-06 | xfail |
| regression | `test_golden_field[fregress_functional_temp_on_region-beta3]` | `fregress` | 1.00e-08 | 2.01e-04 | 8.31e-07 | xfail |
| regression | `test_golden_field[fregress_concurrent_logprecip_on_temp-yhatfdobj_coefs]` | `fregress` | 1.00e-08 | 1.22e-01 | 6.22e-03 | xfail |
| regression | `test_golden_field[fregress_concurrent_logprecip_on_temp-Cmat]` | `fregress` | 1.00e-08 | 8.00e+00 | 1.07e-03 | xfail |
| regression | `test_golden_field[fregress_concurrent_logprecip_on_temp-Dmat]` | `fregress` | 1.00e-08 | 1.43e+01 | 3.44e-03 | xfail |
| regression | `test_golden_field[fregress_concurrent_logprecip_on_temp-beta0]` | `fregress` | 1.00e-08 | 1.20e-01 | 1.19e-02 | xfail |
| regression | `test_golden_field[fregress_concurrent_logprecip_on_temp-beta1]` | `fregress` | 1.00e-08 | 8.68e-03 | 6.61e-03 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-df]` | `fregress` | 1.00e-08 | 7.82e-06 | 1.51e-06 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-yhatfdobj]` | `fregress` | 1.00e-08 | 8.38e-05 | 3.06e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-gcv]` | `fregress` | 1.00e-08 | 1.68e-09 | 1.96e-06 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-OCV]` | `fregress` | 1.00e-08 | 1.11e-05 | 3.02e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-Cmat]` | `fregress` | 1.00e-08 | 4.08e-04 | 1.16e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-Dmat]` | `fregress` | 1.00e-08 | 3.66e-04 | 1.01e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-beta0]` | `fregress` | 1.00e-08 | 2.17e-05 | 1.07e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-beta1]` | `fregress` | 1.00e-08 | 7.11e-06 | 1.41e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-beta2]` | `fregress` | 1.00e-08 | 1.23e-05 | 3.86e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional-beta3]` | `fregress` | 1.00e-08 | 1.10e-04 | 7.03e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-df]` | `fregress` | 1.00e-08 | 1.72e-13 | 1.57e-14 | pass |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-yhatfdobj]` | `fregress` | 1.00e-08 | 1.93e-04 | 7.00e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-gcv]` | `fregress` | 1.00e-08 | 1.94e-09 | 1.20e-06 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-OCV]` | `fregress` | 1.00e-08 | 2.37e-04 | 3.63e-04 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-Cmat]` | `fregress` | 1.00e-08 | 4.08e-04 | 1.16e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-Dmat]` | `fregress` | 1.00e-08 | 3.66e-04 | 1.01e-05 | xfail |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-beta0]` | `fregress` | 1.00e-08 | 5.33e-15 | 2.65e-15 | pass |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-beta1]` | `fregress` | 1.00e-08 | 8.66e-15 | 1.69e-14 | pass |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-beta2]` | `fregress` | 1.00e-08 | 7.22e-15 | 2.25e-14 | pass |
| regression | `test_golden_field[fregress_synthetic_2scalar_1functional_lambda0-beta3]` | `fregress` | 1.00e-08 | 4.51e-03 | 2.08e-03 | xfail |
| regression | `test_null_golden_cases_really_are_empty[fregress_stderr_functional_temp_on_region]` | `fregress` | — | — | — | pass |
| regression | `test_null_golden_cases_really_are_empty[fregress_stderr_scalar_precip_on_temp]` | `fregress` | — | — | — | pass |
| regression | `test_every_defect_names_a_real_field` | `fregress` | — | — | — | pass |
| regression | `test_scalar_design_matches_r[Cmat]` | `fregress` | 1.00e-08 | 0 | 0 | pass |
| regression | `test_scalar_design_matches_r[Dmat]` | `fregress` | 1.00e-08 | 0 | 0 | pass |
| regression | `test_scalar_design_matches_r[beta]` | `fregress` | 1.00e-08 | 1.78e-15 | 9.45e-16 | pass |
| regression | `test_scalar_design_matches_r[betastderr]` | `fregress` | 1.00e-08 | 7.22e-16 | 3.97e-15 | pass |
| regression | `test_scalar_design_matches_r[bvar]` | `fregress` | 1.00e-08 | 9.02e-17 | 2.73e-15 | pass |
| regression | `test_scalar_design_matches_r[df]` | `fregress` | 1.00e-08 | 4.44e-16 | 1.48e-16 | pass |
| regression | `test_scalar_design_matches_r[errfd_cv]` | `fregress` | 1.00e-08 | 3.91e-14 | 4.93e-14 | pass |
| regression | `test_scalar_design_matches_r[gcv]` | `fregress` | 1.00e-08 | 8.67e-18 | 2.53e-15 | pass |
| regression | `test_scalar_design_matches_r[ocv]` | `fregress` | 1.00e-08 | 5.77e-15 | 1.86e-15 | pass |
| regression | `test_scalar_design_matches_r[predicted]` | `fregress` | 1.00e-08 | 2.00e-15 | 6.44e-16 | pass |
| regression | `test_scalar_design_matches_r[sigma2]` | `fregress` | 1.00e-08 | 2.36e-16 | 2.55e-15 | pass |
| regression | `test_scalar_design_matches_r[sse_cv]` | `fregress` | 1.00e-08 | 6.00e-14 | 1.93e-14 | pass |
| regression | `test_scalar_design_matches_r[yhat]` | `fregress` | 1.00e-08 | 5.33e-15 | 1.25e-15 | pass |
| regression | `test_scalar_design_constant_basis_is_what_r_uses` | `fregress` | — | — | — | pass |
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
| stats | `test_r_generator_matches_r` | `FData` | — | — | — | pass |
| stats | `test_var_fd` | `stats.cov` | 1.00e-08 | 4.41e-13 | 4.70e-15 | pass |
| stats | `test_cor_fd` | `stats.cor` | 1.00e-08 | 3.22e-15 | 4.03e-15 | pass |
| stats | `test_mean_sd_fd[mean_coefs]` | `FData` | 1.00e-08 | 8.53e-14 | 4.41e-16 | pass |
| stats | `test_mean_sd_fd[mean_fitted]` | `FData` | 1.00e-08 | 2.84e-14 | 1.67e-15 | pass |
| stats | `test_mean_sd_fd[sd_coefs]` | `FData` | 1.00e-08 | 4.26e-14 | 3.40e-16 | pass |
| stats | `test_mean_sd_fd[sd_fitted]` | `FData` | 1.00e-08 | 1.60e-14 | 1.65e-15 | pass |
| stats | `test_fbplot_mbd` | `stats.boxplot` | 1.00e-08 | 7.33e-15 | 1.38e-14 | pass |
| stats | `test_fbplot_input_is_the_smoothed_weather` | `stats.boxplot` | 1.00e-08 | 7.82e-14 | 2.33e-15 | pass |
| stats | `test_fdepth_fm[prof]` | `stats.depth` | 1.00e-08 | 2.22e-16 | 2.35e-16 | pass |
| stats | `test_fdepth_fm[median]` | `stats.depth` | 1.00e-08 | 0 | 0 | pass |
| stats | `test_fdepth_fm[lmed]` | `stats.depth` | — | — | — | pass |
| stats | `test_fdepth_fm[mtrim]` | `stats.depth` | 1.00e-08 | 7.11e-15 | 4.01e-16 | pass |
| stats | `test_fdepth_fm[ltrim]` | `stats.depth` | — | — | — | pass |
| stats | `test_fdepth_accepts_fdata` | `stats.depth` | 1.00e-08 | 2.22e-16 | 2.35e-16 | pass |
| stats | `test_permutation_tests[tperm_fd_weather_atlantic_vs_pacific-Tobs]` | `stats.t_test` | 1.00e-08 | 6.22e-15 | 1.20e-15 | pass |
| stats | `test_permutation_tests[tperm_fd_weather_atlantic_vs_pacific-Tnull]` | `stats.t_test` | 1.00e-08 | 3.38e-14 | 4.69e-15 | pass |
| stats | `test_permutation_tests[tperm_fd_weather_atlantic_vs_pacific-pval]` | `stats.t_test` | 1.00e-08 | 0 | 0 | pass |
| stats | `test_permutation_tests[tperm_fd_weather_atlantic_vs_pacific-qval]` | `stats.t_test` | 1.00e-08 | 8.44e-15 | 1.44e-14 | pass |
| stats | `test_permutation_tests[tperm_fd_weather_atlantic_vs_pacific-argvals]` | `stats.t_test` | 1.00e-08 | 0 | 0 | pass |
| stats | `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-Fobs]` | `stats.f_test` | 1.00e-08 | 4.36e-05 | 9.85e-05 | xfail |
| stats | `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-Fnull]` | `stats.f_test` | 1.00e-08 | 1.67e-05 | 4.65e-05 | xfail |
| stats | `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-pval]` | `stats.f_test` | 1.00e-08 | 0 | 0 | pass |
| stats | `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-qval]` | `stats.f_test` | 1.00e-08 | 6.58e-07 | 1.11e-04 | xfail |
| stats | `test_permutation_tests[fperm_fd_weather_temp_atlantic_dummy-argvals]` | `stats.f_test` | 1.00e-08 | 0 | 0 | pass |
| stats | `test_fperm_through_a_fitted_model` | `FData` | 1.00e-08 | 0 | 0 | pass |
