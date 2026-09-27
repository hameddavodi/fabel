# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow SemVer.

## [1.0.0] - Unreleased

First public release: a clean-room Python rewrite of R `fda` 6.3.0 with
golden-file parity (`rtol = 1e-8`, `1e-5` for iterative fits).

### Added
- Basis systems (`fabel.basis`): `BSpline`, `Fourier`, `Monomial`, `Exponential`,
  `Power`, `Constant`, `Polygonal`, with evaluation and derivatives, roughness
  penalties for any `LDO`, cached Gram matrices and exact basis products.
- Core objects (`fabel.core`): `FData` (callable curves, exact derivatives,
  arithmetic, `mean` / `std` / `center` / `cov`, indexing, `@` inner product),
  `BiFData`, `LDO` (including the harmonic accelerator) and `inprod`.
- Smoothing (`fabel.smoothing`): `smooth()` with GCV or degrees-of-freedom
  selection of λ, positive / monotone / morph constraints, observation weights
  and irregular per-curve designs; `SmoothResult`; the scikit-learn `Smoother`;
  `gcv_curve`, `lambda_to_df`, `df_to_lambda`.
- Decomposition (`fabel.decomposition`): `FPCA` (with roughness-penalised
  harmonics and varimax rotation) and `FCCA`, both scikit-learn estimators.
- Regression (`fabel.regression`): `fregress()` for scalar and functional
  responses with scalar and functional covariates (model type read from the
  arguments, or an R-style formula string with treatment-coded factors), with
  `predict()`, `stderr()` and `cv()`; the scikit-learn `FRegress` estimator.
- Registration (`fabel.registration`): continuous registration `register()`
  (Newton with the exact Hessian, optional periodic shift), landmark registration
  (`landmarks=` or `landmark_register()`), the amplitude/phase decomposition
  `RegistrationResult.decompose()` and the scikit-learn `Registrator`.
- Dynamics (`fabel.dynamics`): principal differential analysis `PDA` for single
  equations and coupled systems, with an ODE `solve()` and `plot_overlay()`, and
  the `phase_plane()` plot.
- Statistics (`fabel.stats`): `cov`, `cor`, functional depth (MBD, BD2, FM),
  the functional `boxplot`, and the permutation tests `t_test` and `f_test`.
  `f_test` takes either raw inputs `(y, x, basis=, lam=, penalty=)` like R
  `Fperm.fd`, or a fitted `fregress` model: `f_test(model, n_perm=, q=, t=,
  random_state=)`.
- scikit-learn: `check_estimator` passes for `Smoother`, `FPCA`, `FRegress` and
  `Registrator` with no exemptions. `FRegress` validates `y` as scikit-learn does
  (finite, column vectors flattened with a warning, at least 2 samples).
  `Registrator` has `n_iter_`, and its default basis for an n-column coefficient
  matrix is `BSpline(n_basis=n, order=min(4, n))`, the same rule as `FPCA`.
- PyTorch layers (`fabel.nn`, optional `fabel[torch]` extra): `BasisLayer`,
  `SmoothingLayer` (learnable λ) and `FDataDataset`. `import fabel` does not
  import PyTorch; `fabel.nn` loads on first use.
- Top-level exports: `smooth`, `Smoother`, `SmoothResult`, `FPCA`, `FCCA`,
  `fregress`, `FRegress`, `register`, `landmark_register`, `Registrator`, `PDA`,
  `phase_plane`, and the `stats` module.
- Datasets (`fabel.datasets`): 14 loaders for the FDA book datasets; `growth`,
  `gait` and `pinch` ship in the package, the others download once with SHA-256
  verification and a `FABEL_DATA_DIR`-overridable cache.
- I/O (`fabel.io`): `from_pandas`, `to_pandas`, `to_xarray` and `read_rds` for R
  `fd` / `bifd` / `basisfd` objects.
- Array API backend: NumPy and PyTorch inputs, with gradients flowing through
  evaluation, products and inner products.
- Golden-file parity suite against R `fda` 6.3.0, and `tools/parity_report.py`,
  which measures every parity check and writes `PARITY_REPORT.md`.
- Documentation site (MkDocs Material + mkdocstrings): quickstart, R migration
  table, API reference for every module, and six tutorials (smoothing, FPCA,
  registration, regression, dynamics, machine learning). Tests run every
  tutorial block and check that every public symbol is rendered.
- `notebooks/book_figures.ipynb`: figures of Ramsay, Hooker & Graves (2009),
  assembled from `notebooks/book/ch*.py` by `tools/build_book_notebook.py`
  (`make book`) and run in CI with `nbmake`.
- Packaging: typed (`py.typed`) wheel and sdist for Python 3.10 to 3.13,
  BSD-3-Clause license, `CITATION.cff`.

### Fixed
- `FData.std` samples the pointwise standard deviation on
  `max(201, 10 * n_basis + 1)` points, matching R `sd.fd` to 1.8e-15.

### Known differences from R `fda`
- Where R `fda` 6.3.0 is demonstrably less accurate (for example `deriv.fd`,
  `times.fd`, B-spline penalties with no interior knots, `smooth.pos` stopping
  one iteration short), Fabel returns the exact value. Each case is a strict
  expected failure in the test suite with its measured error; see
  `PARITY_REPORT.md`.
