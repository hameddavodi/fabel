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
  table and API reference.
- Packaging: typed (`py.typed`) wheel and sdist for Python 3.10 to 3.13,
  BSD-3-Clause license, `CITATION.cff`.

### Known differences from R `fda`
- Where R `fda` 6.3.0 is demonstrably less accurate (for example `deriv.fd`,
  `times.fd`, B-spline penalties with no interior knots, `smooth.pos` stopping
  one iteration short), Fabel returns the exact value. Each case is a strict
  expected failure in the test suite with its measured error; see
  `PARITY_REPORT.md`.
