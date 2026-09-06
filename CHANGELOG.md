# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow SemVer.

## [1.0.0] - 2026-09-06

### Added
- Core objects: `FData`, `BiFData`, `LDO`; basis systems `BSpline`, `Fourier`, `Monomial`,
  `Exponential`, `Power`, `Constant`, `Polygonal`.
- `smooth()` with GCV / df-based lambda selection, monotone and positive constraints,
  irregular per-curve sampling and weights.
- `FPCA`, `FCCA`, `fregress()` / `FRegress`, `register()`, `landmark_register()`, `PDA`,
  `phase_plane()`, functional statistics (`cov`, `cor`, `depth`, `boxplot`, `f_test`, `t_test`).
- Datasets from the FDA books with lazy download and checksum verification.
- Array API backend (NumPy / PyTorch), `fabel.nn` layers, scikit-learn estimator API.
- Golden-file parity suite against R `fda` 6.3.0.
