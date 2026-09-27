# API reference

One page per module of Fabel's public API.

| Module | Contents |
|---|---|
| [Basis systems](basis.md) | `Basis`, `BSpline`, `Fourier`, `Monomial`, `Exponential`, `Power`, `Constant`, `Polygonal` |
| [Core objects](core.md) | `FData`, `BiFData`, `LDO`, `inprod` |
| [Smoothing](smoothing.md) | `smooth`, `SmoothResult`, `Smoother`, `gcv_curve`, `lambda_to_df`, `df_to_lambda` |
| [Decomposition](decomposition.md) | `FPCA`, `FCCA` |
| [Regression](regression.md) | `fregress`, `FRegress`, `FRegressResult`, `FRegressStderr`, `FRegressCV` |
| [Registration](registration.md) | `register`, `landmark_register`, `RegistrationResult`, `AmpPhaseDecomposition`, `Registrator` |
| [Dynamics](dynamics.md) | `PDA`, `phase_plane` |
| [Statistics](stats.md) | `cov`, `cor`, `depth`, `boxplot`, `t_test`, `f_test` and their result types |
| [PyTorch layers](nn.md) | `BasisLayer`, `SmoothingLayer`, `FDataDataset` (optional extra) |
| [Datasets](datasets.md) | `load_*` loaders and their result types |
| [Input and output](io.md) | `from_pandas`, `to_pandas`, `to_xarray`, `read_rds` |
