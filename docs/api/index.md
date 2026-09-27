# API reference

One page per module of Fabel's public API. Every name in `fabel.__all__`, every
public function and class of `fabel.stats`, `fabel.datasets`, `fabel.io` and
`fabel.nn`, and every result class is documented here. For worked examples,
see the [tutorials](../tutorials/index.md).

| Module | Contents |
|---|---|
| [Basis systems](basis.md) | `Basis`, `BSpline`, `Fourier`, `Monomial`, `Exponential`, `Power`, `Constant`, `Polygonal` |
| [Core objects](core.md) | `FData`, `BiFData`, `LDO`, `inprod` |
| [Smoothing](smoothing.md) | `smooth`, `SmoothResult`, `Smoother`, `gcv_curve`, `lambda_to_df`, `df_to_lambda` |
| [Decomposition](decomposition.md) | `FPCA`, `FCCA` |
| [Regression](regression.md) | `fregress`, `FRegress`, `FRegressResult`, `FRegressStderr`, `FRegressCV` |
| [Registration](registration.md) | `register`, `landmark_register`, `RegistrationResult`, `AmpPhaseDecomposition`, `Registrator` |
| [Dynamics](dynamics.md) | `PDA`, `phase_plane` |
| [Statistics](stats.md) | `cov`, `cor`, `depth`, `boxplot`, `t_test`, `f_test`, `DepthResult`, `BoxplotResult`, `PermutationTestResult` |
| [PyTorch layers](nn.md) | `BasisLayer`, `SmoothingLayer`, `FDataDataset` (optional extra) |
| [Datasets](datasets.md) | 14 `load_*` loaders and the dataclasses they return (`Growth`, `CanadianWeather`, `Lip`, ...) |
| [Input and output](io.md) | `from_pandas`, `LongData`, `to_pandas`, `to_xarray`, `read_rds` |
