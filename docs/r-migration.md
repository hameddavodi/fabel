# Migrating from R `fda`

Fabel replaces the 515 functions of R's `fda` 6.3.0 with about 40 Python
symbols. The ideas stay the same: a curve is a set of basis coefficients, a
roughness penalty is a linear differential operator, and every estimator works
on those coefficients. What changes is the surface:

- **Objects are callable.** `fd(t)` evaluates a curve. There is no `eval.*` family.
- **Operators are math.** `fd1 + fd2`, `fd * 2`, `fd1 @ fd2` (inner product), `basis1 * basis2`.
- **One entry point per task.** One `smooth()`, one `register()`, one `.plot()`.
- **Immutable objects.** Every method returns a new object, so chaining is safe.
- **scikit-learn estimators.** `Smoother`, `FPCA`, `FCCA`, `FRegress` and
  `Registrator` implement `fit` / `transform` / `predict`.

In the tables below, `fb` is `import fabel as fb`.

## Functional data objects (`fd` → `FData`)

| R `fda` | Fabel |
|---|---|
| `fd(coef, basisobj)` | `fb.FData(coefs, basis)` |
| `eval.fd`, `eval.monfd`, `eval.posfd` | `fd(t, deriv=0)` |
| `eval.fd(t, fdobj, Lfdobj)` | `fd(t, operator)` with an `LDO` |
| `deriv.fd` | `fd.derivative(n=1)` |
| `mean.fd`, `sd.fd`, `stddev.fd`, `center.fd` | `fd.mean()`, `fd.std()`, `fd.center()` |
| `arithmetic.fd`, `sum.fd` (`+`, `-`, `*`, `^`) | `fd + fd`, `fd * 2`, `fd ** 2` |
| `inprod`, `inprod.bspline` | `fd1 @ fd2`, `fb.inprod(a, b, lfd1=, lfd2=)` |
| `var.fd` | `fd.cov()` |
| `plot.fd`, `lines.fd`, `plotfit.fd` | `fd.plot()`, `fd.plot_fit(y, t)` |
| `[.fd`, `subscript.fd` | `fd[i]`, `fd[i:j]` |
| `fdlabels`, `norder`, `nbasis` | `len(fd)`, `fd.n_curves`, `fd.domain`, `fd.basis` |
| (none) | `fd.to_numpy(t)`, `fd.to_torch(t)` |
| `bifd`, `eval.bifd` | `fb.BiFData(coefs, s_basis, t_basis)`, `bifd(s, t)` |

## Basis systems (`create.*.basis` → `Basis` subclasses)

| R `fda` | Fabel |
|---|---|
| `create.bspline.basis` | `fb.BSpline(domain=, n_basis=, order=, breaks=)` |
| `create.fourier.basis` | `fb.Fourier(domain=, n_basis=, period=)` |
| `create.monomial.basis` | `fb.Monomial(domain=, n_basis=)` |
| `create.exponential.basis` | `fb.Exponential(domain=, rates=)` |
| `create.power.basis` | `fb.Power(domain=, exponents=)` |
| `create.constant.basis` | `fb.Constant(domain=)` |
| `create.polygonal.basis` | `fb.Polygonal(argvals)` |
| `eval.basis`, `getbasismatrix`, `bsplineS`, `fourier`, `expon`, `monomial`, `polyg`, `powerbasis` | `basis(t, deriv=0)` |
| `getbasispenalty`, `eval.penalty`, `bsplinepen`, `fourierpen`, `polygpen`, `powerpen`, `polynompen` | `basis.penalty(op=2)` |
| `basisfd.product` | `basis1 * basis2` |
| internal `inprod` calls on a basis | `basis.gram()` (cached) |

## Linear differential operators (`Lfd` → `LDO`)

| R `fda` | Fabel |
|---|---|
| `int2Lfd(m)` | `fb.LDO(m)` |
| `vec2Lfd`, `Lfd` | `fb.LDO(weights=[w0, w1, ...])` |
| harmonic acceleration `vec2Lfd(c(0, (2*pi/T)^2, 0), rng)` | `fb.LDO.harmonic(period=T)` |

## Smoothing (`smooth.*` → `smooth()`)

| R `fda` | Fabel |
|---|---|
| `smooth.basis`, `smooth.basisPar`, `Data2fd` | `fb.smooth(y, t, basis=, lam=, penalty=)` |
| `smooth.monotone` | `fb.smooth(y, t, constraint="monotone")` |
| `smooth.pos` | `fb.smooth(y, t, constraint="positive")` |
| `smooth.morph` | `fb.smooth(y, t, constraint="morph")` |
| `lambda2gcv`, GCV search over λ | `fb.smooth(y, t, lam="gcv")`, `gcv_curve(...)` |
| `df2lambda` | `fb.smooth(y, t, lam="df=12")`, `df_to_lambda(...)` |
| `lambda2df` | `lambda_to_df(...)`, `result.df` |
| `smooth.basis.sparse` (irregular sampling per curve) | `fb.smooth([y1, y2, ...], [t1, t2, ...])` |
| `smooth.bibasis` | `fb.smooth()` with a 2-D `y` |
| the returned list (`fd`, `df`, `gcv`, `SSE`, `penmat`, `y2cMap`) | `SmoothResult`: `.fd`, `.df`, `.gcv`, `.sse`, `.penalty_matrix`, `.y2c_map` |
| (none) | `fb.Smoother()`: scikit-learn transformer |

## Decomposition (`pca.fd`, `cca.fd` → `FPCA`, `FCCA`)

| R `fda` | Fabel |
|---|---|
| `pca.fd(fdobj, nharm, harmfdPar)` | `fb.FPCA(n=4, lam="gcv").fit(fd)` |
| `pca.fd(...)$harmonics` / `$scores` / `$values` / `$varprop` | `pca.harmonics`, `pca.scores`, `pca.values`, `pca.varprop` |
| scores of new curves | `pca.transform(fd)` |
| `plot.pca.fd` | `pca.harmonics.plot()` |
| `varmx.pca.fd` | `pca.rotate("varimax")` |
| `cca.fd` | `fb.FCCA(n=3).fit(fd1, fd2)` |

## Regression (`fRegress` family → `fregress()` / `FRegress`)

| R `fda` | Fabel |
|---|---|
| `fRegress` (scalar or functional response) | `fb.fregress(y, x_list)` |
| `fRegress.formula` | `fb.fregress("temp ~ region + precip", data)` |
| `predict.fRegress` | `model.predict(x_new)` |
| `fRegress.stderr` | `model.stderr()` |
| `fRegress.CV` | `model.cv()` |
| `Fperm.fd` | `fb.stats.f_test(model, n_perm=1000)` |
| (none) | `fb.FRegress`: scikit-learn estimator |

## Registration (`register.fd`, `landmarkreg` → `register()`)

| R `fda` | Fabel |
|---|---|
| `register.fd` (continuous registration to the mean) | `fb.register(fd)` |
| `landmarkreg` | `fb.register(fd, landmarks=lm)`, `fb.landmark_register(fd, landmarks)` |
| `$regfd`, `$warpfd` | `result.registered`, `result.warp` |
| `AmpPhaseDecomp` | `result.decompose()` → amplitude MSE, phase MSE, R² |
| (none) | `fb.Registrator`: scikit-learn transformer |

## Dynamics (`pda.fd` → `PDA`)

| R `fda` | Fabel |
|---|---|
| `pda.fd` | `fb.PDA(order=2).fit(fd)` |
| `pda.overlay` | `pda.plot_overlay()` |
| `phaseplanePlot` | `fb.phase_plane(fd)` |

## Statistics (`fabel.stats`)

| R `fda` | Fabel |
|---|---|
| `var.fd` | `fb.stats.cov(fd)` |
| `cor.fd` | `fb.stats.cor(fd1, fd2)` |
| `fdepth` | `fb.stats.depth(fd)` |
| `fbplot`, `boxplot.fd` | `fb.stats.boxplot(fd)` |
| `Fperm.fd` | `fb.stats.f_test(model, n_perm=)` |
| `tperm.fd` | `fb.stats.t_test(fd1, fd2, n_perm=)` |

## Datasets (`data(...)` → `fabel.datasets`)

| R `fda` | Fabel |
|---|---|
| `CanadianWeather` | `fb.datasets.load_canadian_weather()` |
| `growth` | `fb.datasets.load_growth()` |
| `gait` | `fb.datasets.load_gait()` |
| `handwrit`, `handwritTime` | `fb.datasets.load_handwriting()` |
| `pinch`, `pinchtime`, `pinchraw` | `fb.datasets.load_pinch()` |
| `melanoma` | `fb.datasets.load_melanoma()` |
| `refinery` | `fb.datasets.load_refinery()` |
| `seabird` | `fb.datasets.load_seabird()` |
| `ReginaPrecip` | `fb.datasets.load_regina_precip()` |
| `MontrealTemp` | `fb.datasets.load_montreal_temp()` |
| `daily` | `fb.datasets.load_daily()` |
| `infantGrowth` | `fb.datasets.load_infant_growth()` |
| `nondurables` | `fb.datasets.load_nondurables()` |
| `lip`, `lipmarks`, `liptime` | `fb.datasets.load_lip()` |

`growth`, `gait` and `pinch` ship inside the package. The others are
downloaded once from the project's GitHub release, checked against a SHA-256
checksum, and cached in `~/.cache/fabel` (override with `FABEL_DATA_DIR`).

## Ecosystem I/O and deep learning

| R `fda` | Fabel |
|---|---|
| `readRDS` on a saved `fd`, `bifd` or `basisfd` | `fb.read_rds("weather.rds")` |
| (none) | `fb.from_pandas(df, id_col, t_col, y_col)` (long format) |
| (none) | `fd.to_pandas(t)`, `fd.to_xarray(t)` |
| (none) | `fb.nn.BasisLayer(basis)`: differentiable `fd(t)` for PyTorch |
| (none) | `fb.nn.FDataDataset(fd, labels)`: `DataLoader`-ready dataset |

## Not carried over

- R plotting quirks (`axisIntervals`, `matplot`): Fabel uses matplotlib defaults.
- `readHMD` (Human Mortality Database scraper).
- Internal helpers (`wtcheck`, `symsolve`, `polintmat`, `lnsrch`, ...): private
  in Fabel or replaced by SciPy.

## Behaviour differences from R

| Topic | R `fda` 6.3.0 | Fabel |
|---|---|---|
| Weights in `fRegress` | `fRegress(y, xfdlist, betalist, wt = w)` fits by weighted least squares, scalar and functional response alike. The argument is `wt`: `wtvec = w` falls into `...` and is ignored without a warning, so R returns the unweighted fit. | `fb.fregress(y, x, weights=w)` fits by weighted least squares and agrees with R's `wt = w` to the accuracy of R's integration. Drop `weights` when porting a script that passed `wtvec` to `fRegress`. |
| Weights in `Fperm.fd` | `Fperm.fd(..., wt = w)` accepts weights but ignores them: `Fobs` and the null distribution are the same as with no weights. | `fb.stats.f_test(model)` refits a weighted model with its weights under every permutation, so the statistic changes with the weights. Test a model fitted without `weights` to reproduce R. The raw form `f_test(y, x)` takes no weights. |

With a roughness penalty, only the relative size of the weights and `lam`
matters: multiplying every weight by `c` has the same effect as dividing `lam`
by `c`. With unit weights both functions give exactly the unweighted result.

## Where Fabel and R disagree on purpose

Fabel matches R to `rtol = 1e-8` (`1e-5` for iterative fits) on every golden
case except a documented set where R is the less accurate side, for example
`deriv.fd` re-projecting a spline derivative (1.3% error, Fabel is exact) or
`times.fd` using a basis too smooth to hold the product (12.8% error). Each of
these is a strict expected failure in the test suite with its measured reason;
the full list is in `PARITY_REPORT.md` in the repository.
