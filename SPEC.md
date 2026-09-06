# Fabel — API Specification v0.1

**Goal:** Python rewrite of R `fda` (Ramsay et al., v6.3.0).
**Promise:** 515 R functions → ~40 public symbols. Same math. Zero boilerplate.

---

## 1. Design Principles

1. **One obvious way.** One `smooth()`, one `.plot()`, one `register()`.
2. **Objects are callable.** `fd(t)` evaluates. No `eval.*` functions.
3. **Operators are math.** `fd1 + fd2`, `fd1 @ fd2` (inner product), `basis1 * basis2`.
4. **Smart defaults.** Basis size, λ (GCV), knots — all auto unless overridden.
5. **Backend-agnostic.** NumPy / PyTorch / JAX via Array API. Tensors in → tensors out. Gradients flow.
6. **Immutable objects.** Every method returns a new object → safe chaining.
7. **sklearn-compatible.** Estimators implement `fit / transform / predict`.
8. **Internals are private.** Parity is tested on public API + book figures only.

---

## 2. Package Layout

```
fabel/
├── core.py          # FData, BiFData, LDO
├── basis.py         # Basis + 7 subclasses
├── smoothing.py     # smooth(), Smoother
├── decomposition.py # FPCA, FCCA
├── regression.py    # fregress(), FRegress
├── registration.py  # register(), landmark_register()
├── dynamics.py      # PDA, phase_plane()
├── stats.py         # cov, cor, depth, boxplot, f_test, t_test
├── datasets.py      # load_*() — all 15+ book datasets
├── nn.py            # PyTorch: BasisLayer, FDataDataset
├── _backend.py      # array-api dispatch (private)
└── _linalg.py       # banded/sparse solvers, caching (private)
```

---

## 3. Core Objects

### 3.1 `FData` — replaces `fd` + 20 R functions

```python
fd = fb.FData(coefs, basis)          # explicit construction
y  = fd(t)                           # eval        ← eval.fd
y1 = fd(t, deriv=1)                  # eval deriv  ← eval.fd(Lfdobj=1)
```

| Python | Replaces (R) |
|---|---|
| `fd(t, deriv=0)` | `eval.fd`, `eval.monfd`, `eval.posfd` |
| `fd.derivative(n=1)` | `deriv.fd` |
| `fd.mean()` / `fd.std()` / `fd.center()` | `mean.fd`, `sd.fd`, `stddev.fd`, `center.fd` |
| `fd + fd`, `fd * 2`, `fd ** 2` | `arithmetic.fd`, `sum.fd` |
| `fd1 @ fd2` | `inprod`, `inprod.bspline` |
| `fd.cov()` | `var.fd` |
| `fd.plot()` / `fd.plot_fit(y, t)` | `plot.fd`, `plotfit.fd`, `lines.fd` |
| `fd[i]`, `fd[i:j]` | `subscript.fd`, `[.fd` |
| `fd.to_numpy(t)` / `fd.to_torch(t)` | — (new) |
| `len(fd)`, `fd.n_curves`, `fd.domain` | `fdlabels`, `norder`, etc. |

### 3.2 `Basis` — replaces 30+ R functions

```python
b = fb.BSpline(domain=(0, 365), n_basis=65)   # ← create.bspline.basis
b = fb.Fourier(domain=(0, 365), n_basis=65)   # ← create.fourier.basis
```

Subclasses: `BSpline, Fourier, Monomial, Exponential, Power, Constant, Polygonal`

| Python | Replaces (R) |
|---|---|
| `b(t, deriv=0)` | `eval.basis`, `getbasismatrix`, `bsplineS`, `fourier`, `expon`, `monomial`, `polyg`, `powerbasis` |
| `b.penalty(op=2)` | `getbasispenalty`, `eval.penalty`, `bsplinepen`, `fourierpen`, `polygpen`, `powerpen`, `polynompen` |
| `b1 * b2` | `basisfd.product` |
| `b.gram()` (cached) | internal `inprod` calls |

### 3.3 `LDO` — linear differential operator, replaces `Lfd`

```python
L = fb.LDO(2)                      # ← int2Lfd(2)
L = fb.LDO(weights=[w0, w1])       # ← vec2Lfd / Lfd
harm = fb.LDO.harmonic(period=365) # ← the book's harmonic accelerator
```

---

## 4. High-Level API (the one-liners)

### 4.1 `smooth()` — replaces 12 R functions

```python
fd = fb.smooth(y, t)                              # auto basis, auto λ (GCV)
fd = fb.smooth(y, t, basis=b, lam=1e2)            # full control
fd = fb.smooth(y, t, constraint="monotone")       # ← smooth.monotone
fd = fb.smooth(y, t, constraint="positive")       # ← smooth.pos
fd = fb.smooth(y, t, lam="gcv", penalty=fb.LDO(2))
```

| Argument | Replaces (R) |
|---|---|
| default call | `smooth.basis`, `smooth.basisPar`, `Data2fd` |
| `constraint=` | `smooth.monotone`, `smooth.pos`, `smooth.morph` |
| `lam="gcv"` \| `"df=12"` | `df2lambda`, `lambda2df`, `lambda2gcv` |
| sparse/irregular `t` per curve | `smooth.basis.sparse` |
| 2-D `y` | `smooth.bibasis` |

Returns `SmoothResult`: `.fd`, `.df`, `.gcv`, `.sse`, `.penalty_matrix` — everything R scatters across list elements.

### 4.2 Decomposition

```python
pca = fb.FPCA(n=4, lam="gcv").fit(fd)     # ← pca.fd
pca.harmonics.plot()                       # ← plot.pca.fd
scores = pca.transform(fd)                 # sklearn-style
pca.rotate("varimax")                      # ← varmx
cca = fb.FCCA(n=3).fit(fd1, fd2)          # ← cca.fd
```

### 4.3 Regression — replaces `fRegress` family (7 functions)

```python
m = fb.fregress(y_fd, x_list)             # scalar/functional mix auto-detected
m = fb.fregress("temp ~ region + precip", data)   # ← fRegress.formula
m.predict(x_new)                           # ← predict.fRegress
m.stderr()                                 # ← fRegress.stderr
m.cv()                                     # ← fRegress.CV
fb.stats.f_test(m, n_perm=1000)            # ← Fperm.fd
```

### 4.4 Registration — replaces 5 functions

```python
res = fb.register(fd)                      # ← register.fd (to mean)
res = fb.register(fd, landmarks=lm)        # ← landmarkreg
res.registered, res.warp                   # warped curves + warping fns
res.decompose()                            # ← AmpPhaseDecomp → (amp_mse, phase_mse, R²)
```

### 4.5 Dynamics

```python
pda = fb.PDA(order=2).fit(fd)              # ← pda.fd
pda.plot_overlay()                          # ← pda.overlay
fb.phase_plane(fd)                          # ← phaseplanePlot
```

---

## 5. ML Integration

### 5.1 scikit-learn (native, no wrappers)

`Smoother`, `FPCA`, `FCCA`, `Registrator`, `FRegress` all implement the estimator API:

```python
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV

pipe = Pipeline([
    ("smooth", fb.Smoother()),
    ("fpca",  fb.FPCA(n=5)),
    ("clf",   LogisticRegression()),
])
GridSearchCV(pipe, {"fpca__n": [3, 5, 8]}).fit(X, y)
```

### 5.2 PyTorch (`fabel.nn`)

```python
layer = fb.nn.BasisLayer(fb.BSpline(domain, 32))  # differentiable fd(t)
ds    = fb.nn.FDataDataset(fd, labels)            # → DataLoader ready
```

- All core ops accept `torch.Tensor` → return `torch.Tensor`
- Gradients flow through smoothing, eval, inner products (autodiff registration!)
- `.to("cuda")` supported via Array API dispatch

### 5.3 Ecosystem I/O

```python
fb.from_pandas(df, id_col, t_col, y_col)   # long-format → FData
fd.to_xarray() / fd.to_pandas(t)
fb.read_rds("weather.rds")                  # ingest R fda objects directly
```

---

## 6. Performance Spec

| Bottleneck in R | Fabel fix |
|---|---|
| Dense penalty matrices | scipy.sparse banded (B-spline Gram is banded) |
| Repeated basis re-evaluation | LRU cache keyed on (knots, order, t-hash) |
| Loops over curves | einsum-batched, one BLAS call |
| Cholesky per λ in GCV search | factorize once, rank-1 λ updates |
| Single-core | optional Numba `parallel=True`; torch → GPU |

**Targets:** ≥10× R on smoothing 1k curves; ≥50× on GCV grid search.

---

## 7. Parity & Testing

- **Golden files:** Docker + rpy2 runs R fda v6.3.0 on fixed inputs → JSON.
- **pytest:** every public symbol vs golden, `rtol=1e-8` (1e-5 for iterative: monotone, register).
- **Acceptance:** reproduce all 76 figures from the 2009 book.
- **CI gate:** parity job must pass to merge.
- **License note:** clean-room from books + docs. No R source translation → BSD-3 safe.

---

## 8. Datasets (`fabel.datasets`)

`load_canadian_weather()`, `load_growth()`, `load_gait()`, `load_handwriting()`,
`load_pinch()`, `load_melanoma()`, `load_refinery()`, `load_seabird()`,
`load_cstr()`, `load_chinese_script()`, `load_regina_precip()` — lazy-download, cached, returns `FData` + metadata.

---

## 9. Team Split & Phases

| Phase | Owner | Deliverable | Depends on |
|---|---|---|---|
| 0 (wk 1) | All | Freeze this spec | — |
| 1 | Eng 1 | `basis` + `core` + `_backend` | — |
| 1 | Eng 2 | Parity harness + golden files | — |
| 2 | Eng 3 | `smoothing` + `LDO` + `stats` | Phase 1 |
| 2 | Eng 4 | `datasets` + I/O | Phase 1 |
| 3 | Eng 3 | `decomposition` + `regression` | Phase 2 |
| 3 | Eng 4 | `nn` + sklearn compliance tests | Phase 2 |
| 4 | Eng 1+4 | `registration` + `dynamics` (hardest) | Phase 3 |
| 5 | All | Docs, 76-figure notebook, PyPI | Phase 4 |

---

## 10. Explicitly NOT in v1

- R plotting quirks (`axisIntervals`, `matplot`) — replaced by matplotlib/plotly defaults
- `readHMD` (Human Mortality DB scraper) — out of scope
- Deprecated/internal R helpers (`wtcheck`, `symsolve`, `polintmat`, `lnsrch`...) — private or scipy
