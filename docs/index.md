# fdatools

**Functional data analysis for Python.** fdatools is a clean-room Python rewrite of
R's `fda` package (Ramsay, Hooker & Graves, version 6.3.0). It covers basis
expansions, penalised smoothing, functional PCA and CCA, functional regression,
curve registration and principal differential analysis.

## Why fdatools

- **Same numbers as R.** Every public function is tested against golden output
  from R `fda` 6.3.0 at `rtol = 1e-8` (`1e-5` for iterative fits such as
  `smooth.monotone` and `register.fd`). Where R is demonstrably less accurate,
  the difference is measured and documented, never hidden.
- **Small surface.** About 40 public symbols replace 515 R functions. A curve
  is callable (`fd(t)`), arithmetic is plain Python (`fd1 + fd2`, `fd1 @ fd2`
  for the inner product), and every method returns a new object.
- **Fits the Python stack.** Estimators follow the scikit-learn API, so they go
  straight into a `Pipeline` or `GridSearchCV`. Arrays go through the Array API,
  so NumPy in gives NumPy out and PyTorch in gives PyTorch out, with gradients.

## Install

```bash
uv add "fdatools @ git+https://github.com/hameddavodi/fdatools@v1.1.0"
```

fdatools is installed from GitHub with [uv](https://docs.astral.sh/uv/). With an
extra, put it in brackets: `uv add "fdatools[plot] @ git+https://github.com/hameddavodi/fdatools@v1.1.0"`.

| Extra | Adds |
|---|---|
| `fdatools[plot]` | matplotlib, for `.plot()` |
| `fdatools[pandas]` | pandas and xarray, for `from_pandas` / `to_pandas` / `to_xarray` |
| `fdatools[io]` | `rdata`, for `read_rds` (load R `fd` objects) |
| `fdatools[torch]` | PyTorch backend and `fdatools.nn` |

fdatools supports Python 3.10 to 3.13.

## A first look

```python
import numpy as np
import fdatools as fdt

growth = fdt.datasets.load_growth()
basis = fdt.BSpline(domain=(1.0, 18.0), n_basis=12)
coefs, *_ = np.linalg.lstsq(basis(growth.age), growth.hgtf, rcond=None)
girls = fdt.FData(coefs, basis)
girls.derivative().mean()(np.array([5.0, 12.0]))   # mean growth speed at 5 and 12
```

## Where to go next

- [Quickstart](quickstart.md): the core objects and a full smoothing + FPCA workflow.
- [Migrating from R `fda`](r-migration.md): every R function and its fdatools equivalent.
- [Tutorials](tutorials/index.md): six end-to-end analyses (smoothing, FPCA, registration,
  regression, dynamics, machine learning).
- [API reference](api/index.md): every public class and function.

## References

- J. O. Ramsay and B. W. Silverman, *Functional Data Analysis*, 2nd ed., Springer, 2005.
- J. O. Ramsay, G. Hooker and S. Graves, *Functional Data Analysis with R and MATLAB*, Springer, 2009.
