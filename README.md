# Fabel

[![CI](https://github.com/AISMAsrl/fabel/actions/workflows/ci.yml/badge.svg)](https://github.com/AISMAsrl/fabel/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/fabel.svg)](https://pypi.org/project/fabel/)
[![Python versions](https://img.shields.io/pypi/pyversions/fabel.svg)](https://pypi.org/project/fabel/)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/AISMAsrl/fabel/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-aismasrl.github.io%2Ffabel-blue.svg)](https://aismasrl.github.io/fabel)

**Functional data analysis for Python.** Fabel is a clean-room Python rewrite of
R's [`fda`](https://cran.r-project.org/package=fda) package (Ramsay, Hooker &
Graves): basis expansions, penalised smoothing, functional PCA and CCA,
functional regression, curve registration and principal differential analysis.

- **Same numbers as R.** Every public function is tested against golden output
  from R `fda` 6.3.0 at `rtol = 1e-8` (`1e-5` for iterative fits). See the
  [parity report](https://github.com/AISMAsrl/fabel/blob/main/PARITY_REPORT.md).
- **About 40 symbols instead of 515 functions.** Curves are callable (`fd(t)`),
  arithmetic is plain Python (`fd1 + fd2`, `fd1 @ fd2`), and there is one
  `smooth()`, one `register()`, one `.plot()`.
- **Fits the Python stack.** scikit-learn estimators, NumPy or PyTorch arrays
  through the Array API (gradients flow), pandas and xarray I/O.

## Install

```bash
pip install fabel
```

Optional extras: `fabel[plot]` (matplotlib), `fabel[pandas]` (pandas + xarray),
`fabel[io]` (read R `.rds` files), `fabel[torch]` (PyTorch). Python 3.10 to 3.13.

## Quickstart

```python
import numpy as np
import fabel as fb

growth = fb.datasets.load_growth()                   # Berkeley growth study, ships with Fabel
basis = fb.BSpline(domain=(1.0, 18.0), n_basis=12)   # cubic B-splines on ages 1-18
coefs, *_ = np.linalg.lstsq(basis(growth.age), growth.hgtf, rcond=None)
girls = fb.FData(coefs, basis)                       # 54 height curves in one object
speed = girls.derivative()                           # growth speed in cm/year, exact
print(speed.mean()(np.array([5.0, 12.0])))           # curves are callable: evaluate anywhere
print(fb.inprod(girls[0], girls[0]))                 # L2 inner product of one curve
```

Next steps: penalised smoothing with automatic λ (`smooth`), functional PCA
(`FPCA`) and more in the [quickstart guide](https://aismasrl.github.io/fabel/quickstart/).
Coming from R? The [migration table](https://aismasrl.github.io/fabel/r-migration/)
maps every `fda` function to its Fabel equivalent.

## Documentation

Full documentation, API reference and tutorials: <https://aismasrl.github.io/fabel>.

## Citing

If you use Fabel in published work, please cite it together with Ramsay &
Silverman (2005); see [`CITATION.cff`](https://github.com/AISMAsrl/fabel/blob/main/CITATION.cff).

## License

BSD 3-Clause. Fabel is a clean-room implementation written from the published
mathematics (Ramsay & Silverman 2005; Ramsay, Hooker & Graves 2009) and the
public behaviour of R `fda`; it contains no R source code.
