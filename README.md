# Fabel

**Repo type:** Shared library / package (PyPI) — AISMA Engineering Manual §8.5.

Fabel is a clean-room Python rewrite of R's `fda` package (Ramsay, Hooker & Graves) for functional data analysis:
basis expansions, penalised smoothing, functional PCA/CCA, functional regression, curve registration and
principal differential analysis. ~40 public symbols replace 515 R functions with identical numerics
(golden-file parity against R `fda` 6.3.0).

```python
import fabel as fb

w = fb.datasets.load_canadian_weather()
fd = fb.smooth(w.temp, w.t, basis=fb.Fourier(domain=(0, 365), n_basis=65))
fd.mean().plot()
pca = fb.FPCA(n=4).fit(fd)
```

See `SPEC.md` for the API surface and `docs/` for tutorials.
