# Functional PCA and CCA: the main shapes in the weather

You have 35 temperature curves, one per weather station. How do they differ?
*Functional principal component analysis* (FPCA) answers this. It finds a few
typical shapes, called *harmonics*, so that every station's curve is close to

mean curve + score₁ × harmonic₁ + score₂ × harmonic₂ + ...

Each *score* says how much of a harmonic one station has. The *share of
variance* of a harmonic says how much of the differences between stations it
explains.

At the end we use *functional canonical correlation analysis* (FCCA) to ask a
second question: which temperature patterns go together with which
precipitation patterns?

This tutorial needs the Canadian weather data. It is downloaded on first use
and cached (see [Datasets](../api/datasets.md)).

## 1. Smooth the data

We smooth temperature and log precipitation with a 65-function Fourier basis
and the harmonic accelerator penalty, as in the
[smoothing tutorial](smoothing.md).

```python
# requires-data: canadian_weather
import matplotlib.pyplot as plt
import numpy as np
import fabel as fb
from fabel.decomposition import FCCA, FPCA

weather = fb.datasets.load_canadian_weather()
day = weather.t
fourier = fb.Fourier(domain=(0.0, 365.0), n_basis=65)
harmonic = fb.LDO.harmonic(period=365.0)

temp = fb.smooth(weather.temp, day, basis=fourier, penalty=harmonic, lam=1e2).fd
precip = fb.smooth(weather.log10precip, day, basis=fourier, penalty=harmonic, lam=1e4).fd
temp.n_curves, precip.n_curves          # 35 stations each
```

## 2. Fit FPCA

`FPCA` is a scikit-learn *estimator*: you create it with its settings, then
call `.fit()`.

```python
# requires-data: canadian_weather
pca = FPCA(n=4).fit(temp)

pca.varprop.round(4)          # share of variance of each harmonic
pca.harmonics.n_curves        # 4 harmonic curves
pca.scores.shape              # (35, 4): one row per station
pca.mean_fd.n_curves          # the mean curve
```

The first harmonic explains about 89% of the variation, the second about 8%.
Four harmonics together explain more than 99%.

`pca.values` holds the *full* list of eigenvalues (one per basis function),
as R's `pca.fd` reports it. Only the first `n` harmonics are kept.

## 3. Look at the harmonics

A harmonic is easiest to read as a change of the mean curve: the mean plus and
minus a multiple of the harmonic. `pca.plot()` draws that, one panel per
harmonic.

```python
# requires-data: canadian_weather
pca.plot()

fig, ax = plt.subplots()
ax.plot(day, pca.harmonics(day))
ax.legend([f"harmonic {k + 1}" for k in range(4)])
ax.set(xlabel="day of year", title="Temperature harmonics")
```

- **Harmonic 1** is positive all year and largest in winter. A station with a
  high score is warmer overall, and much warmer in winter: a mild, coastal
  climate.
- **Harmonic 2** has opposite signs in summer and winter. It measures the
  size of the seasonal swing: continental stations have hot summers and cold
  winters.

## 4. Scores: where each station sits

The scores are ordinary numbers, one row per station, ready for a scatter plot
or any other statistical tool.

```python
# requires-data: canadian_weather
fig, ax = plt.subplots()
for region in sorted(set(weather.region)):
    rows = [i for i, r in enumerate(weather.region) if r == region]
    ax.scatter(pca.scores[rows, 0], pca.scores[rows, 1], label=region)
ax.legend()
ax.set(xlabel="score 1 (warm overall)", ylabel="score 2 (seasonal swing)")

# transform() gives the same scores, and works on new curves too
bool(np.allclose(pca.transform(temp), pca.scores))
```

`inverse_transform` goes the other way, from scores back to curves. With four
harmonics the rebuilt curves are very close to the originals:

```python
# requires-data: canadian_weather
rebuilt = pca.inverse_transform(pca.scores)
error = float(np.max(np.abs(rebuilt(day) - temp(day))))
round(error, 1)               # a few degrees at most, on curves that span 50 degrees
```

## 5. A smoother set of harmonics

Harmonics can be noisy when the curves are. Pass `lam` (and `penalty`) to
penalise their roughness, just like in smoothing. `lam="gcv"` chooses it by
leave-one-curve-out cross-validation, which is slower.

```python
# requires-data: canadian_weather
smooth_pca = FPCA(n=4, lam=1e4, penalty=harmonic).fit(temp)
smooth_pca.varprop.round(4)
```

## 6. Varimax rotation

Harmonics are built one after the other, each explaining as much variance as
possible. That often makes the later ones hard to read. A *varimax rotation*
mixes the kept harmonics so that each one is large on a short part of the
year and near zero elsewhere. The total variance explained stays the same.

```python
# requires-data: canadian_weather
rotated = pca.rotate("varimax")          # a new fitted FPCA; pca is unchanged
rotated.varprop.round(4)                 # now spread more evenly
bool(np.isclose(rotated.varprop.sum(), pca.varprop.sum()))   # True
rotated.rotation.shape                   # (4, 4) orthogonal matrix

fig, ax = plt.subplots()
ax.plot(day, rotated.harmonics(day))
ax.set(xlabel="day of year", title="Rotated harmonics")
```

After rotation the harmonics read like seasons: each one is largest in one
season (winter, spring, autumn or summer) and small in the others.

## 7. FCCA: temperature against precipitation

*Canonical correlation* looks for a weight function for temperature and one
for precipitation such that the two resulting station scores are as
correlated as possible. Without a penalty this is meaningless with 35
stations and 65 basis functions (every correlation would be 1), so we
penalise the weights with `lam1` and `lam2`.

```python
# requires-data: canadian_weather
cca = FCCA(n=3, lam1=1e6, lam2=1e6, penalty=harmonic).fit(temp, precip)

cca.correlations[:3].round(3)            # canonical correlations, largest first
cca.scores1.shape, cca.scores2.shape     # (35, 3) each

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(day, cca.weights1(day)[:, 0], label="temperature")
axes[0].plot(day, cca.weights2(day)[:, 0], label="log precipitation")
axes[0].legend()
axes[0].set(xlabel="day of year", title="First pair of weight functions")
axes[1].scatter(cca.scores1[:, 0], cca.scores2[:, 0])
axes[1].set(xlabel="temperature score", ylabel="precipitation score")
fig.tight_layout()
```

The first pair has a correlation near 0.98. The temperature weight is largest
in winter, so the first temperature score is high for stations with mild
winters (St. Johns, Victoria) and low for cold inland and Arctic stations
(Dawson, Resolute). The coastal stations are also the wet ones. Canonical correlations are easy to
over-read with only 35 stations, so treat them as a description of this
sample, not a test.

## R equivalent

```r
library(fda)
daybasis <- create.fourier.basis(c(0, 365), 65)
harmLfd  <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
tempfd   <- smooth.basis(day.5, CanadianWeather$dailyAv[, , "Temperature.C"],
                         fdPar(daybasis, harmLfd, 1e2))$fd
precfd   <- smooth.basis(day.5, CanadianWeather$dailyAv[, , "log10precip"],
                         fdPar(daybasis, harmLfd, 1e4))$fd

temppca  <- pca.fd(tempfd, nharm = 4)
temppca$varprop; temppca$scores; plot(temppca)
rotpca   <- varmx.pca.fd(temppca)

ccafd    <- cca.fd(tempfd, precfd, ncan = 3,
                   fdPar(daybasis, harmLfd, 1e6), fdPar(daybasis, harmLfd, 1e6))
ccafd$ccacorr
```

Differences to know:

- Fabel computes the scores and the Gram matrix with exact integrals. R uses
  numerical integration that is good to about 4 or 5 digits, so R's scores can
  differ from Fabel's in the third or fourth digit.
- R's `varmx` stops a little before the varimax optimum. Fabel iterates to the
  optimum, so rotated harmonics can differ slightly.
- Signs: every harmonic has a positive coefficient sum, the rule R follows.
  For FCCA the two members of a pair share one sign so their scores stay
  positively correlated (R has no rule).
- `FPCA(lam="gcv")` is a Fabel addition; R's `pca.fd` has no automatic
  choice.
