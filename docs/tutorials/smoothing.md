# Smoothing: from raw numbers to curves

Functional data analysis starts with a step called *smoothing*. You have noisy
numbers measured at a few time points. You want a smooth curve that follows
them, and that you can evaluate and differentiate anywhere.

This tutorial uses two datasets:

- the **Berkeley growth study**: the height of 54 girls, measured 31 times
  between age 1 and 18 (ships inside fdatools);
- the **Canadian weather** data: the average temperature of every day of the
  year at 35 weather stations (downloaded once, then cached).

You will learn how `fdt.smooth()` picks a basis and a smoothing level for you,
how to control both, and how to force a curve to only go up.

## 1. The growth data

```python
import matplotlib.pyplot as plt
import numpy as np
import fdatools as fdt

growth = fdt.datasets.load_growth()
age = growth.age            # 31 ages, not equally spaced
heights = growth.hgtf       # shape (31, 54): one column per girl
age.shape, heights.shape
```

Each column of `heights` is one girl. In fdatools, observations always have the
time axis first: `y` has shape `(n_points, n_curves)`.

## 2. One call, everything automatic

```python
fit = fdt.smooth(heights, age)

fit.fd.basis                # the basis chosen for you
fit.lam                     # the smoothing parameter chosen by GCV
fit.df                      # the equivalent degrees of freedom
```

Without arguments `smooth()` does three things for you:

1. **Basis.** A *basis* is a fixed set of simple functions; each curve is a
   weighted sum of them. The default is a cubic *B-spline* (piecewise
   polynomials joined smoothly) with `min(n_points, 40) + 2` functions. Here
   that is 33 functions for 31 ages.
2. **Penalty.** With that many functions the curve could pass through every
   point, noise included. So the fit adds a *roughness penalty*: λ (lambda)
   times the integral of the squared second derivative. A large λ gives a
   stiff, smooth curve; a small λ follows the data closely.
3. **λ by GCV.** *Generalised cross-validation* (GCV) estimates how well the
   fit would predict new data. `smooth()` tries λ from 10⁻⁸ to 10⁸ and keeps
   the one with the lowest mean GCV score.

The *degrees of freedom* `fit.df` tell you how flexible the chosen fit is:
roughly, how many free parameters it really uses. For the growth data GCV
picks a very small λ and about 30 degrees of freedom for 31 points: the
heights were measured with little error, so the curve may follow them closely.

## 3. What is inside a `SmoothResult`

`smooth()` returns a `SmoothResult`. It holds everything R's `smooth.basis`
scatters over a list:

```python
fit.fd                      # FData: the 54 fitted curves
fit.gcv.shape               # one GCV score per girl: (54,)
fit.sse                     # residual sum of squares over all points and girls
fit.penalty_matrix.shape    # the roughness matrix R, (33, 33)
fit.y2c_map.shape           # the linear map from data to coefficients, (33, 31)
fit.constraint is None      # True: no shape constraint
```

`y2c_map` is the matrix that turns raw data into coefficients. You need it
later for standard errors in [regression](regression.md).

The curves are callable. Evaluate them, or their derivatives, at any point:

```python
t = np.linspace(1.0, 18.0, 200)
height = fit.fd(t)                 # (200, 54)
velocity = fit.fd(t, deriv=1)      # growth speed, cm per year
acceleration = fit.fd(t, deriv=2)  # change in growth speed

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(age, heights[:, :5], "o", markersize=3)
axes[0].plot(t, height[:, :5])
axes[0].set(xlabel="age (years)", ylabel="height (cm)", title="Five girls")
axes[1].plot(t, velocity[:, :5])
axes[1].set(xlabel="age (years)", ylabel="cm / year", title="Growth speed")
fig.tight_layout()
```

## 4. Choosing λ yourself

Three ways to set the smoothing level:

```python
basis = fdt.BSpline(domain=(1.0, 18.0), n_basis=20, order=6)

by_value = fdt.smooth(heights, age, basis=basis, lam=1e-1, penalty=3)
by_df = fdt.smooth(heights, age, basis=basis, df=8.0, penalty=3)
by_gcv = fdt.smooth(heights, age, basis=basis, lam="gcv", penalty=3)

round(by_df.df, 6), by_df.lam, by_gcv.lam
```

- `lam=1e-1` uses that number directly.
- `df=8.0` (or `lam="df=8"`) finds the λ that gives 8 degrees of freedom.
- `lam="gcv"` is the default.

We used an order-6 B-spline (piecewise polynomials of degree 5) and
`penalty=3`, the third derivative. A penalty on `D³` keeps the *second*
derivative smooth, which matters when you study acceleration.

`gcv_curve` shows the whole GCV curve, and `lambda_to_df` / `df_to_lambda`
convert between the two scales:

```python
from fdatools.smoothing import df_to_lambda, gcv_curve, lambda_to_df

lambdas = 10.0 ** np.arange(-6.0, 3.0, 0.5)
scores = gcv_curve(heights, age, basis, lambdas, penalty=3)   # (n_lambda, 54)

fig, ax = plt.subplots()
ax.loglog(lambdas, scores.mean(axis=1), marker="o")
ax.set(xlabel="lambda", ylabel="mean GCV")

lam_8 = df_to_lambda(age, basis, 8.0, penalty=3)
round(lambda_to_df(age, basis, lam_8, penalty=3), 6)          # 8.0
```

## 5. Heights never go down: the monotone constraint

A plain smooth can wiggle down a little between two measurements. Height
cannot. `constraint="monotone"` fits curves of the form

x(t) = β₀ + β₁ ∫ exp W(u) du,

where `W` is an unconstrained curve. Because `exp W` is always positive, x(t)
can only increase.

```python
girls = heights[:, :4]
monotone = fdt.smooth(
    girls, age, basis=fdt.BSpline(domain=(1.0, 18.0), n_basis=15), lam=1e-1,
    constraint="monotone",
)
monotone.constraint, monotone.beta.shape      # "monotone", (2, 4)

speed = monotone(t, deriv=1)                  # call the result, not .fd
bool(np.all(speed > 0))                       # True: growth speed is always positive
```

For a constrained fit, `monotone.fd` holds the latent curve `W`. Call the
result itself, `monotone(t)`, to get heights. `constraint="positive"` works
the same way and fits `x = exp W`, a curve that stays above zero.

## 6. Temperature: a Fourier basis and the harmonic accelerator

Weather repeats every year, so a *Fourier* basis (sines and cosines) is the
natural choice. We use 65 functions over the 365 days.

```python
# requires-data: canadian_weather
weather = fdt.datasets.load_canadian_weather()
day = weather.t                               # 1, 2, ..., 365
fourier = fdt.Fourier(domain=(0.0, 365.0), n_basis=65)
weather.temp.shape                            # (365, 35)
```

Now the penalty. The second derivative `D²` punishes every curve that bends.
But a pure yearly cycle `a + b sin(ωt) + c cos(ωt)` is exactly what we
expect, and we do not want to punish it. The *harmonic accelerator*

L = ω² D + D³, with ω = 2π / 365,

gives exactly zero on those curves. A *linear differential operator* (LDO) is
a weighted sum of derivatives like this one; `fdt.LDO.harmonic` builds it.

```python
# requires-data: canadian_weather
harmonic = fdt.LDO.harmonic(period=365.0)
temp = fdt.smooth(weather.temp, day, basis=fourier, penalty=harmonic, lam="gcv")
temp.lam, round(temp.df, 2)

stations = ["Montreal", "Edmonton", "Pr. Rupert", "Resolute"]
columns = [weather.stations.index(name) for name in stations]
fig, ax = plt.subplots()
ax.plot(day, weather.temp[:, columns], ".", markersize=2, alpha=0.4)
ax.set_prop_cycle(None)
ax.plot(day, temp.fd(day)[:, columns])
ax.legend(stations)
ax.set(xlabel="day of year", ylabel="deg C", title="Temperature, smoothed")
```

Check the claim: the harmonic accelerator kills a pure yearly cycle.

```python
# requires-data: canadian_weather
omega = 2 * np.pi / 365.0
cycle = fdt.FData(np.linalg.lstsq(fourier(day), 3 + 5 * np.sin(omega * day), rcond=None)[0], fourier)
float(np.max(np.abs(cycle(day, harmonic))))   # about 0
```

## R equivalent

```r
library(fda)
# growth, automatic GCV is a loop in R: smooth.basis at each lambda + gcv
hgtbasis <- create.bspline.basis(c(1, 18), 20, norder = 6)
hgtPar   <- fdPar(hgtbasis, Lfdobj = 3, lambda = 1e-1)
hgtfit   <- smooth.basis(growth$age, growth$hgtf, hgtPar)
hgtfit$df; hgtfit$gcv; hgtfit$SSE; hgtfit$y2cMap
lambda   <- df2lambda(growth$age, hgtbasis, Lfdobj = 3, df = 8)

# monotone
Wfd  <- fdPar(create.bspline.basis(c(1, 18), 15), 2, 1e-1)
mono <- smooth.monotone(growth$age, growth$hgtf[, 1], Wfd)

# weather, harmonic accelerator
daybasis <- create.fourier.basis(c(0, 365), 65)
harmLfd  <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
tempfit  <- smooth.basis(day.5, CanadianWeather$dailyAv[, , "Temperature.C"],
                         fdPar(daybasis, harmLfd, 1e2))
```

Differences to know:

- `fdt.smooth()` picks the basis and λ for you; R needs both.
- `fit.gcv` is a score per curve, like R. When the fit interpolates the data
  (df equal to the number of points), GCV is 0/0: fdatools returns `inf`, R
  returns `NULL`.
- One `smooth()` call with `constraint=` replaces `smooth.monotone`,
  `smooth.pos` and `smooth.morph`, and it smooths all curves at once.
- R's day grid is `day.5` (0.5, 1.5, ...); `load_canadian_weather().t` is
  1, ..., 365. Both lie inside the domain (0, 365).
