# Functional regression: rain, temperature and region

*Functional regression* is linear regression where the response, a
covariate, or both are curves. One function, `fregress()`, covers every case.
It reads the kind of model from what you pass in:

| Response `y` | Covariate `x` | Model |
|---|---|---|
| numbers | numbers | ordinary regression |
| numbers | curves | y_i = α + ∫ x_i(t) β(t) dt (*scalar-on-function*) |
| curves | numbers | y_i(t) = β₀(t) + z_i β₁(t) (*function-on-scalar*) |
| curves | curves | y_i(t) = β₀(t) + x_i(t) β₁(t) (*concurrent*) |

The coefficients β are curves, so they need a basis and, usually, a roughness
penalty λ, as in smoothing.

This tutorial uses the Canadian weather data (downloaded once, then cached).

## 1. Smooth the temperatures

```python
# requires-data: canadian_weather
import matplotlib.pyplot as plt
import numpy as np
import fdatools as fdt
from fdatools.regression import fregress
from fdatools.stats import f_test

weather = fdt.datasets.load_canadian_weather()
day = weather.t
fourier = fdt.Fourier(domain=(0.0, 365.0), n_basis=65)
harmonic = fdt.LDO.harmonic(period=365.0)

temp_fit = fdt.smooth(weather.temp, day, basis=fourier, penalty=harmonic, lam=1e2)
temp = temp_fit.fd                                    # 35 temperature curves
```

## 2. Scalar response: total rain from the temperature curve

Can the temperature profile over the year predict how much it rains in a
year? The response is one number per station, the log of the total annual
precipitation:

```python
# requires-data: canadian_weather
log_precip = np.log10(weather.precip.sum(axis=0))     # 35 numbers
log_precip.shape
```

The model is

log_precip_i = α + ∫ temp_i(t) β(t) dt + error.

With 35 stations we cannot estimate a free curve β(t) with 65 coefficients,
so we give β a 35-function Fourier basis and a strong harmonic-accelerator
penalty. A covariate given as a number (`1.0`) is the intercept.

```python
# requires-data: canadian_weather
beta_basis = fdt.Fourier(domain=(0.0, 365.0), n_basis=35)
model = fregress(
    log_precip,
    {"const": 1.0, "temp": temp},
    beta={"temp": (beta_basis, 10**12.5, harmonic)},  # (basis, lambda, penalty)
)

model.names                          # ('const', 'temp')
round(model.df, 2)                   # equivalent degrees of freedom of the fit
model.gcv, model.ocv                 # GCV and leave-one-out CV scores

r_squared = 1 - np.sum((log_precip - model.fitted) ** 2) / np.sum(
    (log_precip - log_precip.mean()) ** 2
)
round(r_squared, 2)                  # about 0.75
```

`model.beta` holds one `FData` per term. The intercept lives in a `Constant`
basis; the temperature coefficient is a curve:

```python
# requires-data: canadian_weather
alpha = model.beta[0].coefs[0, 0]
beta_t = model.beta[1]

fig, ax = plt.subplots()
ax.plot(day, beta_t(day))
ax.axhline(0.0, color="grey", linewidth=0.5)
ax.set(xlabel="day of year", ylabel="beta(t)", title="Effect of temperature on log rain")
```

To choose λ, compare `model.gcv` or `model.cv().sse` for a few values and keep
the smallest.

## 3. Predictions, standard errors and cross-validation

`predict` takes new covariates laid out like the ones you fitted with.
`stderr()` gives the standard error of every coefficient, and `cv()` refits
the model leaving out one station at a time.

```python
# requires-data: canadian_weather
model.predict({"const": 1.0, "temp": temp[0:3]})      # first three stations

se = model.stderr()                  # residual variance estimated from the fit
se.cov.shape                         # covariance of all coefficients, (36, 36)
band = 2 * se.beta[1](day)[:, 0]     # pointwise +- 2 standard errors of beta(t)

fig, ax = plt.subplots()
ax.plot(day, beta_t(day), color="black")
ax.plot(day, beta_t(day)[:, 0] + band, "--", color="grey")
ax.plot(day, beta_t(day)[:, 0] - band, "--", color="grey")
ax.set(xlabel="day of year", title="beta(t) with 2 standard errors")

cv = model.cv()
cv.errors.shape, round(cv.sse, 3)    # leave-one-out errors and their sum of squares
```

The standard errors ignore the uncertainty in choosing λ, as in R. Read them
as a rough guide, not an exact confidence band.

## 4. Functional response: temperature by climate region

Now the response is the temperature curve itself and the covariate is the
station's climate region: Arctic, Atlantic, Continental or Pacific.

The formula interface builds the design for you. A covariate given as text
labels becomes *indicator* (0/1) variables, one per region, with the first
region in alphabetical order (Arctic) left out as the reference. `const` is
the intercept.

```python
# requires-data: canadian_weather
by_region = fregress("temp ~ region", {"temp": temp_fit, "region": weather.region})
by_region.names       # ('const', 'region.Atlantic', 'region.Continental', 'region.Pacific')

fig, ax = plt.subplots()
for name, beta in zip(by_region.names, by_region.beta):
    ax.plot(day, beta(day), label=name)
ax.legend()
ax.set(xlabel="day of year", ylabel="deg C", title="Region effects")
```

`const` is the mean Arctic curve. Each other curve is how much warmer that
region is than the Arctic on each day: much warmer in winter than in summer.

We passed `temp_fit` (the `SmoothResult`), not only `temp_fit.fd`. That keeps
the smoothing map `y2c_map`, which `stderr()` needs for a functional response.
It also needs the covariance of the residuals across days, `sigma_e`:

```python
# requires-data: canadian_weather
residuals = weather.temp - by_region.fitted(day)      # (365, 35)
sigma_e = residuals @ residuals.T / residuals.shape[1]
region_se = by_region.stderr(sigma_e=sigma_e)
len(region_se.beta)                                   # one standard-error curve per term

region_cv = by_region.cv()
round(region_cv.sse)                                  # summed integrated squared errors
```

## 5. Is the effect real? A permutation F test

A *permutation test* asks: if the covariates had nothing to do with the
response, how large would the F statistic be by chance? It shuffles the
responses many times, refits each time, and compares the real F with the
shuffled ones.

`fdt.stats.f_test` takes the raw regression inputs, like R's `Fperm.fd`. It
does not add an intercept: pass a column of ones yourself.

```python
# requires-data: canadian_weather
regions = sorted(set(weather.region))                 # Arctic first
dummies = [
    np.array([1.0 if r == name else 0.0 for r in weather.region]) for name in regions[1:]
]
test = f_test(
    temp,
    [np.ones(35), *dummies],
    basis=fdt.Fourier(domain=(0.0, 365.0), n_basis=11),
    n_perm=100,
    random_state=1,
)
round(test.statistic, 2), test.pvalue, round(test.critical_value, 2)

fig, ax = plt.subplots()
ax.plot(test.t, test.pointwise, label="observed F(t)")
ax.plot(test.t, test.pointwise_critical_value, "--", label="pointwise 95% of null")
ax.axhline(test.critical_value, linestyle=":", label="95% of null maximum")
ax.legend()
ax.set(xlabel="day of year", title="Permutation F test")
```

None of the 100 shuffles came close: `pvalue` is 0, so region clearly matters.
Use more permutations (the default is 200) for a finer
p-value.

The same function tests a scalar response. Then F is a single number:

```python
# requires-data: canadian_weather
scalar_test = f_test(
    log_precip,
    [np.ones(35), temp],
    basis=[None, beta_basis],
    lam=[0.0, 10**12.5],
    penalty=[2, harmonic],
    n_perm=100,
    random_state=1,
)
scalar_test.pvalue                                    # 0.0
```

## R equivalent

```r
library(fda)
daybasis <- create.fourier.basis(c(0, 365), 65)
harmLfd  <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
tempSmooth <- smooth.basis(day.5, CanadianWeather$dailyAv[, , "Temperature.C"],
                           fdPar(daybasis, harmLfd, 1e2))
tempfd   <- tempSmooth$fd
annprec  <- log10(apply(CanadianWeather$dailyAv[, , "Precipitation.mm"], 2, sum))

# scalar response
xfdlist  <- list(const = rep(1, 35), temp = tempfd)
betalist <- list(const = fdPar(create.constant.basis(c(0, 365))),
                 temp  = fdPar(create.fourier.basis(c(0, 365), 35), harmLfd, 10^12.5))
fit      <- fRegress(annprec, xfdlist, betalist)
fit$df; fit$gcv; fit$OCV
cv       <- fRegress.CV(annprec, xfdlist, betalist)

# functional response with a formula
region   <- factor(CanadianWeather$region)
regfit   <- fRegress(tempfd ~ region)

# permutation F test
regionX  <- list(const = rep(1, 35),
                 atl = as.numeric(region == "Atlantic"),
                 con = as.numeric(region == "Continental"),
                 pac = as.numeric(region == "Pacific"))
betaR    <- rep(list(fdPar(create.fourier.basis(c(0, 365), 11))), 4)
Fres     <- Fperm.fd(tempfd, regionX, betaR, nperm = 100)
```

Differences to know:

- `fregress()` replaces `fRegress`, `fRegress.formula`, `fRegress.fd` and
  `fRegress.double`. `.predict()`, `.stderr()` and `.cv()` on the result
  replace `predict.fRegress`, `fRegress.stderr` and `fRegress.CV`.
- fdatools computes every design integral exactly. R approximates them, so
  coefficients can differ from R's around the sixth digit.
- `ocv` is the *sum* of squared leave-one-out residuals, as in R, while `gcv`
  is SSE / (n - df)².
- R's `predict` fails on a functional-response model; fdatools' works.
- `f_test` takes the raw inputs, like `Fperm.fd`. It does not yet accept a
  fitted `fregress` result.
