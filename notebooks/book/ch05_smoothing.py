# %% [markdown]
# # Chapter 5: Smoothing: Computing Curves from Noisy Data
#
# Figures 5.1 to 5.8 of Ramsay, Hooker and Graves (2009): choosing the smoothing
# parameter by GCV, roughness-penalised smoothing of the growth and log
# precipitation data, and constrained fits (positive, monotone, density).

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import fdatools as fdt

plt.rcParams["figure.max_open_warning"] = 0

growth = fdt.datasets.load_growth()
weather = fdt.datasets.load_canadian_weather()
age = growth.age
hgtf = growth.hgtf
day5 = weather.t - 0.5  # the book's day.5: mid-day points 0.5, ..., 364.5
harmonic = fdt.LDO.harmonic(period=365.0)
daybasis = fdt.Fourier(domain=(0.0, 365.0), n_basis=365)

# %% [markdown]
# ### Figure 5.1
# GCV criterion against log10(lambda) for smoothing the 54 Berkeley girls' heights
# with order-6 B-splines (knots at the ages) and a fourth-derivative penalty. The
# criterion is summed over girls.

# %%
growthbasis = fdt.BSpline(domain=(1.0, 18.0), order=6, breaks=age)
loglam = np.arange(-6.0, 0.01, 0.25)
gcvsave = np.array(
    [
        float(np.sum(fdt.smooth(hgtf, age, basis=growthbasis, lam=10.0**ll, penalty=4).gcv))
        for ll in loglam
    ]
)
best_loglam = float(loglam[np.argmin(gcvsave)])
print(f"GCV minimum at log10(lambda) = {best_loglam:.2f}")

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(loglam, gcvsave, "o-", color="black", markersize=4)
ax.set_xlabel(r"$\log_{10}(\lambda)$")
ax.set_ylabel(r"GCV($\lambda$)")
ax.set_title("GCV for the girls' growth data")
fig

# %% [markdown]
# ### Figure 5.2
# Acceleration (second derivative) curves of the first ten girls with lambda = 0.1
# (the GCV minimum of Figure 5.1 fits the heights well but leaves the second
# derivative rough); the heavy dashed line is the mean acceleration of all 54 girls.

# %%
growthfit = fdt.smooth(hgtf, age, basis=growthbasis, lam=0.1, penalty=4)
print(f"growth: lambda = 0.1, df = {growthfit.df:.2f}")
agefine = np.linspace(1.0, 18.0, 401)
accel = np.asarray(growthfit.fd(agefine, deriv=2))

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(agefine, accel[:, :10], linewidth=1)
ax.plot(agefine, accel.mean(axis=1), "k--", linewidth=3, label="mean (54 girls)")
ax.axhline(0.0, color="grey", linestyle=":")
ax.set_ylim(-4.0, 2.0)
ax.set_xlabel("Age (years)")
ax.set_ylabel(r"Height acceleration (cm/yr$^2$)")
ax.legend()
fig

# %% [markdown]
# ### Figure 5.3
# Log10 daily precipitation at the 35 Canadian weather stations, smoothed with a
# 365-term Fourier basis and the harmonic acceleration penalty; lambda is chosen by
# GCV over 1e4, ..., 1e9 (the book's choice, 1e6, is the GCV minimum).

# %%
logprec = weather.log10precip
precgcv = np.array(
    [
        float(np.sum(fdt.smooth(logprec, day5, basis=daybasis, lam=10.0**ll, penalty=harmonic).gcv))
        for ll in range(4, 10)
    ]
)
prec_lam = 10.0 ** (4 + int(np.argmin(precgcv)))
logprecfit = fdt.smooth(logprec, day5, basis=daybasis, lam=prec_lam, penalty=harmonic)
print(f"log precipitation: lambda = {prec_lam:.0e}, df = {logprecfit.df:.2f}")

fig, ax = plt.subplots(figsize=(7, 4))
logprecfit.fd.plot(ax=ax, linewidth=1)
ax.set_xlabel("Day (July 1 to June 30)")
ax.set_ylabel(r"$\log_{10}$ precipitation (mm)")
ax.set_title("Smoothed log precipitation, 35 stations")
fig

# %% [markdown]
# ### Figure 5.4
# Vancouver's average daily precipitation with a positive smooth x(t) = exp W(t),
# W expanded in the 365-term Fourier basis with a harmonic acceleration penalty.

# %%
vancouver = weather.stations.index("Vancouver")
vanprec = weather.precip[:, vancouver]
vanfit = fdt.smooth(vanprec, day5, basis=daybasis, lam=1e4, penalty=harmonic, constraint="positive")
print(f"Vancouver positive smooth: lambda = 1e4, df = {vanfit.df:.2f}")
dayfine = np.linspace(0.0, 365.0, 731)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(day5, vanprec, "o", markersize=2, fillstyle="none", color="grey", label="data")
ax.plot(dayfine, np.asarray(vanfit(dayfine)), "k-", linewidth=2, label="positive smooth")
ax.set_xlabel("Day (July 1 to June 30)")
ax.set_ylabel("Precipitation (mm)")
ax.set_title("Vancouver precipitation")
ax.legend()
fig

# %% [markdown]
# ### Figure 5.5
# Tibia lengths of one newborn infant over its first 40 days, with a monotone smooth
# x(t) = b0 + b1 int exp W (order-6 B-splines, knots at every day, D^3 penalty on W).

# %%
infant = fdt.datasets.load_infant_growth()
tibiabasis = fdt.BSpline(domain=(1.0, 40.0), order=6, breaks=infant.day)
tibiafit = fdt.smooth(
    infant.tibia_length,
    infant.day,
    basis=tibiabasis,
    lam=1e-1,
    penalty=3,
    constraint="monotone",
)
dayinf = np.linspace(1.0, 40.0, 391)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(infant.day, infant.tibia_length, "o", fillstyle="none", color="grey", label="data")
ax.plot(dayinf, np.asarray(tibiafit(dayinf)), "k-", linewidth=2, label="monotone smooth")
ax.set_xlabel("Day")
ax.set_ylabel("Tibia length (mm)")
ax.legend()
fig

# %% [markdown]
# ### Figure 5.6
# Velocity of tibia growth from the monotone fit of Figure 5.5; the fit's derivative
# b1 exp W(t) is positive everywhere, showing the growth spurts.

# %%
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(dayinf, np.asarray(tibiafit(dayinf, deriv=1)), "k-", linewidth=2)
ax.axhline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Day")
ax.set_ylabel("Tibia velocity (mm/day)")
fig

# %% [markdown]
# ### Figure 5.7
# Height velocity of the first ten Berkeley girls from monotone smoothing (knots at
# the ages, order 6, D^3 penalty on W): unlike an unconstrained fit, the fitted
# heights can never decrease.

# %%
girlsfit = fdt.smooth(
    hgtf[:, :10], age, basis=growthbasis, lam=1e-1, penalty=3, constraint="monotone"
)
velocity = np.asarray(girlsfit(agefine, deriv=1))

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(agefine, velocity, linewidth=1)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Height velocity (cm/yr)")
ax.set_title("Monotone smooth, first ten girls")
fig

# %% [markdown]
# ### Figure 5.8
# Estimated density of June daily precipitation in Regina (days with 2 to 45 mm),
# p(x) = exp W(x) / int exp W. fdatools' smooth() has no density constraint, so the
# penalised log-likelihood is maximised here from the public Basis primitives.

# %%
regina = fdt.datasets.load_regina_precip().value
rain = np.sort(regina[(regina > 2.0) & (regina <= 45.0)])
densbasis = fdt.BSpline(domain=(2.0, 45.0), n_basis=13)
densgrid = np.linspace(2.0, 45.0, 861)
phi_data = np.asarray(densbasis(rain))
phi_grid = np.asarray(densbasis(densgrid))
penmat = np.asarray(densbasis.penalty(2))
quad = np.full(densgrid.size, densgrid[1] - densgrid[0])
quad[[0, -1]] *= 0.5
dens_lam = 1e-1


def neg_loglik(coefs):
    w = phi_grid @ coefs
    shift = w.max()
    mass = np.exp(w - shift)
    total = quad @ mass
    value = -phi_data.sum(axis=0) @ coefs + rain.size * (np.log(total) + shift)
    value += dens_lam * coefs @ penmat @ coefs
    grad = -phi_data.sum(axis=0) + rain.size * (phi_grid.T @ (quad * mass)) / total
    return value, grad + 2.0 * dens_lam * penmat @ coefs


opt = minimize(neg_loglik, np.zeros(densbasis.n_basis), jac=True, method="BFGS")
wgrid = phi_grid @ opt.x
density = np.exp(wgrid - wgrid.max())
density /= quad @ density

fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(rain, bins=30, density=True, color="lightgrey", edgecolor="grey")
ax.plot(densgrid, density, "k-", linewidth=2)
ax.plot(rain, np.zeros_like(rain), "|", color="black", markersize=8)
ax.set_xlabel("Precipitation (mm)")
ax.set_ylabel("Probability density")
ax.set_title(f"Regina June precipitation ({rain.size} days)")
fig
