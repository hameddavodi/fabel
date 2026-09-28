# %% [markdown]
# # Chapter 6: Descriptions of Functional Data
#
# Figures 6.1 to 6.8 of Ramsay, Hooker and Graves (2009): covariance and
# correlation surfaces, the residual variance, phase-plane plots (a sinusoid, the
# nondurable goods index, girls' growth) and pointwise confidence limits.

# %%
import matplotlib.pyplot as plt
import numpy as np

import fdatools as fdt

plt.rcParams["figure.max_open_warning"] = 0

weather = fdt.datasets.load_canadian_weather()
day5 = weather.t - 0.5  # the book's day.5: mid-day points 0.5, ..., 364.5
harmonic = fdt.LDO.harmonic(period=365.0)
daybasis = fdt.Fourier(domain=(0.0, 365.0), n_basis=365)
logprec = weather.log10precip
logprecfit = fdt.smooth(logprec, day5, basis=daybasis, lam=1e6, penalty=harmonic)
logprecfd = logprecfit.fd
tempfd = fdt.smooth(weather.temp, day5, basis=daybasis, lam=1e4, penalty=harmonic).fd
dayfine = np.linspace(0.0, 365.0, 121)

# %% [markdown]
# ### Figure 6.1
# The variance-covariance surface v(s, t) of the smoothed log precipitation curves,
# as a perspective plot (left) and a contour plot (right).

# %%
logprecvar = fdt.stats.cov(logprecfd)
surface = np.asarray(logprecvar(dayfine, dayfine))
grid_s, grid_t = np.meshgrid(dayfine, dayfine, indexing="ij")

fig = plt.figure(figsize=(11, 4.5))
ax3d = fig.add_subplot(1, 2, 1, projection="3d")
ax3d.plot_surface(grid_s, grid_t, surface, cmap="viridis", linewidth=0)
ax3d.set_xlabel("Day (s)")
ax3d.set_ylabel("Day (t)")
ax3d.set_zlabel("v(s, t)")
ax2 = fig.add_subplot(1, 2, 2)
contours = ax2.contour(dayfine, dayfine, surface.T, levels=12, cmap="viridis")
ax2.clabel(contours, fontsize=7)
ax2.set_xlabel("Day (s)")
ax2.set_ylabel("Day (t)")
ax2.set_aspect("equal")
fig.suptitle("Variance-covariance surface of log precipitation")
fig

# %% [markdown]
# ### Figure 6.2
# Cross-correlation corr(temperature(s), log precipitation(t)) across the 35
# stations, from smoothed temperature and log precipitation curves.

# %%
crosscor = np.asarray(fdt.stats.cor(tempfd, logprecfd, s=dayfine, t=dayfine))

fig, ax = plt.subplots(figsize=(5.5, 4.5))
filled = ax.contourf(dayfine, dayfine, crosscor.T, levels=np.linspace(-1, 1, 21), cmap="RdBu_r")
ax.contour(dayfine, dayfine, crosscor.T, levels=[0.0], colors="black", linewidths=1)
fig.colorbar(filled, ax=ax, label="correlation")
ax.set_xlabel("Temperature day (s)")
ax.set_ylabel("Log precipitation day (t)")
fig

# %% [markdown]
# ### Figure 6.3
# Standard deviation of the log precipitation residuals across stations for each
# day (points) with its smooth, the diagonal of the residual covariance Sigma_e.
# The log variance is smoothed so the fitted standard deviation stays positive.

# %%
logprecres = logprec - np.asarray(logprecfd(day5))
logprecvar1 = np.sum(logprecres**2, axis=1) / (logprec.shape[1] - 1)
logvarfit = fdt.smooth(
    np.log(logprecvar1),
    day5,
    basis=fdt.Fourier(domain=(0.0, 365.0), n_basis=365),
    penalty=harmonic,
)
print(f"log variance smooth: GCV lambda = {logvarfit.lam:.3g}, df = {logvarfit.df:.2f}")
sd_smooth = np.exp(np.asarray(logvarfit.fd(day5)) / 2.0)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(day5, np.sqrt(logprecvar1), "o", markersize=2, fillstyle="none", color="grey")
ax.plot(day5, sd_smooth, "k-", linewidth=2)
ax.set_xlabel("Day (July 1 to June 30)")
ax.set_ylabel("Standard deviation of log10 precipitation")
fig

# %% [markdown]
# ### Figure 6.4
# Phase-plane plot of the simple harmonic x(t) = sin(2 pi t): acceleration against
# velocity traces an ellipse; kinetic energy peaks on the horizontal axis and
# potential energy on the vertical axis.

# %%
unitbasis = fdt.Fourier(domain=(0.0, 1.0), n_basis=3)
sinscale = 1.0 / float(np.asarray(unitbasis(np.array([0.25])))[0, 1])
sinefd = fdt.FData(np.array([[0.0], [sinscale], [0.0]]), unitbasis)

fig, ax = plt.subplots(figsize=(5, 5))
fdt.phase_plane(
    sinefd,
    np.linspace(0.0, 1.0, 201),
    labels={0.0: "t=0", 0.25: "t=1/4", 0.5: "t=1/2", 0.75: "t=3/4"},
    ax=ax,
    color="black",
)
ax.axhline(0.0, color="grey", linestyle=":")
ax.axvline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Velocity")
ax.set_ylabel("Acceleration")
ax.set_title(r"Phase plane of $\sin(2\pi t)$")
fig

# %% [markdown]
# ### Figure 6.5
# Log10 of the US nondurable goods index, 1960 to 1970, with a smooth by order-6
# B-splines with a knot at every month and a D^4 penalty.

# %%
nondur = fdt.datasets.load_nondurables()
years = 1919.0 + np.arange(nondur.value.size) / 12.0
window = (years >= 1960.0) & (years < 1971.0)
ndyears = years[window]
lognondur = np.log10(nondur.value[window])
ndbasis = fdt.BSpline(domain=(ndyears[0], ndyears[-1]), order=6, breaks=ndyears)
ndfit = fdt.smooth(lognondur, ndyears, basis=ndbasis, lam=1e-7, penalty=4)
ndfine = np.linspace(ndyears[0], ndyears[-1], 1201)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(ndyears, lognondur, "o", markersize=2, fillstyle="none", color="grey")
ax.plot(ndfine, np.asarray(ndfit.fd(ndfine)), "k-", linewidth=1.5)
ax.set_xlabel("Year")
ax.set_ylabel(r"$\log_{10}$ nondurable goods index")
fig

# %% [markdown]
# ### Figure 6.6
# Phase-plane plot of the nondurable goods index for 1964: acceleration against
# velocity of the smooth in Figure 6.5, labelled with the months.

# %%
months = "JFMAMJJASOND"
year64 = np.linspace(1964.0, 1965.0, 241)

fig, ax = plt.subplots(figsize=(5.5, 5))
fdt.phase_plane(
    ndfit.fd,
    year64,
    labels={1964.0 + (m + 0.5) / 12.0: months[m] for m in range(12)},
    ax=ax,
    color="black",
)
ax.axhline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Velocity")
ax.set_ylabel("Acceleration")
ax.set_title("Nondurable goods index, 1964")
fig

# %% [markdown]
# ### Figure 6.7
# Phase-plane plots of height acceleration against velocity for the first ten
# Berkeley girls from age 3 to 18 (order-6 B-splines, D^4 penalty, lambda = 1); the
# pubertal spurt is the large loop, and age 11.5 is marked on each curve.

# %%
growth = fdt.datasets.load_growth()
growthbasis = fdt.BSpline(domain=(1.0, 18.0), order=6, breaks=growth.age)
girlsfd = fdt.smooth(growth.hgtf, growth.age, basis=growthbasis, lam=1.0, penalty=4).fd[:10]

fig, ax = plt.subplots(figsize=(6, 5))
fdt.phase_plane(girlsfd, np.linspace(3.0, 18.0, 301), labels={11.5: "o"}, ax=ax, linewidth=1)
ax.axhline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Velocity (cm/yr)")
ax.set_ylabel(r"Acceleration (cm/yr$^2$)")
fig

# %% [markdown]
# ### Figure 6.8
# Log precipitation at Prince Rupert with its smooth and 95% pointwise confidence
# limits, from the linear data-to-coefficient map and the smoothed residual
# variance of Figure 6.3 (Sigma_e taken as diagonal).

# %%
rupert = weather.stations.index("Pr. Rupert")
sigma_e = np.exp(np.asarray(logvarfit.fd(day5)))
y2c = np.asarray(logprecfit.y2c_map)
phi = np.asarray(daybasis(day5))
coef_var = (y2c * sigma_e) @ y2c.T
fit_se = np.sqrt(np.einsum("ij,jk,ik->i", phi, coef_var, phi))
rupertfit = np.asarray(logprecfd(day5))[:, rupert]

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(day5, logprec[:, rupert], "o", markersize=2, fillstyle="none", color="grey")
ax.plot(day5, rupertfit, "k-", linewidth=2, label="smooth")
ax.plot(day5, rupertfit + 2.0 * fit_se, "k--", linewidth=1, label="95% limits")
ax.plot(day5, rupertfit - 2.0 * fit_se, "k--", linewidth=1)
ax.set_xlabel("Day (July 1 to June 30)")
ax.set_ylabel(r"$\log_{10}$ precipitation (mm)")
ax.set_title("Prince Rupert")
ax.legend()
fig
