# %% [markdown]
# # Chapter 7: Exploring Variation: Functional Principal and Canonical Components Analysis
#
# Figures 7.1 to 7.6 of Ramsay, Hooker and Graves (2009): principal components of
# the log precipitation curves (unrotated and varimax rotated) and their scores,
# the bivariate PCA of the handwriting data, and the canonical correlation analysis
# of temperature and log precipitation.

# %%
import matplotlib.pyplot as plt
import numpy as np

import fdatools as fdt

plt.rcParams["figure.max_open_warning"] = 0

weather = fdt.datasets.load_canadian_weather()
day5 = weather.t - 0.5  # the book's day.5: mid-day points 0.5, ..., 364.5
harmonic = fdt.LDO.harmonic(period=365.0)
daybasis = fdt.Fourier(domain=(0.0, 365.0), n_basis=365)
logprecfd = fdt.smooth(weather.log10precip, day5, basis=daybasis, lam=1e6, penalty=harmonic).fd
tempfd = fdt.smooth(weather.temp, day5, basis=daybasis, lam=1e4, penalty=harmonic).fd
stations = weather.stations

# %% [markdown]
# ### Figure 7.1
# The first four principal components of the log precipitation curves, each shown as
# the mean curve (solid) plus (dashed) and minus (dotted) a multiple of the harmonic.

# %%
logprecpca = fdt.FPCA(n=4).fit(logprecfd)
print("proportion of variance:", np.round(np.asarray(logprecpca.varprop), 4))

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
logprecpca.plot(ax=axes)
for ax in axes[1]:
    ax.set_xlabel("Day (July 1 to June 30)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\log_{10}$ precipitation")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 7.2
# The four varimax-rotated principal components of the log precipitation curves,
# again as perturbations of the mean; rotation concentrates each component on a
# season.

# %%
logprecrot = logprecpca.rotate("varimax")
print("rotated proportion of variance:", np.round(np.asarray(logprecrot.varprop), 4))

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
logprecrot.plot(ax=axes)
for ax in axes[1]:
    ax.set_xlabel("Day (July 1 to June 30)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\log_{10}$ precipitation")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 7.3
# Scores of the 35 weather stations on the first two varimax-rotated components of
# log precipitation, labelled by station.

# %%
rotscores = np.asarray(logprecrot.scores)

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(rotscores[:, 0], rotscores[:, 1], s=10, color="black")
for name, (x, y) in zip(stations, rotscores[:, :2], strict=True):
    ax.annotate(name, (x, y), fontsize=7, xytext=(2, 2), textcoords="offset points")
ax.axhline(0.0, color="grey", linestyle=":")
ax.axvline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Rotated component I score")
ax.set_ylabel("Rotated component II score")
fig

# %% [markdown]
# ### Figure 7.4
# Bivariate PCA of the 20 handwritten "fda" samples: for each of the first three
# components, the mean (x(t), y(t)) trajectory with arrows to mean + 2 sd x harmonic,
# showing how each component distorts the whole script.

# %%
handwrit = fdt.datasets.load_handwriting()
hwbasis = fdt.BSpline(domain=(0.0, 2300.0), n_basis=105, order=6)
hwfit = fdt.smooth(handwrit.value, handwrit.t, basis=hwbasis, penalty=4)
print(f"handwriting: GCV lambda = {hwfit.lam:.3g}, df = {hwfit.df:.1f}")
hwpca = fdt.FPCA(n=3).fit(hwfit.fd)
hwvarprop = np.asarray(hwpca.varprop)
hwtime = np.linspace(0.0, 2300.0, 401)
hwmean = np.asarray(hwpca.mean_fd(hwtime))[:, 0, :]
hwharm = np.asarray(hwpca.harmonics(hwtime))
hwsd = np.sqrt(np.mean(np.asarray(hwpca.scores) ** 2, axis=0))
arrows = slice(0, None, 8)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for j, ax in enumerate(axes):
    shifted = hwmean + 2.0 * hwsd[j] * hwharm[:, j, :]
    ax.plot(hwmean[:, 0], hwmean[:, 1], "k-", linewidth=1.5)
    ax.plot(shifted[:, 0], shifted[:, 1], color="tab:blue", linewidth=0.8)
    ax.quiver(
        hwmean[arrows, 0],
        hwmean[arrows, 1],
        (shifted - hwmean)[arrows, 0],
        (shifted - hwmean)[arrows, 1],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=0.003,
        color="tab:red",
    )
    ax.set_aspect("equal")
    ax.set_title(f"PC {j + 1} ({100 * hwvarprop[j]:.1f}%)")
    ax.set_xlabel("x (m)")
axes[0].set_ylabel("y (m)")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 7.5
# First pair of canonical weight functions for temperature (solid) and log
# precipitation (dashed), harmonic acceleration penalty lambda = 1e10 on both
# (smaller values give rough weights and canonical correlations near one).

# %%
weathercca = fdt.FCCA(n=3, lam1=1e10, lam2=1e10, penalty=harmonic).fit(tempfd, logprecfd)
canocor = np.asarray(weathercca.correlations)[:3]
print("canonical correlations:", np.round(canocor, 4))
dayfine = np.linspace(0.0, 365.0, 366)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(dayfine, np.asarray(weathercca.weights1(dayfine))[:, 0], "k-", label="temperature")
ax.plot(dayfine, np.asarray(weathercca.weights2(dayfine))[:, 0], "k--", label="log precipitation")
ax.axhline(0.0, color="grey", linestyle=":")
ax.set_xlabel("Day (July 1 to June 30)")
ax.set_ylabel("Canonical weight function")
ax.set_title(f"First canonical pair (r = {canocor[0]:.3f})")
ax.legend()
fig

# %% [markdown]
# ### Figure 7.6
# Scores of the 35 stations on the first canonical variable for log precipitation
# against those for temperature, labelled by station.

# %%
cscore1 = np.asarray(weathercca.scores1)[:, 0]
cscore2 = np.asarray(weathercca.scores2)[:, 0]

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(cscore1, cscore2, s=10, color="black")
for name, x, y in zip(stations, cscore1, cscore2, strict=True):
    ax.annotate(name, (x, y), fontsize=7, xytext=(2, 2), textcoords="offset points")
ax.set_xlabel("Temperature canonical score")
ax.set_ylabel("Log precipitation canonical score")
fig
