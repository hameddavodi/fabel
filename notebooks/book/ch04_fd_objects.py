# %% [markdown]
# # Chapter 4: How to build functional data objects
#
# Figures 4.1 to 4.3 of Ramsay, Hooker & Graves (2009). A functional data object
# is a coefficient array plus a basis. This chapter builds such objects for the
# Canadian weather data, by smoothing and by plain regression, and applies a
# linear differential operator to them.

# %%
import matplotlib.pyplot as plt
import numpy as np

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0

# %% [markdown]
# ### Figure 4.1
# Mean daily temperature curves for the 35 Canadian weather stations, estimated
# by smoothing the daily data with a 65-term Fourier basis.

# %%
weather = fb.datasets.load_canadian_weather()
day = np.asarray(weather.t, dtype=float) - 0.5
temp_basis = fb.Fourier(domain=(0.0, 365.0), n_basis=65)
harmonic = fb.LDO.harmonic(365.0)
temp_fd = fb.smooth(weather.temp, day, basis=temp_basis, lam=1e4, penalty=harmonic).fd

fig, ax = plt.subplots(figsize=(8, 5))
temp_fd.plot(ax=ax, color="k", linewidth=0.6)
ax.set_xlabel("Day (January 1 to December 31)")
ax.set_ylabel("Mean temperature (deg C)")
fig

# %% [markdown]
# ### Figure 4.2
# Montreal's daily mean temperatures (circles) and the curve fitted by ordinary
# least-squares regression on a 13-term Fourier basis (no roughness penalty).

# %%
montreal = weather.stations.index("Montreal")
regression_basis = fb.Fourier(domain=(0.0, 365.0), n_basis=13)
design = regression_basis(day)
montreal_coefs, *_ = np.linalg.lstsq(design, weather.temp[:, montreal], rcond=None)
montreal_fd = fb.FData(montreal_coefs, regression_basis)

fig, ax = plt.subplots(figsize=(8, 5))
montreal_fd.plot_fit(weather.temp[:, montreal], day, ax=ax, color="k", linewidth=1.5)
ax.set_xlabel("Day (January 1 to December 31)")
ax.set_ylabel("Temperature (deg C)")
ax.set_title("Montreal")
fig

# %% [markdown]
# ### Figure 4.3
# The harmonic acceleration L x = (2 pi / 365)^2 D x + D^3 x of the 35 smoothed
# temperature curves. L removes pure sinusoids of period one year, so what is
# left is the part of each curve that is not a simple annual cycle.

# %%
year_grid = np.linspace(0.0, 365.0, 366)
accelerated = temp_fd(year_grid, harmonic)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(year_grid, accelerated, color="k", linewidth=0.6)
ax.axhline(0.0, color="0.5", linestyle=":", linewidth=0.8)
ax.set_xlabel("Day (January 1 to December 31)")
ax.set_ylabel("L-temperature")
fig
