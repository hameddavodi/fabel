# %% [markdown]
# # Chapter 1: Introduction to functional data analysis
#
# Figures 1.1 to 1.10 of Ramsay, Hooker & Graves (2009). The chapter shows the
# kinds of data the book works with: growth curves, an economic index, refinery
# input/output records, gait cycles, handwriting and weather records.

# %%
# ruff: noqa: B018 -- each figure cell ends with a bare `fig` so the notebook shows it
import matplotlib.pyplot as plt
import numpy as np

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0

# %% [markdown]
# ### Figure 1.1
# The heights of the first 10 girls of the Berkeley growth study at 31 ages
# (circles), with the monotone smooth of each record (lines).

# %%
growth = fb.datasets.load_growth()
age = growth.age
heights = growth.hgtf[:, :10]
growth_basis = fb.BSpline(domain=(1.0, 18.0), breaks=age, order=6)
girls = fb.smooth(heights, age, basis=growth_basis, lam=0.01, penalty=3, constraint="monotone")
fine_age = np.linspace(1.0, 18.0, 401)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(fine_age, girls(fine_age), color="0.3", linewidth=1.0)
ax.plot(age, heights, "o", color="k", markersize=3, fillstyle="none")
ax.set_xlabel("Age (years)")
ax.set_ylabel("Height (cm)")
fig

# %% [markdown]
# ### Figure 1.2
# The estimated accelerations of height of the same 10 girls (cm/yr^2); the heavy
# dashed line is their mean. The pubertal growth spurt is the dip near age 11-13.
# A monotone fit is x = b0 + b1 * integral of exp(W), so x'' = x' * W'.

# %%
acceleration = girls(fine_age, deriv=1) * girls.fd(fine_age, deriv=1)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(fine_age, acceleration, color="0.5", linewidth=0.8)
ax.plot(fine_age, acceleration.mean(axis=1), "k--", linewidth=2.5)
ax.axhline(0.0, color="k", linestyle=":", linewidth=0.8)
ax.set_xlim(3.0, 18.0)
ax.set_ylim(-4.0, 2.0)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Acceleration (cm/yr$^2$)")
fig

# %% [markdown]
# ### Figure 1.3
# The monthly index of nondurable goods manufacturing in the United States.
# Dates follow the dataset metadata: monthly values starting January 1919.

# %%
nondurables = fb.datasets.load_nondurables()
index = np.asarray(nondurables.value, dtype=float)
year = 1919.0 + np.arange(index.size) / 12.0

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(year, index, color="k", linewidth=0.8)
ax.set_xlabel("Year")
ax.set_ylabel("Nondurable goods index")
fig

# %% [markdown]
# ### Figure 1.4
# Refinery data: the amount of product on tray 47 of a distillation column (top)
# and the reflux flow into that tray (bottom) during an experiment.

# %%
refinery = fb.datasets.load_refinery()

fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
top.plot(refinery.time, refinery.tray47, "o", color="k", markersize=2.5)
top.set_ylabel("Tray 47 level")
bottom.plot(refinery.time, refinery.reflux, "o", color="k", markersize=2.5)
bottom.set_ylabel("Reflux flow")
bottom.set_xlabel("Time (min)")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 1.5
# Hip and knee angles in the sagittal plane for 39 children over one gait cycle
# (time is the proportion of the cycle). Curves are Fourier smooths of the data.

# %%
gait = fb.datasets.load_gait()
gait_basis = fb.Fourier(domain=(0.0, 1.0), n_basis=21)
gait_fd = fb.smooth(
    np.asarray(gait.value), gait.t, basis=gait_basis, lam=1e-11, penalty=fb.LDO.harmonic(1.0)
).fd
cycle = np.linspace(0.0, 1.0, 201)
gait_values = gait_fd(cycle)

fig, (hip_ax, knee_ax) = plt.subplots(1, 2, figsize=(10, 4.5))
hip_ax.plot(cycle, gait_values[:, :, 0], linewidth=0.7)
hip_ax.set_title("Hip")
hip_ax.set_ylabel("Hip angle (degrees)")
knee_ax.plot(cycle, gait_values[:, :, 1], linewidth=0.7)
knee_ax.set_title("Knee")
knee_ax.set_ylabel("Knee angle (degrees)")
for axis in (hip_ax, knee_ax):
    axis.set_xlabel("Time (proportion of gait cycle)")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 1.6
# Knee angle against hip angle over one cycle for one child (solid) and for the
# mean child (dashed). Letters A-E mark five equally spaced moments of the cycle.

# %%
marks = np.linspace(0.0, 0.8, 5)
child = gait_fd[4](cycle)[:, 0, :]
average = gait_fd.mean()(cycle)[:, 0, :]
child_marks = gait_fd[4](marks)[:, 0, :]
mean_marks = gait_fd.mean()(marks)[:, 0, :]

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(child[:, 0], child[:, 1], "k-", linewidth=1.2)
ax.plot(average[:, 0], average[:, 1], "k--", linewidth=1.2)
for letter, point, centre in zip("ABCDE", child_marks, mean_marks, strict=True):
    ax.text(point[0], point[1], letter, fontsize=12, fontweight="bold")
    ax.text(centre[0], centre[1], letter, fontsize=10, color="0.4")
ax.set_xlabel("Hip angle (degrees)")
ax.set_ylabel("Knee angle (degrees)")
fig

# %% [markdown]
# ### Figure 1.7
# Twenty samples of the script "fda" written by one person: X against Y
# coordinates of the pen tip (meters), each record centered.

# %%
handwriting = fb.datasets.load_handwriting()
pen = np.asarray(handwriting.value)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(pen[:, :, 0], pen[:, :, 1], color="k", linewidth=0.5)
ax.set_aspect("equal")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
fig

# %% [markdown]
# ### Figure 1.8
# Mean daily temperature at four Canadian weather stations: Montreal, Edmonton,
# Prince Rupert and Resolute, smoothed with a 65-term Fourier basis.

# %%
weather = fb.datasets.load_canadian_weather()
chosen = ["Montreal", "Edmonton", "Pr. Rupert", "Resolute"]
columns = [weather.stations.index(name) for name in chosen]
day = np.asarray(weather.t, dtype=float) - 0.5
temp_fd = fb.smooth(
    weather.temp[:, columns],
    day,
    basis=fb.Fourier(domain=(0.0, 365.0), n_basis=65),
    lam=1e4,
    penalty=fb.LDO.harmonic(365.0),
).fd
year_grid = np.linspace(0.0, 365.0, 366)

fig, ax = plt.subplots(figsize=(8, 5))
for name, style, curve in zip(chosen, ["-", "--", "-.", ":"], temp_fd(year_grid).T, strict=True):
    ax.plot(year_grid, curve, "k", linestyle=style, label=name)
ax.set_xlabel("Day")
ax.set_ylabel("Mean temperature (deg C)")
ax.legend(loc="lower center")
fig

# %% [markdown]
# ### Figure 1.9
# The log10 nondurable goods index for 1964-1967, a period of steady growth
# (circles), with a smooth fit (solid) and the straight-line trend (dashed).

# %%
log_index = np.log10(index)
window = (year >= 1962.0) & (year < 1970.0)
window_year = year[window]
window_log = log_index[window]
monthly_breaks = np.linspace(1962.0, 1970.0, 97)
index_fd = fb.smooth(
    window_log,
    window_year,
    basis=fb.BSpline(domain=(1962.0, 1970.0), breaks=monthly_breaks, order=6),
    lam=1e-6,
    penalty=4,
).fd
shown = (window_year >= 1964.0) & (window_year < 1968.0)
fine_year = np.linspace(1964.0, 1968.0, 481)
slope, intercept = np.polyfit(window_year[shown], window_log[shown], 1)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(window_year[shown], window_log[shown], "o", color="k", markersize=3, fillstyle="none")
ax.plot(fine_year, index_fd(fine_year)[:, 0], "k-", linewidth=1.0)
ax.plot(fine_year, intercept + slope * fine_year, "k--", linewidth=1.0)
ax.set_xlabel("Year")
ax.set_ylabel("log10 nondurable goods index")
fig

# %% [markdown]
# ### Figure 1.10
# Phase-plane plot for 1964: acceleration against velocity of the smoothed log10
# nondurable goods index, labelled by month. Seasonal cycles make loops.

# %%
months = "JFMAMJJASOND"
month_labels = {1964.0 + (m + 0.5) / 12.0: months[m] for m in range(12)}
year_1964 = np.linspace(1964.0, 1965.0, 201)

fig, ax = plt.subplots(figsize=(6, 6))
fb.phase_plane(index_fd, year_1964, labels=month_labels, ax=ax, color="k")
ax.axhline(0.0, color="k", linestyle=":", linewidth=0.8)
ax.axvline(0.0, color="k", linestyle=":", linewidth=0.8)
ax.set_xlabel("Velocity")
ax.set_ylabel("Acceleration")
fig
