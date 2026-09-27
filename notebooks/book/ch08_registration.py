# %% [markdown]
# # Chapter 8: Registration: aligning features for samples of curves
#
# Figures 8.1 to 8.7. The Berkeley growth data (first ten girls) are smoothed, their
# height acceleration curves are aligned first by a landmark (the age of the pubertal
# growth spurt) and then by continuous registration, and the variation is split into
# amplitude and phase parts.

# %%
# ruff: noqa: B018
import matplotlib.pyplot as plt
import numpy as np

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["figure.autolayout"] = True

growth = fb.datasets.load_growth()
age = growth.age
heights = growth.hgtf[:, :10]
height_basis = fb.BSpline(domain=(1.0, 18.0), breaks=age, order=6)
height_fit = fb.smooth(heights, age, basis=height_basis, lam=1e-1, penalty=4)
accel = height_fit.fd.derivative(2)

# Landmark: the age of the pubertal growth spurt (peak velocity, where the
# acceleration crosses zero from above).
spurt_grid = np.linspace(8.0, 16.0, 1601)
pgs_age = spurt_grid[np.argmax(height_fit.fd(spurt_grid, 1), axis=0)]

t_plot = np.linspace(3.0, 18.0, 401)
t_full = np.linspace(1.0, 18.0, 401)
landmark_fit = fb.landmark_register(accel, pgs_age[:, None])
continuous_fit = fb.register(
    landmark_fit.registered,
    landmark_fit.registered.mean(),
    warp_basis=fb.BSpline(domain=(1.0, 18.0), n_basis=7),
    lam=1e-2,
)

# %% [markdown]
# ### Figure 8.1
# Height acceleration of the first ten girls of the Berkeley growth study with the
# cross-sectional mean (heavy dashed). The mean understates the pubertal spurt because
# the spurts happen at different ages. The book smooths with a monotone fit; a
# D^4-penalised order-6 spline is used here (same curves up to small boundary effects).

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t_plot, accel(t_plot), color="0.45", linewidth=1)
ax.plot(t_plot, accel.mean()(t_plot)[:, 0], "k--", linewidth=2.5, label="cross-sectional mean")
ax.plot(pgs_age, np.zeros_like(pgs_age), "o", color="C3", label="PGS landmark")
ax.axhline(0.0, color="0.7", linewidth=0.5)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Height acceleration (cm/yr$^2$)")
ax.set_ylim(-4.0, 2.0)
ax.legend(loc="lower left")
fig

# %% [markdown]
# ### Figure 8.2
# A time-warping function h(t) for one girl (left, the diagonal is no warping) and the
# effect of warping on her acceleration curve (right): the registered curve has its
# pubertal spurt moved to the mean spurt age.

# %%
girl = int(np.argmin(pgs_age))
h = landmark_fit.warp_values(t_full)[:, girl]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
ax1.plot(t_full, h, "C0", linewidth=2)
ax1.plot([1, 18], [1, 18], "k:", linewidth=1)
ax1.plot([pgs_age.mean()], [pgs_age[girl]], "o", color="C3")
ax1.set_xlabel("Clock time t (years)")
ax1.set_ylabel("Warped time h(t)")
ax1.set_title(f"Warping function, girl {girl + 1}")
ax2.plot(t_plot, accel[girl](t_plot)[:, 0], "C0--", label="unregistered")
ax2.plot(t_plot, landmark_fit.registered[girl](t_plot)[:, 0], "C0", label="registered")
ax2.plot(t_plot, landmark_fit.registered.mean()(t_plot)[:, 0], "k", linewidth=2, label="mean")
ax2.axhline(0.0, color="0.7", linewidth=0.5)
ax2.set_ylim(-4.0, 2.0)
ax2.set_xlabel("Age (years)")
ax2.set_ylabel("Acceleration (cm/yr$^2$)")
ax2.legend(loc="lower left")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 8.3
# The ten acceleration curves after landmark registration on the pubertal spurt, with
# their mean (heavy line). The spurt now lines up and the mean keeps its depth.

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t_plot, landmark_fit.registered(t_plot), color="0.45", linewidth=1)
ax.plot(t_plot, landmark_fit.registered.mean()(t_plot)[:, 0], "k", linewidth=2.5)
ax.plot(t_plot, accel.mean()(t_plot)[:, 0], "k--", linewidth=1.5)
ax.axvline(pgs_age.mean(), color="C3", linestyle=":")
ax.axhline(0.0, color="0.7", linewidth=0.5)
ax.set_ylim(-4.0, 2.0)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Height acceleration (cm/yr$^2$)")
ax.set_title("Landmark-registered curves (solid mean) and unregistered mean (dashed)")
fig

# %% [markdown]
# ### Figure 8.4
# The landmark warping functions h_i(t) (left) and their deformations h_i(t) - t
# (right). A girl with an early spurt has a warp below the diagonal around puberty.

# %%
warps = landmark_fit.warp_values(t_full)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
ax1.plot(t_full, warps)
ax1.plot([1, 18], [1, 18], "k:", linewidth=1)
ax1.set_xlabel("Age (years)")
ax1.set_ylabel("h(t)")
ax1.set_title("Landmark warping functions")
ax2.plot(t_full, warps - t_full[:, None])
ax2.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax2.set_xlabel("Age (years)")
ax2.set_ylabel("h(t) - t")
ax2.set_title("Deformations")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 8.5
# Continuous registration (minimum-eigenvalue criterion) of the landmark-registered
# curves to their mean, with a 7-function cubic B-spline basis for W and lambda = 0.01.
# The mean of the registered curves is the heavy line.

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t_plot, continuous_fit.registered(t_plot), color="0.45", linewidth=1)
ax.plot(t_plot, continuous_fit.registered.mean()(t_plot)[:, 0], "k", linewidth=2.5)
ax.axhline(0.0, color="0.7", linewidth=0.5)
ax.set_ylim(-4.0, 2.0)
ax.set_xlabel("Age (years)")
ax.set_ylabel("Height acceleration (cm/yr$^2$)")
ax.set_title("Continuously registered acceleration curves")
fig

# %% [markdown]
# ### Figure 8.6
# Warping functions of the continuous registration (left) and their deformations
# h(t) - t (right). They are small because the landmark step already removed most of
# the phase variation.

# %%
warps_c = continuous_fit.warp_values(t_full)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
ax1.plot(t_full, warps_c)
ax1.plot([1, 18], [1, 18], "k:", linewidth=1)
ax1.set_xlabel("Age (years)")
ax1.set_ylabel("h(t)")
ax1.set_title("Continuous warping functions")
ax2.plot(t_full, warps_c - t_full[:, None])
ax2.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax2.set_xlabel("Age (years)")
ax2.set_ylabel("h(t) - t")
ax2.set_title("Deformations")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 8.7
# Amplitude and phase decomposition (Kneip and Ramsay 2008): the three cross-sectional
# means (unregistered, landmark, landmark + continuous) and, per registration, the
# share of the variation that was due to phase (RSQ). The continuous step starts from the
# landmark-registered curves, so its phase part is small.

# %%
landmark_parts = landmark_fit.decompose(domain=(3.0, 18.0))
continuous_parts = continuous_fit.decompose(domain=(3.0, 18.0))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [2, 1]})
ax1.plot(t_plot, accel.mean()(t_plot)[:, 0], "k--", label="unregistered")
ax1.plot(t_plot, landmark_fit.registered.mean()(t_plot)[:, 0], "C0", label="landmark")
ax1.plot(t_plot, continuous_fit.registered.mean()(t_plot)[:, 0], "C3", label="continuous")
ax1.axhline(0.0, color="0.7", linewidth=0.5)
ax1.set_xlabel("Age (years)")
ax1.set_ylabel("Mean acceleration (cm/yr$^2$)")
ax1.legend(loc="lower left")
names = ["landmark", "continuous\n(after landmark)"]
position = np.arange(2)
amp = [landmark_parts.amp_mse, continuous_parts.amp_mse]
phase = [landmark_parts.phase_mse, continuous_parts.phase_mse]
ax2.bar(position - 0.2, amp, width=0.4, color="C0", label="amplitude MSE")
ax2.bar(position + 0.2, phase, width=0.4, color="C1", label="phase MSE")
for k, parts in enumerate([landmark_parts, continuous_parts]):
    ax2.text(k, max(amp[k], phase[k]) * 1.05, f"RSQ={parts.rsq:.2f}", ha="center")
ax2.set_xticks(position, names)
ax2.axhline(0.0, color="k", linewidth=0.5)
ax2.set_ylabel("Mean squared error")
ax2.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig
