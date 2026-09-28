# %% [markdown]
# # Chapter 3: How to specify basis systems for building functions
#
# Figures 3.1 to 3.4 of Ramsay, Hooker & Graves (2009). A function is a weighted
# sum of basis functions; this chapter looks at B-spline bases: how a spline is
# built from its basis, how the order controls smoothness, and how well a spline
# basis approximates a known function and its derivatives.

# %%
import matplotlib.pyplot as plt
import numpy as np

import fdatools as fdt

plt.rcParams["figure.max_open_warning"] = 0

# %% [markdown]
# ### Figure 3.1
# A cubic (order 4) spline function on [0, 10] with nine equally spaced interior
# knots (dotted lines). The coefficients are a fixed smooth sequence chosen for
# illustration, as in the book, where the spline is built from its coefficients.

# %%
spline_basis = fdt.BSpline(domain=(0.0, 10.0), n_basis=13, order=4)
knot_points = np.asarray(spline_basis.breaks)
spline_coefs = np.array([0.0, 1.5, 2.8, 1.2, -1.0, -2.2, -0.5, 1.8, 2.6, 1.0, -0.8, -1.5, 0.2])
spline_fd = fdt.FData(spline_coefs, spline_basis)
x_grid = np.linspace(0.0, 10.0, 501)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(x_grid, spline_fd(x_grid)[:, 0], "k-", linewidth=1.5)
for knot in knot_points[1:-1]:
    ax.axvline(knot, color="0.5", linestyle=":", linewidth=0.8)
ax.set_xlabel("t")
ax.set_ylabel("x(t)")
fig

# %% [markdown]
# ### Figure 3.2
# The thirteen B-spline basis functions of order 4 with nine equally spaced
# interior knots (dashed lines) on [0, 10].

# %%
phi = spline_basis(x_grid)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(x_grid, phi, "k-", linewidth=1.0)
for knot in knot_points[1:-1]:
    ax.axvline(knot, color="0.5", linestyle="--", linewidth=0.8)
ax.set_xlabel("t")
ax.set_ylabel("B-spline basis functions")
fig

# %% [markdown]
# ### Figure 3.3
# B-spline bases of orders 1 (step functions) to 4 (cubic) on the same breaks.
# Each order adds one continuous derivative at the knots. fdatools requires distinct
# breaks, so the book's view of smoothness at a knot is shown through the order.

# %%
breaks = np.linspace(0.0, 1.0, 5)
unit_grid = np.linspace(0.0, 1.0, 401)

fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
for order, axis in zip(range(1, 5), axes.ravel(), strict=True):
    basis = fdt.BSpline(domain=(0.0, 1.0), breaks=breaks, order=order)
    axis.plot(unit_grid, basis(unit_grid), linewidth=1.0)
    for knot in breaks[1:-1]:
        axis.axvline(knot, color="0.5", linestyle=":", linewidth=0.8)
    axis.set_title(f"Order {order}: {basis.n_basis} basis functions")
for axis in axes[1]:
    axis.set_xlabel("t")
fig.tight_layout()
fig

# %% [markdown]
# ### Figure 3.4
# Approximating sin(t) on [0, 2 pi] by least squares with 13 order-4 B-splines:
# the error of the approximation (top) and of its first (middle) and second
# (bottom) derivatives. Errors grow with each derivative taken.

# %%
circle = np.linspace(0.0, 2.0 * np.pi, 201)
sine_basis = fdt.BSpline(domain=(0.0, 2.0 * np.pi), n_basis=13, order=4)
sine_coefs, *_ = np.linalg.lstsq(sine_basis(circle), np.sin(circle), rcond=None)
sine_fd = fdt.FData(sine_coefs, sine_basis)
exact = [np.sin(circle), np.cos(circle), -np.sin(circle)]
titles = ["sin(t)", "D sin(t)", "D$^2$ sin(t)"]

fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
for deriv, axis in enumerate(axes):
    error = sine_fd(circle, deriv=deriv)[:, 0] - exact[deriv]
    axis.plot(circle, error, "k-", linewidth=1.0)
    axis.axhline(0.0, color="0.5", linestyle=":", linewidth=0.8)
    axis.set_ylabel(f"Error in {titles[deriv]}")
axes[-1].set_xlabel("t")
fig.tight_layout()
fig
