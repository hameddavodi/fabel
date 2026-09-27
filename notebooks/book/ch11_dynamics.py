# %% [markdown]
# # Chapter 11: Functional models and dynamics
#
# Figures 11.1 to 11.9. Principal differential analysis (PDA) of the lip movement data,
# a coupled third-order model of handwriting, a first-order forced model of the oil
# refinery data and gradient matching for a simulated stirred-tank chemical reactor.

# %%
# ruff: noqa: B018
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["figure.autolayout"] = True

lip = fb.datasets.load_lip()
lip_domain = (float(lip.t[0]), float(lip.t[-1]))
lip_basis = fb.BSpline(domain=lip_domain, breaks=lip.t, order=6)
lip_fd = fb.smooth(lip.value, lip.t, basis=lip_basis, lam=1e-12, penalty=4).fd
lip_grid = np.linspace(lip_domain[0], lip_domain[1], 201)
lip_pda = fb.PDA(order=2, weight_basis=fb.BSpline(domain=lip_domain, n_basis=7)).fit(lip_fd)

# %% [markdown]
# ### Figure 11.1
# Lower lip position during 20 repetitions of the syllable "bob" (left) and the lip
# acceleration (right), from an order-6 spline with a D^4 roughness penalty.

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(lip_grid, lip_fd(lip_grid), linewidth=0.8)
axes[0].set_ylabel("Lip position (mm)")
axes[1].plot(lip_grid, lip_fd(lip_grid, 2), linewidth=0.8)
axes[1].set_ylabel("Lip acceleration (mm/s$^2$)")
for ax in axes:
    ax.set_xlabel("Time (s)")
fig

# %% [markdown]
# ### Figure 11.2
# Phase-plane plot of the mean lip curve: acceleration against velocity. The loops
# show the exchange between kinetic and potential energy as the lip opens and closes.

# %%
fig, ax = plt.subplots(figsize=(6, 5.5))
fb.phase_plane(
    lip_fd.mean(),
    lip_grid,
    labels={0.05: "0.05 s", 0.15: "0.15 s", 0.25: "0.25 s", 0.33: "0.33 s"},
    ax=ax,
    color="C0",
    linewidth=2,
)
ax.axhline(0.0, color="0.6", linewidth=0.5)
ax.axvline(0.0, color="0.6", linewidth=0.5)
ax.set_xlabel("Velocity (mm/s)")
ax.set_ylabel("Acceleration (mm/s$^2$)")
fig

# %% [markdown]
# ### Figure 11.3
# The PDA weight functions of D^2 x + beta1(t) Dx + beta0(t) x = 0 estimated from the
# lip curves (7 cubic B-spline functions each): the stiffness beta0 and damping beta1.

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(lip_grid, lip_pda.weights_[0](lip_grid)[:, 0], "C0", linewidth=2)
axes[0].set_title("beta0(t)")
axes[1].plot(lip_grid, lip_pda.weights_[1](lip_grid)[:, 0], "C1", linewidth=2)
axes[1].set_title("beta1(t)")
for ax in axes:
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_xlabel("Time (s)")
fig

# %% [markdown]
# ### Figure 11.4
# Stability diagram of the lip model: the path (beta1(t), beta0(t)) against the parabola
# beta0 = beta1^2 / 4. Above it the system oscillates; negative beta1 amplifies.

# %%
fig, ax = plt.subplots(figsize=(6.5, 5.5))
lip_pda.plot_overlay(
    ax=ax, labels={0.0: "0 s", 0.1: "0.1 s", 0.2: "0.2 s", 0.3: "0.3 s"}, color="C0", linewidth=2
)
ax.set_title("PDA stability diagram, lip data")
fig

# %% [markdown]
# ### Figure 11.5
# Solutions of the estimated equation L x = 0 started from each curve's initial
# position and velocity (dashed) compared with the smoothed lip curves (solid), for
# four repetitions.

# %%
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True)
for ax, rep in zip(axes.ravel(), [0, 5, 10, 15], strict=True):
    start = [lip_fd(lip_grid[:1])[0, rep], lip_fd(lip_grid[:1], 1)[0, rep]]
    ax.plot(lip_grid, lip_fd(lip_grid)[:, rep], "C0", linewidth=2, label="data")
    ax.plot(lip_grid, lip_pda.solve(lip_grid, start), "C3--", linewidth=2, label="L x = 0")
    ax.set_title(f"Repetition {rep + 1}")
axes[0, 0].legend(fontsize=8)
for ax in axes[1]:
    ax.set_xlabel("Time (s)")
for ax in axes[:, 0]:
    ax.set_ylabel("Lip position (mm)")
fig

# %% [markdown]
# ### Figure 11.6
# A coupled third-order PDA of handwriting, D^3 x_i = -sum_k sum_j beta_ikj(t) D^j x_k,
# fitted to 20 samples and solved from the mean curve's initial state. The book uses
# Chinese script, which Fabel does not ship; the cursive "fda" data stand in for it.

# %%
writing = fb.datasets.load_handwriting()
pen_domain = (0.0, 2300.0)
pen_basis = fb.BSpline(domain=pen_domain, n_basis=105, order=6)
pen_fit = fb.smooth(writing.value[::2], writing.t[::2], basis=pen_basis, lam=1e2, penalty=4)
pen_pda = fb.PDA(order=3, weight_basis=fb.BSpline(domain=pen_domain, n_basis=43), n_grid=1001).fit(
    pen_fit.fd
)
pen_grid = np.linspace(pen_domain[0], pen_domain[1], 461)
pen_mean = pen_fit.fd.mean()
mean_values = pen_mean(pen_grid)[:, 0, :]
pen_start = np.array([[pen_mean(pen_grid[:1], j)[0, 0, v] for j in range(3)] for v in range(2)])
pen_solution = pen_pda.solve(pen_grid, pen_start)
fig, axes = plt.subplots(1, 3, figsize=(13, 4), gridspec_kw={"width_ratios": [1.4, 1, 1]})
axes[0].plot(mean_values[:, 0], mean_values[:, 1], "C0", linewidth=2, label="mean script")
axes[0].plot(pen_solution[:, 0], pen_solution[:, 1], "C3--", linewidth=1.5, label="L x = 0")
axes[0].set_aspect("equal")
axes[0].legend(fontsize=8)
axes[0].set_title("Script")
for k, name in enumerate(["X", "Y"]):
    axes[k + 1].plot(pen_grid, mean_values[:, k], "C0", linewidth=2)
    axes[k + 1].plot(pen_grid, pen_solution[:, k], "C3--", linewidth=1.5)
    axes[k + 1].set_title(f"{name} coordinate")
    axes[k + 1].set_xlabel("Time (ms)")
fig

# %% [markdown]
# ### Figure 11.7
# The refinery data: a step drop in reflux flow (bottom) and the tray 47 level (top),
# with the solution of the forced first-order model Dx = -beta x + alpha u. The
# constant coefficients come from a concurrent regression of Dx on x and u with
# constant coefficient bases (PDA forcing functions are not in Fabel).

# %%
refinery = fb.datasets.load_refinery()
ref_domain = (0.0, 193.0)
reflux_fd = fb.smooth(
    refinery.reflux,
    refinery.time,
    basis=fb.BSpline(ref_domain, order=1, breaks=[0.0, 67.0, 193.0]),
    lam=0.0,
).fd
tray_breaks = np.unique(np.r_[np.linspace(0.0, 67.0, 5), np.linspace(70.0, 193.0, 30)])
tray_basis = fb.BSpline(ref_domain, breaks=tray_breaks, order=4)
tray_fd = fb.smooth(refinery.tray47, refinery.time, basis=tray_basis, lam=1.0).fd
forced = fb.fregress(
    tray_fd.derivative(), {"tray": tray_fd, "reflux": reflux_fd}, beta=fb.Constant(ref_domain)
)
beta_hat = -float(np.asarray(forced.beta[0].coefs)[0, 0])
alpha_hat = float(np.asarray(forced.beta[1].coefs)[0, 0])
tray_solution = solve_ivp(
    lambda s, x: [-beta_hat * x[0] + alpha_hat * float(reflux_fd(np.array([s]))[0, 0])],
    ref_domain,
    [float(tray_fd(np.array([0.0]))[0, 0])],
    t_eval=refinery.time,
    max_step=1.0,
).y[0]
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
axes[0].plot(refinery.time, refinery.tray47, "o", color="0.5", markersize=3, label="data")
axes[0].plot(refinery.time, tray_solution, "C3", linewidth=2, label="fitted ODE")
axes[0].set_ylabel("Tray 47 level")
axes[0].set_title(f"beta = {beta_hat:.4f}, alpha = {alpha_hat:.3f}")
axes[0].legend(fontsize=8)
axes[1].plot(refinery.time, refinery.reflux, "o", color="0.5", markersize=3)
axes[1].plot(refinery.time, reflux_fd(refinery.time)[:, 0], "C0", linewidth=2)
axes[1].set_ylabel("Reflux flow")
axes[1].set_xlabel("Time (min)")
fig

# %% [markdown]
# ### Figure 11.8
# A continuously stirred tank reactor: input flow and input temperature (bottom) drive
# output concentration C and temperature T (top). The book's CSTR measurements are not
# in Fabel, so a dimensionless CSTR model is simulated (fixed seed) and noise is added.


# %%
def flow(s):
    return 1.0 - 0.3 * (16.0 <= s < 32.0) + 0.25 * (40.0 <= s < 52.0)


def inlet_temp(s):
    return 1.0 + 0.15 * (8.0 <= s < 24.0) - 0.1 * (48.0 <= s < 60.0)


activation, inlet_conc, coolant = 6.0, 1.0, 0.9
true_params = (0.8, 0.6, 1.0)  # rate constant, heat of reaction, heat transfer


def reactor(s, state, rate, heat, transfer):
    conc, temp = state
    k = rate * np.exp(-activation * (1.0 / temp - 1.0))
    return [
        flow(s) * (inlet_conc - conc) - k * conc,
        flow(s) * (inlet_temp(s) - temp) + heat * k * conc - transfer * (temp - coolant),
    ]


cstr_time = np.linspace(0.0, 64.0, 257)
truth = solve_ivp(
    reactor, (0.0, 64.0), [0.5, 1.0], t_eval=cstr_time, args=true_params, max_step=0.1
).y.T
rng = np.random.default_rng(2009)
observed = truth + rng.normal(scale=[0.01, 0.005], size=truth.shape)
fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex=True)
axes[0, 0].plot(cstr_time, observed[:, 0], ".", color="0.5", markersize=3)
axes[0, 0].plot(cstr_time, truth[:, 0], "C0", linewidth=1.5)
axes[0, 0].set_title("Output concentration C")
axes[0, 1].plot(cstr_time, observed[:, 1], ".", color="0.5", markersize=3)
axes[0, 1].plot(cstr_time, truth[:, 1], "C3", linewidth=1.5)
axes[0, 1].set_title("Output temperature T")
axes[1, 0].plot(cstr_time, [flow(s) for s in cstr_time], "k")
axes[1, 0].set_title("Input flow")
axes[1, 1].plot(cstr_time, [inlet_temp(s) for s in cstr_time], "k")
axes[1, 1].set_title("Input temperature")
for ax in axes[1]:
    ax.set_xlabel("Time (min)")
fig

# %% [markdown]
# ### Figure 11.9
# Gradient matching for the reactor: the noisy outputs are smoothed, the smooth curves
# and their derivatives give the three reactor constants by least squares, and the
# model is re-solved with them (dashed) over the data. Fabel has no nonlinear ODE
# estimator (parameter cascading), so this is built from smooth() and derivatives.

# %%
cstr_fd = fb.smooth(observed, cstr_time, basis=fb.BSpline((0.0, 64.0), n_basis=70), lam=1e-2).fd
inner = np.linspace(1.0, 63.0, 621)
conc_s, temp_s = cstr_fd(inner).T
dconc_s, dtemp_s = cstr_fd(inner, 1).T
flows = np.array([flow(s) for s in inner])
inlets = np.array([inlet_temp(s) for s in inner])
reaction = np.exp(-activation * (1.0 / temp_s - 1.0)) * conc_s
# dC + q (C - Cin) = -rate * e(T) C
rate_hat = float(np.linalg.lstsq(-reaction[:, None], dconc_s + flows * (conc_s - inlet_conc))[0][0])
# dT - q (Tin - T) = heat * rate * e(T) C - transfer (T - Tc)
design = np.column_stack([rate_hat * reaction, -(temp_s - coolant)])
heat_hat, transfer_hat = np.linalg.lstsq(design, dtemp_s - flows * (inlets - temp_s))[0]
refit = solve_ivp(
    reactor,
    (0.0, 64.0),
    list(observed[0]),
    t_eval=cstr_time,
    args=(rate_hat, heat_hat, transfer_hat),
    max_step=0.1,
).y.T
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for k, name in enumerate(["Concentration C", "Temperature T"]):
    axes[k].plot(cstr_time, observed[:, k], ".", color="0.5", markersize=3, label="data")
    axes[k].plot(cstr_time, refit[:, k], "C3--", linewidth=2, label="re-solved ODE")
    axes[k].set_title(name)
    axes[k].set_xlabel("Time (min)")
axes[0].legend(fontsize=8)
fig.suptitle(
    f"estimates (truth): rate {rate_hat:.3f} ({true_params[0]}), "
    f"heat {heat_hat:.3f} ({true_params[1]}), transfer {transfer_hat:.3f} ({true_params[2]})"
)
fig
