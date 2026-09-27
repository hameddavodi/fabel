# %% [markdown]
# # Chapter 10: Linear models for functional responses
#
# Figures 10.1 to 10.8. A functional analysis of variance of the Canadian temperature
# curves by climate region, the concurrent model of log precipitation on temperature,
# and the concurrent model of knee angle on hip angle in the gait data.

# %%
# ruff: noqa: B018
import matplotlib.pyplot as plt
import numpy as np

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["figure.autolayout"] = True

weather = fb.datasets.load_canadian_weather()
day = weather.t - 0.5
n_stations = len(weather.stations)
harmonic = fb.LDO.harmonic(365.0)
day_basis = fb.Fourier(domain=(0.0, 365.0), n_basis=65)
temp = fb.smooth(weather.temp, day, basis=day_basis, lam=1e-2, penalty=harmonic).fd
regions = sorted(set(weather.region))
region_dummy = {
    name: np.array([1.0 if r == name else 0.0 for r in weather.region]) for name in regions
}
t_grid = np.linspace(0.0, 365.0, 366)
month_ticks = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
month_names = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
region_colors = dict(zip(regions, ["C0", "C1", "C2", "C3"], strict=True))

# %% [markdown]
# ### Figure 10.1
# Functional ANOVA of temperature by climate region: the Canada-wide mean function
# and the four regional effects. The effects are made to sum to zero by appending one
# pseudo-observation (a zero curve with every region indicator on), as in the book.

# %%
anova_y = fb.FData(np.hstack([np.asarray(temp.coefs), np.zeros((day_basis.n_basis, 1))]), day_basis)
anova_x = {"const": np.r_[np.ones(n_stations), 0.0]}
anova_x |= {name: np.r_[region_dummy[name], 1.0] for name in regions}
anova = fb.fregress(anova_y, anova_x, beta=fb.Fourier((0.0, 365.0), 11))
fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharex=True)
for ax, name, beta in zip(axes.ravel(), anova.names, anova.beta, strict=False):
    ax.plot(t_grid, beta(t_grid)[:, 0], linewidth=2, color=region_colors.get(name, "k"))
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_title("Canada mean" if name == "const" else f"{name} effect")
    ax.set_xticks(month_ticks, month_names)
    ax.set_ylabel("deg C")
axes.ravel()[-1].axis("off")
fig

# %% [markdown]
# ### Figure 10.2
# Temperature curves of the stations (thin) coloured by region, with the fitted mean
# curve of each region (thick): the prediction of the ANOVA model.

# %%
fig, ax = plt.subplots(figsize=(8, 4.5))
for k, region in enumerate(weather.region):
    ax.plot(t_grid, temp[k](t_grid)[:, 0], color=region_colors[region], linewidth=0.5, alpha=0.6)
for name in regions:
    new_x = {"const": 1.0} | {other: float(other == name) for other in regions}
    prediction = anova.predict(new_x)
    ax.plot(t_grid, prediction(t_grid)[:, 0], color=region_colors[name], linewidth=3, label=name)
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("Mean temperature (deg C)")
ax.legend(fontsize=8)
fig

# %% [markdown]
# ### Figure 10.3
# Permutation F test for a region effect: the observed pointwise F statistic (solid),
# the pointwise 95% permutation quantiles (dotted) and the 95% quantile of the maximum
# of F (dashed), from 200 permutations (treatment coding, Arctic as baseline).

# %%
region_test = fb.stats.f_test(
    temp,
    [np.ones(n_stations)] + [region_dummy[name] for name in regions[1:]],
    basis=fb.Fourier((0.0, 365.0), 11),
    n_perm=200,
    random_state=2009,
)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(region_test.t, region_test.pointwise, "C0", linewidth=2, label="observed F(t)")
ax.plot(region_test.t, region_test.pointwise_critical_value, "k:", label="pointwise 0.95")
ax.axhline(region_test.critical_value, color="k", linestyle="--", label="maximum 0.95")
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("F statistic")
ax.set_title(f"Permutation p-value of max F = {region_test.pvalue:.3f}")
ax.legend(fontsize=8)
fig

# %% [markdown]
# ### Figure 10.4
# The concurrent model log10 precipitation(t) = beta0(t) + beta1(t) temperature(t):
# intercept (left) and temperature coefficient (right) with plus and minus two pointwise
# standard errors, from the residual covariance over days and the smoothing map.

# %%
log_precip = fb.smooth(weather.log10precip, day, basis=day_basis, lam=1e5, penalty=harmonic)
concurrent = fb.fregress(
    log_precip, {"const": 1.0, "temp": temp}, beta=fb.Fourier((0.0, 365.0), 11)
)
residuals = weather.log10precip - np.asarray(concurrent.fitted(day))
concurrent_se = concurrent.stderr(sigma_e=np.cov(residuals))
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, name, beta, se in zip(
    axes, ["Intercept", "Temperature coefficient"], concurrent.beta, concurrent_se.beta, strict=True
):
    values = beta(t_grid)[:, 0]
    width = 2.0 * se(t_grid)[:, 0]
    ax.plot(t_grid, values, "C0", linewidth=2)
    ax.plot(t_grid, values + width, "C0--", linewidth=1)
    ax.plot(t_grid, values - width, "C0--", linewidth=1)
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_xticks(month_ticks, month_names)
    ax.set_title(name)
fig

# %% [markdown]
# ### Figure 10.5
# Leave-one-station-out cross-validated error sum of squares of the concurrent model as
# a function of the harmonic acceleration penalty on both coefficient functions. A
# moderate penalty improves the out-of-sample fit slightly; a heavy one hurts it.

# %%
log_lambdas = np.arange(4.0, 14.01, 1.0)
cv_errors = [
    fb.fregress(
        log_precip,
        {"const": 1.0, "temp": temp},
        beta=fb.Fourier((0.0, 365.0), 11),
        lam=10.0**value,
        penalty=harmonic,
    )
    .cv()
    .sse
    for value in log_lambdas
]
fig, ax = plt.subplots(figsize=(7, 4))
best = int(np.argmin(cv_errors))
ax.plot(log_lambdas, cv_errors, "o-", color="C0")
ax.axvline(log_lambdas[best], color="C3", linestyle=":")
ax.set_title(f"Minimum at log10 lambda = {log_lambdas[best]:.0f}")
ax.set_xlabel("log10 lambda (harmonic acceleration penalty)")
ax.set_ylabel("Cross-validated SSE")
fig

# %% [markdown]
# ### Figure 10.6
# The gait data: hip angle (left) and knee angle (right) of 39 boys over one gait
# cycle, smoothed with a 21-function Fourier basis (time rescaled to [0, 1]).

# %%
gait = fb.datasets.load_gait()
gait_basis = fb.Fourier(domain=(0.0, 1.0), n_basis=21)
gait_penalty = fb.LDO.harmonic(1.0)
hip = fb.smooth(gait.hip_angle, gait.t, basis=gait_basis, lam=1e-11, penalty=gait_penalty).fd
knee_fit = fb.smooth(gait.knee_angle, gait.t, basis=gait_basis, lam=1e-11, penalty=gait_penalty)
knee = knee_fit.fd
cycle = np.linspace(0.0, 1.0, 201)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(cycle, hip(cycle), linewidth=0.7)
axes[0].set_title("Hip angle")
axes[1].plot(cycle, knee(cycle), linewidth=0.7)
axes[1].set_title("Knee angle")
for ax in axes:
    ax.set_xlabel("Time (proportion of gait cycle)")
    ax.set_ylabel("degrees")
fig

# %% [markdown]
# ### Figure 10.7
# Concurrent model knee(t) = beta0(t) + beta1(t) hip(t): the intercept, the hip
# coefficient with plus and minus two standard errors, and the squared multiple
# correlation R^2(t) along the cycle.

# %%
gait_model = fb.fregress(
    knee_fit,
    {"const": 1.0, "hip": hip},
    beta=fb.Fourier((0.0, 1.0), 21),
    lam=1e-8,
    penalty=gait_penalty,
)
gait_residuals = gait.knee_angle - np.asarray(gait_model.fitted(gait.t))
gait_se = gait_model.stderr(sigma_e=np.cov(gait_residuals))
knee_values = knee(cycle)
fit_values = gait_model.fitted(cycle)
r_squared = 1.0 - np.sum((knee_values - fit_values) ** 2, axis=1) / np.sum(
    (knee_values - knee_values.mean(axis=1, keepdims=True)) ** 2, axis=1
)
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, title, beta, se in zip(
    axes[:2], ["Intercept", "Hip coefficient"], gait_model.beta, gait_se.beta, strict=True
):
    values = beta(cycle)[:, 0]
    width = 2.0 * se(cycle)[:, 0]
    ax.plot(cycle, values, "C0", linewidth=2)
    ax.plot(cycle, values + width, "C0--", linewidth=1)
    ax.plot(cycle, values - width, "C0--", linewidth=1)
    ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("Time (proportion of gait cycle)")
axes[2].plot(cycle, r_squared, "C2", linewidth=2)
axes[2].set_ylim(0.0, 1.0)
axes[2].set_title("Squared multiple correlation R$^2$(t)")
axes[2].set_xlabel("Time (proportion of gait cycle)")
fig

# %% [markdown]
# ### Figure 10.8
# Observed knee angle (solid) and the knee angle predicted from the hip angle by the
# concurrent model (dashed) for four boys.

# %%
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True)
for ax, boy in zip(axes.ravel(), [0, 9, 19, 29], strict=True):
    ax.plot(cycle, knee_values[:, boy], "C0", linewidth=2, label="observed")
    ax.plot(cycle, fit_values[:, boy], "C3--", linewidth=2, label="predicted")
    ax.plot(gait.t, gait.knee_angle[:, boy], "o", color="C0", markersize=3)
    ax.set_title(f"Boy {boy + 1}")
axes[0, 0].legend(fontsize=8)
for ax in axes[1]:
    ax.set_xlabel("Time (proportion of gait cycle)")
for ax in axes[:, 0]:
    ax.set_ylabel("Knee angle (degrees)")
fig
