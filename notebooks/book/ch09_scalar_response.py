# %% [markdown]
# # Chapter 9: Functional linear models for scalar responses
#
# Figures 9.1 to 9.7. The log10 total annual precipitation at the 35 Canadian weather
# stations is predicted from each station's daily temperature profile:
# y_i = alpha + integral of x_i(t) beta(t) dt + e_i.

# %%
import matplotlib.pyplot as plt
import numpy as np

import fabel as fb

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["figure.autolayout"] = True

weather = fb.datasets.load_canadian_weather()
day = weather.t - 0.5
annual_precip = np.log10(weather.precip.sum(axis=0))
harmonic = fb.LDO.harmonic(365.0)
temp_basis = fb.Fourier(domain=(0.0, 365.0), n_basis=65)
temp = fb.smooth(weather.temp, day, basis=temp_basis, lam=1e-2, penalty=harmonic).fd
covariates = {"const": 1.0, "temp": temp}
t_grid = np.linspace(0.0, 365.0, 366)
month_ticks = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
month_names = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

# %% [markdown]
# ### Figure 9.1
# Estimated regression function beta(t) when beta is expanded in only five Fourier basis
# functions and no roughness penalty is used. The low dimension alone keeps it smooth.

# %%
model_low = fb.fregress(annual_precip, covariates, beta={"temp": fb.Fourier((0.0, 365.0), 5)})
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t_grid, model_low.beta[1](t_grid)[:, 0], "C0", linewidth=2)
ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("beta(t) for temperature")
ax.set_title(f"Five Fourier basis functions, df = {model_low.df:.0f}")
fig

# %% [markdown]
# ### Figure 9.2
# Leave-one-out cross-validation score (sum of squared prediction errors) as a function
# of log10 lambda, for beta in 35 Fourier functions with a harmonic acceleration penalty.

# %%
beta_basis = fb.Fourier((0.0, 365.0), 35)
log_lambdas = np.arange(6.0, 15.01, 0.5)
cv_scores = []
for log_lambda in log_lambdas:
    trial = fb.fregress(
        annual_precip, covariates, beta={"temp": (beta_basis, 10.0**log_lambda, harmonic)}
    )
    cv_scores.append(trial.cv().sse)
best_log_lambda = float(log_lambdas[int(np.argmin(cv_scores))])
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(log_lambdas, cv_scores, "o-", color="C0")
ax.axvline(best_log_lambda, color="C3", linestyle=":")
ax.set_xlabel("log10 smoothing parameter lambda")
ax.set_ylabel("Cross-validation score")
ax.set_title(f"Minimum at log10 lambda = {best_log_lambda:.1f}")
fig

# %% [markdown]
# ### Figure 9.3
# The regression function beta(t) estimated with the roughness penalty at the
# cross-validated lambda: temperature in late autumn and winter predicts precipitation.

# %%
model = fb.fregress(
    annual_precip, covariates, beta={"temp": (beta_basis, 10.0**best_log_lambda, harmonic)}
)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t_grid, model.beta[1](t_grid)[:, 0], "C0", linewidth=2)
ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("beta(t) for temperature")
ax.set_title(f"Harmonic acceleration penalty, df = {model.df:.2f}")
fig

# %% [markdown]
# ### Figure 9.4
# Observed log10 annual precipitation against the values fitted by the penalised model.
# The diagonal is perfect prediction; a few coastal stations are the largest misses.

# %%
fitted = np.asarray(model.fitted)
r_squared = 1.0 - np.sum((annual_precip - fitted) ** 2) / np.sum(
    (annual_precip - annual_precip.mean()) ** 2
)
fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.plot(fitted, annual_precip, "o", color="C0")
limits = [min(fitted.min(), annual_precip.min()), max(fitted.max(), annual_precip.max())]
ax.plot(limits, limits, "k--", linewidth=1)
for k in np.argsort(np.abs(annual_precip - fitted))[-3:]:
    ax.annotate(weather.stations[k], (fitted[k], annual_precip[k]), fontsize=8)
ax.set_xlabel("Fitted log10 annual precipitation")
ax.set_ylabel("Observed log10 annual precipitation")
ax.set_title(f"R$^2$ = {r_squared:.2f}")
fig

# %% [markdown]
# ### Figure 9.5
# The penalised estimate of beta(t) with pointwise 95% confidence limits (plus and minus
# two standard errors), from the error variance SSE / (n - df).

# %%
stderr = model.stderr()
beta_values = model.beta[1](t_grid)[:, 0]
beta_se = stderr.beta[1](t_grid)[:, 0]
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t_grid, beta_values, "C0", linewidth=2)
ax.plot(t_grid, beta_values + 2.0 * beta_se, "C0--", linewidth=1)
ax.plot(t_grid, beta_values - 2.0 * beta_se, "C0--", linewidth=1)
ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("beta(t) for temperature")
fig

# %% [markdown]
# ### Figure 9.6
# Functional principal components regression: log precipitation is regressed on the
# scores of the first four (lightly smoothed) temperature harmonics, and
# beta(t) = sum_k b_k xi_k(t) is shown with plus and minus two standard errors.

# %%
pca = fb.FPCA(n=4, lam=1e5, penalty=harmonic).fit(temp)
scores = np.asarray(pca.scores)
score_terms = {"const": 1.0} | {f"pc{k + 1}": scores[:, k] for k in range(4)}
pcr = fb.fregress(annual_precip, score_terms)
score_coefs = np.asarray(pcr.coefficients)[1:]
score_cov = np.asarray(pcr.stderr().cov)[1:, 1:]
harmonic_values = np.asarray(pca.harmonics(t_grid))
pcr_beta = harmonic_values @ score_coefs
pcr_se = np.sqrt(np.einsum("tk,kl,tl->t", harmonic_values, score_cov, harmonic_values))
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t_grid, pcr_beta, "C2", linewidth=2, label="principal components regression")
ax.plot(t_grid, pcr_beta + 2.0 * pcr_se, "C2--", linewidth=1)
ax.plot(t_grid, pcr_beta - 2.0 * pcr_se, "C2--", linewidth=1)
ax.plot(t_grid, beta_values, "C0", linewidth=1, label="roughness penalty")
ax.axhline(0.0, color="k", linestyle=":", linewidth=1)
ax.set_xticks(month_ticks, month_names)
ax.set_xlabel("Day")
ax.set_ylabel("beta(t) for temperature")
ax.legend(loc="upper left", fontsize=8)
fig

# %% [markdown]
# ### Figure 9.7
# Permutation test of no relation between temperature and precipitation: histogram of
# the F statistic over 200 random permutations of the responses, the observed F (solid)
# and the 95% permutation quantile (dashed).

# %%
test = fb.stats.f_test(
    annual_precip,
    [np.ones(temp.n_curves), temp],
    basis=[None, beta_basis],
    lam=[0.0, 10.0**best_log_lambda],
    penalty=[2, harmonic],
    n_perm=200,
    random_state=2009,
)
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(test.null, bins=25, color="0.7", edgecolor="0.4")
ax.axvline(test.statistic, color="C3", linewidth=2, label=f"observed F = {test.statistic:.2f}")
ax.axvline(test.critical_value, color="k", linestyle="--", label="95% permutation quantile")
ax.set_xlabel("F statistic")
ax.set_ylabel("Frequency")
ax.set_title(f"Permutation p-value = {test.pvalue:.3f}")
ax.legend(fontsize=8)
fig
