# tools/golden_r/stats.R
#
# Golden cases for the `stats` module: var.fd, cor.fd, fbplot, fdepth,
# tperm.fd, Fperm.fd, mean.fd, sd.fd -- all exercised on the real
# CanadianWeather dataset. Driven by tools/make_golden.py, which sources
# common.R first (defines add_case, eval_points, basis_info, finalize,
# GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# smooth.basis, fdPar, var.fd, cor.fd, eval.bifd, fbplot, fdepth, tperm.fd,
# Fperm.fd, fRegress, mean.fd, sd.fd). No fda source is read or copied.
#
# fda API quirks / design decisions discovered while writing this script
# (see PROGRESS.md):
#   - fdepth() requires `data` shaped as list(y = <matrix>) (or
#     list(argvals=, y=)), NOT a bare matrix -- passing a bare matrix fails
#     with "$ operator is invalid for atomic vectors" inside fdepth's body.
#   - Fperm.fd(yfdPar, xfdlist, betalist, ...) with EVERY betalist entry on
#     a constant basis (i.e. every regression coefficient a time-constant
#     scalar, including for a functional response with strong day-of-year
#     structure) produces Fobs/Fvals/Fnull identically 0 at every
#     permutation -- a degenerate result, not a real "no effect" finding
#     (the residual from such a flat model is large: resid_norm ~1446 over
#     365x35 points). Using a genuine functional basis for betalist (same
#     Fourier(25) basis as the response, as in every fRegress example in
#     the public docs) gives a well-behaved, non-degenerate Fobs. We use
#     the latter design (intercept + Atlantic-region dummy, both with
#     functional Fourier betas) for the fRegress/Fperm.fd golden case.

weather_basis_fd <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 25)
  fdobj <- smooth.basis(day, temp, fdPar(basis, 2, 1e2))$fd
  list(day = day, basis = basis, fdobj = fdobj, temp = temp)
}

# ---- var.fd / cor.fd --------------------------------------------------

var_cor_cases <- function(w) {
  day <- w$day
  basis <- w$basis
  temp_fd <- w$fdobj
  precip <- CanadianWeather$dailyAv[, , "Precipitation.mm"]
  precip_fd <- smooth.basis(day, precip, fdPar(basis, 2, 1e4))$fd

  vfd <- var.fd(temp_fd)
  grid <- seq(0, 365, length.out = 30)
  vmat <- eval.bifd(grid, grid, vfd)
  add_case(
    name = "var_fd_weather_temperature",
    r_call = "eval.bifd(grid, grid, var.fd(temp_fd))",
    input = list(grid = grid),
    output = list(varmat = vmat)
  )

  cmat <- cor.fd(grid, temp_fd, grid, precip_fd)
  add_case(
    name = "cor_fd_weather_temp_vs_precip",
    r_call = "cor.fd(grid, temp_fd, grid, precip_fd)",
    input = list(grid = grid),
    output = list(cormat = cmat)
  )
}

# ---- fbplot -------------------------------------------------------------

fbplot_case <- function(w) {
  grid <- seq(0, 365, length.out = 365)
  Y <- eval.fd(grid, w$fdobj)
  fb <- fbplot(Y, plot = FALSE, method = "MBD")
  add_case(
    name = "fbplot_weather_temperature_mbd",
    r_call = "fbplot(eval.fd(seq(0,365,length.out=365), temp_fd), plot=FALSE, method='MBD')",
    input = list(grid = grid, Y = Y),
    output = list(depth = fb$depth, outpoint = fb$outpoint, medcurve = as.integer(fb$medcurve))
  )
}

# ---- fdepth ---------------------------------------------------------------
#
# fdepth() requires `data` as list(y = <matrix>): a bare matrix fails with
# "$ operator is invalid for atomic vectors".

fdepth_case <- function(w) {
  grid <- seq(0, 365, length.out = 365)
  Y <- eval.fd(grid, w$fdobj)
  fdp <- fdepth(list(argvals = grid, y = Y))
  add_case(
    name = "fdepth_weather_temperature_fm",
    r_call = "fdepth(list(argvals=grid, y=Y))  # type='FM' (default)",
    input = list(grid = grid, Y = Y),
    output = list(
      median = fdp$output$median,
      lmed = fdp$output$lmed,
      mtrim = fdp$output$mtrim,
      ltrim = fdp$output$ltrim,
      prof = fdp$output$prof
    )
  )
}

# ---- tperm.fd -------------------------------------------------------------

tperm_case <- function(w) {
  region <- CanadianWeather$region
  atl <- w$fdobj[region == "Atlantic"]
  pac <- w$fdobj[region == "Pacific"]
  set.seed(42)
  tp <- tperm.fd(atl, pac, nperm = 200, q = 0.95, plotres = FALSE)
  add_case(
    name = "tperm_fd_weather_atlantic_vs_pacific",
    r_call = "set.seed(42); tperm.fd(temp_fd[region=='Atlantic'], temp_fd[region=='Pacific'], nperm=200, q=0.95, plotres=FALSE)",
    input = list(seed = 42, nperm = 200, q = 0.95),
    output = list(
      pval = tp$pval,
      qval = tp$qval,
      Tobs = tp$Tobs,
      Tnull = tp$Tnull,
      argvals = tp$argvals
    )
  )
}

# ---- Fperm.fd on a small fRegress setup (weather temp ~ region dummy) -----

fperm_case <- function(w) {
  region <- CanadianWeather$region
  n <- length(region)
  xfdlist <- list(const = rep(1, n), atlantic = as.numeric(region == "Atlantic"))
  betalist <- list(
    const = fdPar(w$basis, 2, 1e2),
    atlantic = fdPar(w$basis, 2, 1e2)
  )
  set.seed(11)
  fp <- Fperm.fd(w$fdobj, xfdlist, betalist, nperm = 200, q = 0.95, plotres = FALSE)
  add_case(
    name = "fperm_fd_weather_temp_atlantic_dummy",
    r_call = "set.seed(11); Fperm.fd(temp_fd, list(const=rep(1,35), atlantic=as.numeric(region=='Atlantic')), list(const=fdPar(basis,2,1e2), atlantic=fdPar(basis,2,1e2)), nperm=200, q=0.95, plotres=FALSE)",
    input = list(seed = 11, nperm = 200, q = 0.95, atlantic = as.numeric(region == "Atlantic")),
    output = list(
      pval = fp$pval,
      qval = fp$qval,
      Fobs = fp$Fobs,
      Fnull = fp$Fnull,
      argvals = fp$argvals
    )
  )
}

# ---- mean.fd / sd.fd on the weather dataset (real-dataset duplicate of ----
# ---- the synthetic core.json cases) ----------------------------------------

mean_sd_case <- function(w) {
  grid <- seq(0, 365, length.out = 50)
  mfd <- mean.fd(w$fdobj)
  sfd <- sd.fd(w$fdobj)
  add_case(
    name = "mean_sd_fd_weather_temperature",
    r_call = "mean.fd(temp_fd); sd.fd(temp_fd)",
    input = list(grid = grid),
    output = list(
      mean_coefs = mfd$coefs,
      mean_fitted = eval.fd(grid, mfd),
      sd_coefs = sfd$coefs,
      sd_fitted = eval.fd(grid, sfd)
    )
  )
}

w <- weather_basis_fd()
var_cor_cases(w)
fbplot_case(w)
fdepth_case(w)
tperm_case(w)
fperm_case(w)
mean_sd_case(w)

finalize("stats", rtol = 1e-8)
