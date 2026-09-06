# tools/golden_r/smoothing.R
#
# Golden cases for the `smoothing` module: smooth.basis / smooth.basisPar /
# Data2fd / smooth.monotone / smooth.pos, plus the lambda<->df/gcv helper
# functions (lambda2gcv, lambda2df, df2lambda). Driven by tools/make_golden.py,
# which sources common.R first (defines add_case, eval_points, basis_info,
# finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Uses real fda datasets (CanadianWeather, growth) plus seeded synthetic
# data built from the RNG stream seeded once in common.R.
#
# Clean-room note: only the public fda API is called (create.*.basis,
# fdPar, vec2Lfd, smooth.basis, smooth.basisPar, Data2fd, smooth.monotone,
# smooth.pos, eval.fd, eval.monfd, eval.posfd, lambda2gcv, lambda2df,
# df2lambda). No fda source is read or copied.
#
# fda API quirks discovered while writing this script (see PROGRESS.md):
#   - lambda2gcv() takes log10(lambda); lambda2df()/df2lambda() take linear
#     lambda/df. Asymmetric across three closely related functions.
#   - Data2fd()'s 4th positional argument is `nderiv`, not `lambda` (despite
#     the task shorthand "Data2fd(argvals, y, basis, lambda)") -- lambda
#     must be passed by name.
#   - smooth.pos() has no `beta` in its return value (unlike smooth.monotone).

# ---- smooth.basis family: shared output recorder -----------------------

#' Record the standard smooth.basis/smooth.basisPar output fields.
sb_output <- function(sm, argvals) {
  list(
    coefs = sm$fd$coefs,
    coefs_dim = dim(sm$fd$coefs),
    df = sm$df,
    gcv = sm$gcv,
    SSE = sm$SSE,
    penmat = sm$penmat,
    y2cMap = sm$y2cMap,
    fitted = eval.fd(argvals, sm$fd)
  )
}

# ---- weather Fourier harmonic-accelerator smoothing ---------------------

weather_fourier_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  Lfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  lambdas <- c(1e-2, 1e2, 1e4, 1e6)
  for (lambda in lambdas) {
    fdParobj <- fdPar(basis, Lfd, lambda)
    sm <- smooth.basis(day, temp, fdParobj)
    add_case(
      name = sprintf("smooth_basis_weather_fourier65_harmonic_lambda%g", lambda),
      r_call = sprintf(
        "smooth.basis(1:365, CanadianWeather$dailyAv[,,'Temperature.C'], fdPar(create.fourier.basis(c(0,365),65), vec2Lfd(c(0,(2*pi/365)^2,0),c(0,365)), %g))",
        lambda
      ),
      input = list(argvals = day, basis = basis_info(basis), lambda = lambda),
      output = sb_output(sm, day)
    )
  }
}

# ---- growth hgtm: B-spline order 6, knots at ages, Lfd 4 ----------------

growth_hgtm_cases <- function() {
  age <- growth$age
  basis <- create.bspline.basis(range(age), breaks = age, norder = 6)
  y <- growth$hgtm
  lambdas <- c(1e-2, 1, 1e2)
  for (lambda in lambdas) {
    fdParobj <- fdPar(basis, 4, lambda)
    sm <- smooth.basis(age, y, fdParobj)
    add_case(
      name = sprintf("smooth_basis_growth_hgtm_bspline6_lfd4_lambda%g", lambda),
      r_call = sprintf(
        "smooth.basis(growth$age, growth$hgtm, fdPar(create.bspline.basis(range(age),breaks=age,norder=6), 4, %g))",
        lambda
      ),
      input = list(argvals = age, basis = basis_info(basis), lambda = lambda),
      output = sb_output(sm, age)
    )
  }
}

# ---- synthetic B-spline order 4, 12 basis, 1 and 3 curves ---------------

synthetic_bspline_cases <- function() {
  dom <- c(0, 10)
  argvals <- seq(dom[1], dom[2], length.out = 50)
  basis <- create.bspline.basis(dom, nbasis = 12, norder = 4)
  true_coef <- c(2, 1, -1, 3, 0, -2, 2, 1, -1, 0, 1, -1)
  true_vals <- eval.fd(argvals, fd(true_coef, basis))
  lambdas <- c(0, 1e-3, 1)
  for (ncurve in c(1, 3)) {
    noise <- matrix(rnorm(length(argvals) * ncurve, sd = 0.3), nrow = length(argvals), ncol = ncurve)
    y <- matrix(true_vals, nrow = length(argvals), ncol = ncurve) + noise
    for (lambda in lambdas) {
      fdParobj <- fdPar(basis, 2, lambda)
      sm <- smooth.basis(argvals, y, fdParobj)
      add_case(
        name = sprintf("smooth_basis_synthetic_bspline4_n12_ncurve%d_lambda%g", ncurve, lambda),
        r_call = sprintf(
          "smooth.basis(argvals, y[,1:%d], fdPar(create.bspline.basis(c(0,10),12,4), 2, %g))",
          ncurve, lambda
        ),
        input = list(argvals = argvals, basis = basis_info(basis), lambda = lambda, y = y),
        output = sb_output(sm, argvals)
      )
    }
  }
}

# ---- weight vector case --------------------------------------------------

wtvec_case <- function() {
  dom <- c(0, 10)
  argvals <- seq(dom[1], dom[2], length.out = 30)
  basis <- create.bspline.basis(dom, nbasis = 12, norder = 4)
  true <- sin(argvals) + 0.1 * argvals
  y <- true + rnorm(length(argvals), sd = 0.2)
  wtvec <- seq(0.5, 2, length.out = length(argvals))
  fdParobj <- fdPar(basis, 2, 1e-2)
  sm <- smooth.basis(argvals, y, fdParobj, wtvec = wtvec)
  add_case(
    name = "smooth_basis_synthetic_wtvec",
    r_call = "smooth.basis(argvals, y, fdPar(basis, 2, 1e-2), wtvec=wtvec)",
    input = list(argvals = argvals, basis = basis_info(basis), lambda = 1e-2, y = y, wtvec = wtvec),
    output = sb_output(sm, argvals)
  )
}

# ---- Lfdobj as bare integer 0 -------------------------------------------

lfd_integer_case <- function() {
  dom <- c(0, 10)
  argvals <- seq(dom[1], dom[2], length.out = 30)
  basis <- create.bspline.basis(dom, nbasis = 10, norder = 4)
  y <- cos(argvals) + rnorm(length(argvals), sd = 0.15)
  fdParobj <- fdPar(basis, 0, 1e-2)
  sm <- smooth.basis(argvals, y, fdParobj)
  add_case(
    name = "smooth_basis_synthetic_lfd_integer0",
    r_call = "smooth.basis(argvals, y, fdPar(basis, 0, 1e-2))",
    input = list(argvals = argvals, basis = basis_info(basis), lambda = 1e-2, y = y),
    output = sb_output(sm, argvals)
  )
}

# ---- 2-D y: bare vector (single curve, no matrix dim) --------------------

y_vector_case <- function() {
  dom <- c(0, 1)
  argvals <- eval_points(dom)
  basis <- create.bspline.basis(dom, nbasis = 8, norder = 4)
  true_coef <- c(1, -1, 2, 0, -2, 1, 1, -1)
  true_vals <- eval.fd(argvals, fd(true_coef, basis))
  y <- as.vector(true_vals) + rnorm(length(argvals), sd = 0.1)
  fdParobj <- fdPar(basis, 2, 1e-3)
  sm <- smooth.basis(argvals, y, fdParobj)
  add_case(
    name = "smooth_basis_synthetic_y_bare_vector",
    r_call = "smooth.basis(argvals, y, fdPar(basis, 2, 1e-3))  # y is a bare vector, not a matrix",
    input = list(argvals = argvals, basis = basis_info(basis), lambda = 1e-3, y = y),
    output = sb_output(sm, argvals)
  )
}

# ---- 3-D y: multivariate (handwrit: 1401 x 20 x 2) -----------------------

y_3d_case <- function() {
  basis <- create.bspline.basis(range(handwritTime), nbasis = 21, norder = 6)
  fdParobj <- fdPar(basis, 2, 1e-4)
  sm <- smooth.basis(handwritTime, handwrit, fdParobj)
  add_case(
    name = "smooth_basis_handwrit_3d_y",
    r_call = "smooth.basis(handwritTime, handwrit, fdPar(create.bspline.basis(range(handwritTime),21,6), 2, 1e-4))",
    input = list(argvals = handwritTime, basis = basis_info(basis), lambda = 1e-4, y_dim = dim(handwrit)),
    output = sb_output(sm, handwritTime)
  )
}

# ---- lambda2gcv / lambda2df / df2lambda + GCV grid -----------------------

lambda_functions_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  Lfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  fdParobj <- fdPar(basis, Lfd, lambda = 1)

  # lambda2gcv: log10-scale lambda -> gcv vector (one per curve)
  log10lambdas <- c(-2, 0, 2, 4, 6)
  for (ll in log10lambdas) {
    gcv <- lambda2gcv(ll, day, temp, fdParobj)
    add_case(
      name = sprintf("lambda2gcv_weather_log10lambda_%s", gsub("-", "m", as.character(ll))),
      r_call = sprintf("lambda2gcv(%g, 1:365, temp, fdParobj)", ll),
      input = list(log10lambda = ll),
      output = list(gcv = gcv)
    )
  }

  # lambda2df: linear-scale lambda -> scalar df
  lambdas_lin <- c(1e-2, 1, 1e2)
  for (lam in lambdas_lin) {
    df <- lambda2df(day, basis, wtvec = rep(1, 365), Lfdobj = Lfd, lambda = lam)
    add_case(
      name = sprintf("lambda2df_weather_lambda%g", lam),
      r_call = sprintf("lambda2df(1:365, basis, wtvec=rep(1,365), Lfdobj=Lfd, lambda=%g)", lam),
      input = list(lambda = lam),
      output = list(df = df)
    )
  }

  # df2lambda: linear-scale df -> scalar lambda
  dfs <- c(10, 20, 30)
  for (df_target in dfs) {
    lam <- df2lambda(day, basis, wtvec = rep(1, 365), Lfdobj = Lfd, df = df_target)
    add_case(
      name = sprintf("df2lambda_weather_df%g", df_target),
      r_call = sprintf("df2lambda(1:365, basis, wtvec=rep(1,365), Lfdobj=Lfd, df=%g)", df_target),
      input = list(df = df_target),
      output = list(lambda = lam)
    )
  }

  # GCV grid: lambdas 10^seq(-4, 8, by=0.25) -> mean gcv per lambda
  log10grid <- seq(-4, 8, by = 0.25)
  gcv_matrix <- sapply(log10grid, function(ll) lambda2gcv(ll, day, temp, fdParobj))
  mean_gcv <- colMeans(gcv_matrix)
  add_case(
    name = "lambda2gcv_weather_grid_mean",
    r_call = "sapply(10^seq(-4,8,by=0.25) |> log10(), function(ll) mean(lambda2gcv(ll, 1:365, temp, fdParobj)))",
    input = list(log10lambda_grid = log10grid),
    output = list(mean_gcv = mean_gcv, argmin_log10lambda = log10grid[which.min(mean_gcv)])
  )
}

# ---- smooth.basisPar ------------------------------------------------------

smooth_basisPar_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis1 <- create.fourier.basis(c(0, 365), nbasis = 65)
  sm1 <- smooth.basisPar(day, temp, fdobj = basis1, Lfdobj = 2, lambda = 1e2)
  add_case(
    name = "smooth_basisPar_weather_fourier",
    r_call = "smooth.basisPar(1:365, temp, fdobj=basis, Lfdobj=2, lambda=1e2)",
    input = list(argvals = day, basis = basis_info(basis1), lambda = 1e2),
    output = sb_output(sm1, day)
  )

  age <- growth$age
  basis2 <- create.bspline.basis(range(age), breaks = age, norder = 6)
  sm2 <- smooth.basisPar(age, growth$hgtm, fdobj = basis2, Lfdobj = 4, lambda = 1)
  add_case(
    name = "smooth_basisPar_growth_hgtm",
    r_call = "smooth.basisPar(growth$age, growth$hgtm, fdobj=basis, Lfdobj=4, lambda=1)",
    input = list(argvals = age, basis = basis_info(basis2), lambda = 1),
    output = sb_output(sm2, age)
  )

  dom <- c(0, 10)
  argvals <- seq(dom[1], dom[2], length.out = 40)
  basis3 <- create.bspline.basis(dom, nbasis = 12, norder = 4)
  true <- sin(argvals * 0.5)
  y <- true + rnorm(length(argvals), sd = 0.2)
  sm3 <- smooth.basisPar(argvals, y, fdobj = basis3, Lfdobj = 2, lambda = 1e-2)
  add_case(
    name = "smooth_basisPar_synthetic",
    r_call = "smooth.basisPar(argvals, y, fdobj=basis, Lfdobj=2, lambda=1e-2)",
    input = list(argvals = argvals, basis = basis_info(basis3), lambda = 1e-2, y = y),
    output = sb_output(sm3, argvals)
  )
}

# ---- Data2fd ---------------------------------------------------------------

data2fd_cases <- function() {
  dom <- c(0, 1)
  argvals <- eval_points(dom)
  basis1 <- create.bspline.basis(dom, nbasis = 10, norder = 4)
  y1 <- sin(2 * pi * argvals) + rnorm(length(argvals), sd = 0.1)
  fdobj1 <- Data2fd(argvals, y1, basis1, lambda = 1e-3)
  add_case(
    name = "data2fd_synthetic_bspline",
    r_call = "Data2fd(argvals, y, basis, lambda=1e-3)  # lambda passed by name: Data2fd's 4th positional arg is nderiv",
    input = list(argvals = argvals, basis = basis_info(basis1), lambda = 1e-3, y = y1),
    output = list(coefs = fdobj1$coefs)
  )

  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis2 <- create.fourier.basis(c(0, 365), nbasis = 65)
  fdobj2 <- Data2fd(day, temp, basis2, lambda = 1e2)
  add_case(
    name = "data2fd_weather_fourier",
    r_call = "Data2fd(1:365, temp, basis, lambda=1e2)",
    input = list(argvals = day, basis = basis_info(basis2), lambda = 1e2),
    output = list(coefs = fdobj2$coefs)
  )

  age <- growth$age
  basis3 <- create.bspline.basis(range(age), breaks = age, norder = 6)
  fdobj3 <- Data2fd(age, growth$hgtf, basis3, lambda = 1)
  add_case(
    name = "data2fd_growth_hgtf",
    r_call = "Data2fd(growth$age, growth$hgtf, basis, lambda=1)",
    input = list(argvals = age, basis = basis_info(basis3), lambda = 1),
    output = list(coefs = fdobj3$coefs)
  )
}

# ---- smooth.monotone -------------------------------------------------------
#
# Design note: the task spec asks for growth hgtf "B-spline order 6, 13
# basis (knots at ages)" -- but create.bspline.basis(breaks=growth$age,
# norder=6) with all 31 raw ages gives nbasis=35 (nbreaks + norder - 2), not
# 13. To hit nbasis=13 as specified, we use a reduced 9-point break sequence
# spanning the same range (own design choice, not derived from any source).

smooth_monotone_cases <- function() {
  age <- growth$age
  breaks9 <- seq(min(age), max(age), length.out = 9)
  basis <- create.bspline.basis(range(age), breaks = breaks9, norder = 6)
  y <- growth$hgtf[, 1:3]
  Wfd0 <- fd(matrix(0, basis$nbasis, ncol(y)), basis)
  WfdParobj <- fdPar(Wfd0, 3, 10^(-0.5))
  sm <- smooth.monotone(age, y, WfdParobj, conv = 1e-10, iterlim = 100, dbglev = 0)
  monvals <- eval.monfd(age, sm$Wfdobj)
  fitted <- sweep(monvals, 2, sm$beta[2, ], `*`)
  fitted <- sweep(fitted, 2, sm$beta[1, ], `+`)
  deriv1 <- eval.monfd(age, sm$Wfdobj, 1)
  add_case(
    name = "smooth_monotone_growth_hgtf_3girls",
    r_call = "smooth.monotone(growth$age, growth$hgtf[,1:3], fdPar(fd(matrix(0,13,3),basis), 3, 10^-0.5), conv=1e-10, iterlim=100)",
    input = list(argvals = age, basis = basis_info(basis), lambda = 10^(-0.5), y = y),
    output = list(
      Wfdobj_coefs = sm$Wfdobj$coefs,
      beta = sm$beta,
      fitted = fitted,
      deriv1 = deriv1
    ),
    rtol = 1e-5
  )

  t <- seq(0, 1, length.out = 41)
  true_mono <- 2 + 3 / (1 + exp(-8 * (t - 0.5)))
  y2 <- true_mono + rnorm(length(t), sd = 0.05)
  basis2 <- create.bspline.basis(c(0, 1), nbasis = 10, norder = 6)
  Wfd0b <- fd(matrix(0, basis2$nbasis, 1), basis2)
  WfdParobj2 <- fdPar(Wfd0b, 2, 1e-4)
  sm2 <- smooth.monotone(t, y2, WfdParobj2, conv = 1e-10, iterlim = 100, dbglev = 0)
  fitted2 <- sm2$beta[1, 1] + sm2$beta[2, 1] * eval.monfd(t, sm2$Wfdobj)
  deriv1b <- eval.monfd(t, sm2$Wfdobj, 1)
  add_case(
    name = "smooth_monotone_synthetic_sigmoid",
    r_call = "smooth.monotone(t, y, fdPar(fd(matrix(0,10,1),basis), 2, 1e-4), conv=1e-10, iterlim=100)",
    input = list(argvals = t, basis = basis_info(basis2), lambda = 1e-4, y = y2),
    output = list(
      Wfdobj_coefs = sm2$Wfdobj$coefs,
      beta = sm2$beta,
      fitted = fitted2,
      deriv1 = deriv1b
    ),
    rtol = 1e-5
  )
}

# ---- smooth.pos -------------------------------------------------------------

smooth_pos_cases <- function() {
  t <- seq(0, 1, length.out = 41)
  basis_true <- create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
  truecoef <- matrix(rnorm(8 * 2, sd = 0.5), nrow = 8, ncol = 2)
  truefd <- fd(truecoef, basis_true)
  y <- exp(eval.fd(t, truefd)) + matrix(rnorm(length(t) * 2, sd = 0.05), nrow = length(t), ncol = 2)
  y <- pmax(y, 0.01)
  basis <- create.bspline.basis(c(0, 1), nbasis = 10, norder = 6)
  Wfd0 <- fd(matrix(0, basis$nbasis, 2), basis)
  WfdParobj <- fdPar(Wfd0, 2, 1e-2)
  sm <- smooth.pos(t, y, WfdParobj, conv = 1e-10, iterlim = 100, dbglev = 0)
  fitted <- eval.posfd(t, sm$Wfdobj)
  add_case(
    name = "smooth_pos_synthetic_2curves",
    r_call = "smooth.pos(t, y, fdPar(fd(matrix(0,10,2),basis), 2, 1e-2), conv=1e-10, iterlim=100)",
    input = list(argvals = t, basis = basis_info(basis), lambda = 1e-2, y = y),
    output = list(Wfdobj_coefs = sm$Wfdobj$coefs, fitted = fitted),
    rtol = 1e-5
  )

  day <- 1:365
  precip <- CanadianWeather$dailyAv[, , "Precipitation.mm"]
  place_names <- dimnames(CanadianWeather$dailyAv)[[2]]
  y2 <- precip[, place_names == "Pr. Rupert"]
  basis2 <- create.fourier.basis(c(0, 365), nbasis = 15)
  Wfd0b <- fd(matrix(0, basis2$nbasis, 1), basis2)
  WfdParobj2 <- fdPar(Wfd0b, 2, 1e2)
  sm2 <- smooth.pos(day, y2, WfdParobj2, conv = 1e-10, iterlim = 100, dbglev = 0)
  fitted2 <- eval.posfd(day, sm2$Wfdobj)
  add_case(
    name = "smooth_pos_weather_precip_pr_rupert",
    r_call = "smooth.pos(1:365, precip[,'Pr. Rupert'], fdPar(fd(matrix(0,15,1),basis), 2, 1e2), conv=1e-10, iterlim=100)",
    input = list(argvals = day, basis = basis_info(basis2), lambda = 1e2, y = y2),
    output = list(Wfdobj_coefs = sm2$Wfdobj$coefs, fitted = fitted2),
    rtol = 1e-5
  )
}

# ---- irregular argvals: 3 separate single-curve fits ----------------------
#
# fda's smooth.basis() has no direct support for per-curve irregular argvals
# in one call (all curves in a `y` matrix share one `argvals` vector); we
# record 3 independent single-curve fits, each with its own argvals, so the
# Python side can exercise its irregular-argvals path against them one at a
# time.

irregular_argvals_cases <- function() {
  basis <- create.bspline.basis(c(0, 10), nbasis = 8, norder = 4)
  ns <- c(20, 15, 25)
  for (i in seq_along(ns)) {
    n <- ns[i]
    argvals <- sort(runif(n, 0, 10))
    true <- sin(argvals * 0.5) + 0.05 * argvals^2
    y <- true + rnorm(n, sd = 0.1)
    fdParobj <- fdPar(basis, 2, 1e-2)
    sm <- smooth.basis(argvals, y, fdParobj)
    add_case(
      name = sprintf("smooth_basis_irregular_argvals_case%d_n%d", i, n),
      r_call = "smooth.basis(argvals, y, fdPar(basis, 2, 1e-2))  # per-case distinct argvals",
      input = list(argvals = argvals, basis = basis_info(basis), lambda = 1e-2, y = y),
      output = sb_output(sm, argvals)
    )
  }
}

weather_fourier_cases()
growth_hgtm_cases()
synthetic_bspline_cases()
wtvec_case()
lfd_integer_case()
y_vector_case()
y_3d_case()
lambda_functions_cases()
smooth_basisPar_cases()
data2fd_cases()
smooth_monotone_cases()
smooth_pos_cases()
irregular_argvals_cases()

finalize("smoothing", rtol = 1e-8)
