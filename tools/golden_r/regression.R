# tools/golden_r/regression.R
#
# Golden cases for the `regression` module: fRegress and its satellites
# (predict.fRegress, fRegress.stderr, fRegress.CV) on the three book case
# studies (scalar-on-function, function-on-scalar, concurrent
# function-on-function) plus two synthetic sanity cases. Driven by
# tools/make_golden.py, which sources common.R first (defines add_case,
# eval_points, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# fdPar, vec2Lfd, smooth.basis, fRegress, predict.fRegress, fRegress.CV,
# fd). No fda source is read or copied.
#
# fda API quirks discovered while writing this script (see PROGRESS.md):
#   - fRegress.stderr(fRegressList, y2cMap, SigmaE), called standalone,
#     requires fRegressList$yfdobj to exist -- it errors ("'data' must be
#     of a vector type, was 'NULL'") for a SCALAR-response fRegress object,
#     because fRegress.double()'s return list has no $yfdobj slot at all
#     (only $yvec). The scalar-response path instead computes stderr by
#     passing y2cMap/SigmaE straight to fRegress(y, xfdlist, betalist,
#     y2cMap=, SigmaE=) -- fRegress.double() accepts these as named args
#     and populates betastderrlist/bvar/c2bMap internally. The standalone
#     fRegress.stderr() call works fine for a FUNCTIONAL-response
#     fRegress object (fRegress.fd() does set $yfdobj).
#   - fRegress.stderr()'s own return list uses different capitalization
#     for the functional-response path (betastderrlist, YhatStderr, Bvar,
#     c2bMap) vs the embedded scalar-response path inside fRegress()
#     (betastderrlist, bvar, c2bMap, no YhatStderr) -- Bvar/bvar are the
#     same quantity under different names depending on which code path
#     populated it.
#   - fRegress.stderr() / fRegress.CV() both require every xfdlist entry
#     to be an `fd` object (even "scalar" covariates) when the covariate
#     is a compile-time constant across curves: passing a bare numeric
#     vector for an intercept-like term makes fRegress.CV/stderr try to
#     eval.fd() a non-fd object and fail deep inside (NAs from
#     getbasismatrix, or "$ operator invalid for atomic vectors"-style
#     errors). fRegress() itself is fine with bare-vector scalar
#     covariates -- only these two satellite functions need the
#     constant-basis fd wrapping. Scalar covariates that vary per curve
#     (not all-constant, e.g. z1/z2 below) work as bare vectors
#     everywhere.
#   - A y value produced via inprod(fdobj, weightfd) is a 1-column
#     matrix, not a plain numeric vector -- fRegress.double() rejects it
#     ("Argument y is not a numeric vector"); must wrap with as.numeric().

#' Record a coefficient array from an fdPar/fd list entry generically.
beta_coefs <- function(betaestlist) {
  lapply(betaestlist, function(b) unname(b$fd$coefs))
}

#' Record the standard fRegress() output fields common to all cases.
fregress_output <- function(fr, functional_response) {
  out <- list(
    betaestlist_coefs = beta_coefs(fr$betaestlist),
    df = fr$df
  )
  if (functional_response) {
    out$yhatfdobj_coefs <- unname(fr$yhatfdobj$coefs)
  } else {
    out$yhatfdobj <- unname(as.numeric(fr$yhatfdobj))
  }
  if (!is.null(fr$GCV)) out$gcv <- fr$GCV
  if (!is.null(fr$OCV)) out$OCV <- fr$OCV
  if (!is.null(fr$Cmat)) out$Cmat <- unname(fr$Cmat)
  if (!is.null(fr$Dmat)) out$Dmat <- unname(fr$Dmat)
  out
}

# ---- (a) scalar response: log annual precipitation on temperature fd ----

scalar_response_case <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  temp_fd <- smooth.basis(day, temp, fdPar(basis, harmaccelLfd, 1e2))$fd

  annualprec <- log10(apply(CanadianWeather$dailyAv[, , "Precipitation.mm"], 2, sum))
  n <- length(annualprec)
  cbasis <- create.constant.basis(c(0, 365))
  fbasis5 <- create.fourier.basis(c(0, 365), nbasis = 5)
  betalist <- list(fdPar(cbasis, 0, 0), fdPar(fbasis5, harmaccelLfd, 1e2))
  const_fd <- fd(matrix(rep(1, n), nrow = 1), cbasis)
  xfdlist <- list(const_fd, temp_fd)

  fr <- fRegress(annualprec, xfdlist, betalist)
  add_case(
    name = "fregress_scalar_precip_on_temp",
    r_call = "fRegress(annualprec, list(const_fd, temp_fd), list(fdPar(cbasis,0,0), fdPar(fbasis5,harmaccelLfd,1e2)))",
    input = list(y = annualprec, basis1 = basis_info(cbasis), basis2 = basis_info(fbasis5), lambda = 1e2),
    output = fregress_output(fr, functional_response = FALSE)
  )

  pred <- predict(fr)
  add_case(
    name = "predict_fregress_scalar_precip_on_temp",
    r_call = "predict(fr)  # fr = fregress_scalar_precip_on_temp",
    input = list(source_case = "fregress_scalar_precip_on_temp"),
    output = list(predicted = unname(as.numeric(pred)))
  )

  # fRegress.stderr for the scalar-response case: pass y2cMap/SigmaE
  # directly to fRegress() (see header note) -- standalone
  # fRegress.stderr() errors here.
  resid <- annualprec - as.numeric(fr$yhatfdobj)
  sigma2 <- sum(resid^2) / (n - fr$df)
  y2cMap <- diag(n)
  SigmaE <- diag(rep(sigma2, n))
  fr_se <- fRegress(annualprec, xfdlist, betalist, y2cMap = y2cMap, SigmaE = SigmaE)
  add_case(
    name = "fregress_stderr_scalar_precip_on_temp",
    r_call = "fRegress(annualprec, xfdlist, betalist, y2cMap=diag(n), SigmaE=diag(rep(sigma2,n)))$betastderrlist",
    input = list(sigma2 = sigma2, source_case = "fregress_scalar_precip_on_temp"),
    output = list(betastderrlist_coefs = beta_coefs(fr_se$betastderrlist))
  )

  cv <- fRegress.CV(annualprec, xfdlist, betalist)
  add_case(
    name = "fregress_cv_scalar_precip_on_temp",
    r_call = "fRegress.CV(annualprec, xfdlist, betalist)",
    input = list(source_case = "fregress_scalar_precip_on_temp"),
    output = list(SSE.CV = cv$SSE.CV, errfd.cv = unname(as.numeric(cv$errfd.cv)))
  )
}

# ---- (b) functional response: temperature fd on climate-region dummies --

region_regression_case <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  sm <- smooth.basis(day, temp, fdPar(basis, harmaccelLfd, 1e2))
  temp_fd <- sm$fd
  y2cMap <- sm$y2cMap

  regions <- unique(CanadianWeather$region)
  xfdlist <- lapply(regions, function(r) as.numeric(CanadianWeather$region == r))
  names(xfdlist) <- regions
  betalist <- lapply(regions, function(r) fdPar(basis, harmaccelLfd, 1e2))

  fr <- fRegress(temp_fd, xfdlist, betalist)
  add_case(
    name = "fregress_functional_temp_on_region",
    r_call = "fRegress(temp_fd, xfdlist, betalist)  # xfdlist: 4 region-indicator dummies (no separate intercept)",
    input = list(regions = regions, basis = basis_info(basis), lambda = 1e2),
    output = fregress_output(fr, functional_response = TRUE)
  )

  fitted <- eval.fd(day, temp_fd)
  yhat <- eval.fd(day, fr$yhatfdobj)
  SigmaE <- cov(t(fitted - yhat))
  stderr <- fRegress.stderr(fr, y2cMap, SigmaE)
  add_case(
    name = "fregress_stderr_functional_temp_on_region",
    r_call = "fRegress.stderr(fr, y2cMap, SigmaE)  # fr = fregress_functional_temp_on_region; SigmaE = cov(t(residuals))",
    input = list(source_case = "fregress_functional_temp_on_region"),
    output = list(betastderrlist_coefs = beta_coefs(stderr$betastderrlist))
  )
}

# ---- (c) concurrent functional-on-functional: log precip on temperature -

concurrent_regression_case <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  precip <- CanadianWeather$dailyAv[, , "Precipitation.mm"]
  logprecip <- log(pmax(precip, 0.01))
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  temp_fd <- smooth.basis(day, temp, fdPar(basis, harmaccelLfd, 1e2))$fd
  logprecip_fd <- smooth.basis(day, logprecip, fdPar(basis, harmaccelLfd, 1e2))$fd

  n <- dim(temp)[2]
  fbasis5 <- create.fourier.basis(c(0, 365), nbasis = 5)
  betalist <- list(fdPar(fbasis5, 2, 1e2), fdPar(fbasis5, 2, 1e2))
  xfdlist <- list(rep(1, n), temp_fd)
  fr <- fRegress(logprecip_fd, xfdlist, betalist)
  add_case(
    name = "fregress_concurrent_logprecip_on_temp",
    r_call = "fRegress(logprecip_fd, list(rep(1,n), temp_fd), list(fdPar(fbasis5,2,1e2), fdPar(fbasis5,2,1e2)))",
    input = list(basis = basis_info(fbasis5), lambda = 1e2),
    output = fregress_output(fr, functional_response = TRUE)
  )
}

# ---- (d) synthetic: 2 scalar covariates + 1 functional covariate --------

synthetic_mixed_case <- function() {
  n <- 20
  sbasis <- create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
  xcoefs <- matrix(rnorm(8 * n), 8, n)
  xfd_s <- fd(xcoefs, sbasis)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  beta_true_fd <- fd(matrix(rnorm(8), 8, 1), sbasis)
  y <- 2 + 0.5 * z1 - 0.3 * z2 + as.numeric(inprod(xfd_s, beta_true_fd)) + rnorm(n, sd = 0.1)

  cbasis <- create.constant.basis(c(0, 1))
  xfdlist <- list(rep(1, n), z1, z2, xfd_s)
  betalist <- list(fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(sbasis, 2, 1e-2))
  fr <- fRegress(y, xfdlist, betalist)
  add_case(
    name = "fregress_synthetic_2scalar_1functional",
    r_call = "fRegress(y, list(rep(1,n), z1, z2, xfd_s), list(fdPar(cbasis,0,0) x3, fdPar(sbasis,2,1e-2)))",
    input = list(y = y, z1 = z1, z2 = z2, basis = basis_info(sbasis), coefs = unname(xcoefs), lambda = 1e-2),
    output = fregress_output(fr, functional_response = FALSE)
  )

  # ---- (e) same design, lambda=0: exact linear-algebra check -----------
  betalist0 <- list(fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(sbasis, 2, 0))
  fr0 <- fRegress(y, xfdlist, betalist0)
  add_case(
    name = "fregress_synthetic_2scalar_1functional_lambda0",
    r_call = "fRegress(y, xfdlist, betalist0)  # same design, lambda=0 on the functional term",
    input = list(y = y, z1 = z1, z2 = z2, basis = basis_info(sbasis), coefs = unname(xcoefs), lambda = 0),
    output = fregress_output(fr0, functional_response = FALSE)
  )
}

scalar_response_case()
region_regression_case()
concurrent_regression_case()
synthetic_mixed_case()

finalize("regression", rtol = 1e-8)
