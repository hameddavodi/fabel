# tools/golden_r/bands.R
#
# Golden cases for pointwise confidence bands (fdatools.stats.confidence_band):
# the pointwise standard error of a smooth, phi(t)' S Sigma S' phi(t) with S
# the smooth.basis y2cMap, and the pointwise standard error of fRegress
# coefficient functions, theta(t)' Bvar_jj theta(t) with Bvar the coefficient
# covariance fRegress returns when given y2cMap and SigmaE.  All on
# CanadianWeather, as in Ramsay, Hooker & Graves (2009), sections 5.5 and 9.4.
# Driven by tools/make_golden.py, which sources common.R first.
#
# Clean-room note: only the public fda API is called (create.*.basis, fdPar,
# vec2Lfd, smooth.basis, eval.basis, eval.fd, fd, fRegress).  No fda source
# is read or copied.

day5 <- (1:365) - 0.5
rng <- c(0, 365)
tgrid <- eval_points(rng)
logprec <- CanadianWeather$dailyAv[, , "log10precip"]
nstation <- ncol(logprec)

#' Pointwise standard error sqrt(diag(A Sigma A')) of a linear map A of the data.
pointwise_se <- function(A, Sigma) {
  sqrt(rowSums((A %*% Sigma) * A))
}

# ---- (a) smooth of log10 precipitation: B-spline(53), D^2, lambda 1e4 ------

smooth_cases <- function() {
  basis <- create.bspline.basis(rng, nbasis = 53, norder = 4)
  lambda <- 1e4
  sm <- smooth.basis(day5, logprec, fdPar(basis, 2, lambda))
  y2c <- sm$y2cMap
  nobs <- length(day5)

  # (a1) Sigma = sigma^2 I with sigma^2 = SSE / (N (n - df)).
  sigma2 <- sm$SSE / (nstation * (nobs - sm$df))
  Sigma <- diag(rep(sigma2, nobs))
  phi0 <- eval.basis(tgrid, basis, 0)
  phi1 <- eval.basis(tgrid, basis, 1)
  se0 <- pointwise_se(phi0 %*% y2c, Sigma)
  se1 <- pointwise_se(phi1 %*% y2c, Sigma)
  z <- qnorm(0.975)
  fit1 <- as.numeric(eval.fd(tgrid, sm$fd[1]))
  add_case(
    name = "band_smooth_logprec_bspline53_sigma_scalar",
    r_call = paste(
      "sm <- smooth.basis(day5, logprec, fdPar(create.bspline.basis(c(0,365),53,4), 2, 1e4));",
      "sqrt(diag(Phi(t) %*% sm$y2cMap %*% (SSE/(N*(n-df)) I) %*% t(sm$y2cMap) %*% t(Phi(t))))"
    ),
    input = list(
      argvals = day5, basis = basis_info(basis), lambda = lambda, lfd = 2,
      t = tgrid, level = 0.95
    ),
    output = list(
      df = sm$df, SSE = sm$SSE, sigma2 = sigma2,
      stderr = se0, stderr_deriv1 = se1, z = z,
      lower_curve1 = fit1 - z * se0, upper_curve1 = fit1 + z * se0
    )
  )

  # (a2) Sigma = diag of the per-day residual variance across stations (the
  # book's pointwise error variance, unsmoothed).
  resid <- logprec - eval.fd(day5, sm$fd)
  varvec <- rowSums(resid^2) / (nstation - 1)
  se_pw <- pointwise_se(phi0 %*% y2c, diag(varvec))
  add_case(
    name = "band_smooth_logprec_bspline53_sigma_pointwise",
    r_call = "sqrt(diag(Phi(t) %*% y2cMap %*% diag(rowSums(resid^2)/(N-1)) %*% t(y2cMap) %*% t(Phi(t))))",
    input = list(
      argvals = day5, basis = basis_info(basis), lambda = lambda, lfd = 2,
      t = tgrid
    ),
    output = list(varvec = varvec, stderr = se_pw)
  )
}

# ---- (b) scalar response on scalar covariates (exact in R) ----------------

annualprec <- log10(apply(CanadianWeather$dailyAv[, , "Precipitation.mm"], 2, sum))

fregress_scalar_case <- function() {
  meantemp <- colMeans(CanadianWeather$dailyAv[, , "Temperature.C"])
  latitude <- as.numeric(CanadianWeather$coordinates[, 1])
  cbasis <- create.constant.basis(rng)
  xfdlist <- list(rep(1, nstation), meantemp, latitude)
  betalist <- list(fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0))
  fr <- fRegress(annualprec, xfdlist, betalist)
  resid <- annualprec - as.numeric(fr$yhatfdobj)
  sigma2 <- sum(resid^2) / (nstation - fr$df)
  fr_se <- fRegress(annualprec, xfdlist, betalist,
    y2cMap = diag(nstation), SigmaE = diag(rep(sigma2, nstation))
  )
  beta <- sapply(fr$betaestlist, function(b) as.numeric(b$fd$coefs))
  se <- sqrt(diag(fr_se$bvar))
  z <- qnorm(0.975)
  add_case(
    name = "band_fregress_logannualprec_on_meantemp_latitude",
    r_call = paste(
      "fRegress(annualprec, list(1, meantemp, latitude), constant betas,",
      "y2cMap=diag(35), SigmaE=diag(rep(sigma2,35)))$bvar"
    ),
    input = list(y = annualprec, meantemp = meantemp, latitude = latitude, level = 0.95),
    output = list(
      beta = beta, sigma2 = sigma2, bvar = unname(fr_se$bvar), stderr = se,
      lower = beta - z * se, upper = beta + z * se
    )
  )
}

# ---- (c) scalar response on the temperature curves (book section 9.4) -----

fregress_functional_case <- function() {
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.bspline.basis(rng, nbasis = 53, norder = 4)
  temp_fd <- smooth.basis(day5, temp, fdPar(basis, 2, 1e2))$fd
  cbasis <- create.constant.basis(rng)
  bbasis <- create.bspline.basis(rng, nbasis = 7, norder = 4)
  const_fd <- fd(matrix(rep(1, nstation), nrow = 1), cbasis)
  xfdlist <- list(const_fd, temp_fd)
  lambda <- 1e5
  betalist <- list(fdPar(cbasis, 0, 0), fdPar(bbasis, 2, lambda))
  fr <- fRegress(annualprec, xfdlist, betalist)
  resid <- annualprec - as.numeric(fr$yhatfdobj)
  sigma2 <- sum(resid^2) / (nstation - fr$df)
  fr_se <- fRegress(annualprec, xfdlist, betalist,
    y2cMap = diag(nstation), SigmaE = diag(rep(sigma2, nstation))
  )
  nb <- bbasis$nbasis
  idx <- 1 + seq_len(nb)
  theta <- eval.basis(tgrid, bbasis, 0)
  se <- pointwise_se(theta, fr_se$bvar[idx, idx])
  add_case(
    name = "band_fregress_logannualprec_on_temp_bspline7",
    r_call = paste(
      "fRegress(annualprec, list(const_fd, temp_fd), list(fdPar(cbasis,0,0),",
      "fdPar(create.bspline.basis(c(0,365),7,4),2,1e5)), y2cMap=diag(35),",
      "SigmaE=diag(rep(sigma2,35))); sqrt(diag(theta(t) Bvar_22 theta(t)'))"
    ),
    input = list(
      y = annualprec, temp_coefs = unname(temp_fd$coefs), temp_basis = basis_info(basis),
      beta_basis = basis_info(bbasis), lambda = lambda, t = tgrid
    ),
    output = list(
      beta1 = as.numeric(eval.fd(tgrid, fr$betaestlist[[2]]$fd)),
      df = fr$df, sigma2 = sigma2, bvar = unname(fr_se$bvar),
      stderr0 = sqrt(fr_se$bvar[1, 1]), stderr1 = se
    )
  )
}

smooth_cases()
fregress_scalar_case()
fregress_functional_case()
finalize("bands", rtol = 1e-8)
