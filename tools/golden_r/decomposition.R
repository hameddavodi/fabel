# tools/golden_r/decomposition.R
#
# Golden cases for the `decomposition` module: pca.fd / varmx.pca.fd / cca.fd.
# Driven by tools/make_golden.py, which sources common.R first (defines
# add_case, eval_points, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Uses real fda datasets (CanadianWeather, growth, gait) plus seeded
# synthetic data built from the RNG stream seeded once in common.R.
#
# Clean-room note: only the public fda API is called (create.*.basis,
# fdPar, vec2Lfd, smooth.basis, pca.fd, varmx.pca.fd, cca.fd, fd). No fda
# source is read or copied.
#
# fda API quirks / notes discovered while writing this script:
#   - pca.fd()$values is the FULL eigenvalue spectrum (length = nbasis, or
#     nbasis*nvar for a multivariate fd), not truncated to nharm -- only
#     `harmonics`/`scores`/`varprop` are truncated to nharm.
#   - cca.fd()$ccacorr is likewise the full-length canonical-correlation
#     spectrum (length = nbasis of the underlying basis), not truncated to
#     ncan -- only ccawtfd1/ccawtfd2/ccavar1/ccavar2 are truncated to ncan.
#   - Eigenvector/canonical-weight sign is arbitrary (depends on LAPACK's
#     internal sign convention); recorded as R gives it. Python parity
#     tests must sign-align columns before comparing harmonics/scores/
#     ccawtfd/ccavar, not try to reproduce R's sign exactly.

#' Record a pca.fd (or varmx.pca.fd) result.
pca_output <- function(pcafd) {
  out <- list(
    values = pcafd$values,
    harmonics_coefs = pcafd$harmonics$coefs,
    harmonics_coefs_dim = dim(pcafd$harmonics$coefs),
    scores = pcafd$scores,
    scores_dim = dim(pcafd$scores),
    varprop = pcafd$varprop,
    meanfd_coefs = pcafd$meanfd$coefs,
    meanfd_coefs_dim = dim(pcafd$meanfd$coefs)
  )
  if (!is.null(pcafd$rotmat)) {
    out$rotmat <- pcafd$rotmat
  }
  out
}

# ---- weather temperature: Fourier(65) harmonic-accelerator smoothed ------

weather_pca_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  temp_fd <- smooth.basis(day, temp, fdPar(basis, harmaccelLfd, 1e2))$fd

  varmx_targets <- list()
  for (nharm in c(2, 4)) {
    for (harm_lambda in c(0, 1e4)) {
      harmfdPar <- fdPar(basis, harmaccelLfd, harm_lambda)
      pcafd <- pca.fd(temp_fd, nharm = nharm, harmfdPar = harmfdPar)
      case_name <- sprintf("pca_fd_weather_nharm%d_harmlambda%g", nharm, harm_lambda)
      add_case(
        name = case_name,
        r_call = sprintf(
          "pca.fd(temp_fd, nharm=%d, harmfdPar=fdPar(basis, harmaccelLfd, %g))",
          nharm, harm_lambda
        ),
        input = list(nharm = nharm, harm_lambda = harm_lambda, basis = basis_info(basis)),
        output = pca_output(pcafd)
      )
      if (nharm == 4) {
        varmx_targets[[length(varmx_targets) + 1]] <- list(name = case_name, pcafd = pcafd)
      }
    }
  }

  # varmx.pca.fd on the nharm=4 weather cases (both harmfdPar lambdas)
  for (vt in varmx_targets) {
    vm <- varmx.pca.fd(vt$pcafd)
    add_case(
      name = sprintf("varmx_pca_fd_%s", vt$name),
      r_call = sprintf("varmx.pca.fd(pcafd)  # pcafd = %s", vt$name),
      input = list(source_case = vt$name),
      output = pca_output(vm)
    )
  }
}

# ---- growth hgtm: B-spline order 6, knots at ages, Lfd 4 -----------------

growth_pca_case <- function() {
  age <- growth$age
  basis <- create.bspline.basis(range(age), breaks = age, norder = 6)
  hgtm_fd <- smooth.basis(age, growth$hgtm, fdPar(basis, 4, 1))$fd
  pcafd <- pca.fd(hgtm_fd, nharm = 3)
  add_case(
    name = "pca_fd_growth_hgtm_nharm3",
    r_call = "pca.fd(hgtm_fd, nharm=3)  # hgtm_fd: bspline6 knots-at-ages, Lfd4, lambda=1",
    input = list(nharm = 3, basis = basis_info(basis)),
    output = pca_output(pcafd)
  )
  vm <- varmx.pca.fd(pcafd)
  add_case(
    name = "varmx_pca_fd_growth_hgtm_nharm3",
    r_call = "varmx.pca.fd(pcafd)  # pcafd = pca_fd_growth_hgtm_nharm3",
    input = list(source_case = "pca_fd_growth_hgtm_nharm3"),
    output = pca_output(vm)
  )
}

# ---- synthetic 10-curve B-spline fd ---------------------------------------

synthetic_pca_case <- function() {
  dom <- c(0, 10)
  basis <- create.bspline.basis(dom, nbasis = 12, norder = 4)
  coefs <- matrix(rnorm(12 * 10), nrow = 12, ncol = 10)
  fdobj <- fd(coefs, basis)
  pcafd <- pca.fd(fdobj, nharm = 3)
  add_case(
    name = "pca_fd_synthetic_bspline_n10curves",
    r_call = "pca.fd(fdobj, nharm=3)  # fdobj: bspline4 12basis, 10 seeded random curves",
    input = list(nharm = 3, basis = basis_info(basis), coefs = unname(coefs)),
    output = pca_output(pcafd)
  )
}

# ---- multivariate: gait (3-D coefs, Hip + Knee angle) ---------------------

gait_pca_case <- function() {
  gaittime <- as.matrix(0:19) + 0.5
  gaitrange <- c(0, 20)
  gaitbasis <- create.fourier.basis(gaitrange, nbasis = 21)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 20)^2, 0), rangeval = gaitrange)
  gaitfd <- smooth.basis(gaittime, gait, fdPar(gaitbasis, harmaccelLfd, 1e-2))$fd
  pcafd <- pca.fd(gaitfd, nharm = 3)
  add_case(
    name = "pca_fd_gait_multivariate_nharm3",
    r_call = "pca.fd(gaitfd, nharm=3)  # gaitfd: fourier21 harmonic-accel smoothed, 3-D coefs (Hip+Knee)",
    input = list(nharm = 3, basis = basis_info(gaitbasis), y_dim = dim(gait)),
    output = pca_output(pcafd)
  )
}

# ---- cca.fd: weather temperature vs log precipitation ---------------------

weather_cca_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  precip <- CanadianWeather$dailyAv[, , "Precipitation.mm"]
  logprecip <- log(pmax(precip, 0.01))
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  temp_fd <- smooth.basis(day, temp, fdPar(basis, 2, 1e2))$fd
  logprecip_fd <- smooth.basis(day, logprecip, fdPar(basis, 2, 1e2))$fd

  for (lam in c(1e2, 1e6)) {
    ccafdPar <- fdPar(basis, 2, lam)
    ccafd <- cca.fd(temp_fd, logprecip_fd, ncan = 3, ccafdPar1 = ccafdPar, ccafdPar2 = ccafdPar)
    add_case(
      name = sprintf("cca_fd_weather_temp_vs_logprecip_lambda%g", lam),
      r_call = sprintf(
        "cca.fd(temp_fd, logprecip_fd, ncan=3, ccafdPar1=fdPar(basis,2,%g), ccafdPar2=fdPar(basis,2,%g))",
        lam, lam
      ),
      input = list(ncan = 3, lambda = lam, basis = basis_info(basis)),
      output = list(
        ccawtfd1_coefs = ccafd$ccawtfd1$coefs,
        ccawtfd2_coefs = ccafd$ccawtfd2$coefs,
        ccacorr = ccafd$ccacorr,
        ccavar1 = ccafd$ccavar1,
        ccavar2 = ccafd$ccavar2
      )
    )
  }
}

weather_pca_cases()
growth_pca_case()
synthetic_pca_case()
gait_pca_case()
weather_cca_cases()

finalize("decomposition", rtol = 1e-8)
