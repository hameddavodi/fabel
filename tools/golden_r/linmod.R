# tools/golden_r/linmod.R
#
# Golden cases for `fdatools.regression.linmod`: R fda's `linmod`, the fully
# functional linear model
#
#     y_i(t) = beta0(t) + int x_i(s) beta1(s, t) ds + e_i(t)
#
# with an intercept function (fdPar) and a bivariate regression surface
# (bifdPar, separate roughness penalties in s and t).  Driven by
# tools/make_golden.py, which sources common.R first (defines add_case,
# basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis, fd,
# fdPar, bifd, bifdPar, vec2Lfd, smooth.basis, linmod, inprod,
# eval.penalty).  No fda source is read or copied.
#
# fda API quirks discovered while writing this script (see PROGRESS.md):
#   - linmod(..., wtvec = w) fails for every weight vector tried (plain
#     numeric, n x 1 matrix): "Error in tcrossprod(as.vector(X),
#     as.vector(Y)): requires numeric/complex matrix/vector arguments".
#     No case uses weights.
#   - There is no predict method for linmod objects.
#
# Besides linmod's own output, the synthetic cases record the integrals and
# penalty matrices R's inprod / eval.penalty return for the same bases
# (`r_integrals`).  They let the parity test separate the model (normal
# equations) from R's approximate quadrature.

#' R's inner-product and penalty matrices for one linmod design.
linmod_integrals <- function(xfd, yfd, abasis, sbasis, tbasis, La, Ls, Lt) {
  ybasis <- yfd$basis
  list(
    Z = unname(inprod(xfd, sbasis)),
    Gaa = unname(inprod(abasis, abasis)),
    Gat = unname(inprod(abasis, tbasis)),
    Gtt = unname(inprod(tbasis, tbasis)),
    Gss = unname(inprod(sbasis, sbasis)),
    Gay = unname(inprod(abasis, ybasis)),
    Gty = unname(inprod(tbasis, ybasis)),
    Ra = unname(eval.penalty(abasis, La)),
    Rs = unname(eval.penalty(sbasis, Ls)),
    Rt = unname(eval.penalty(tbasis, Lt))
  )
}

#' Fit linmod and return the case output.
linmod_fit <- function(xfd, yfd, abasis, sbasis, tbasis, La, Ls, Lt, lams) {
  betaList <- list(
    fdPar(abasis, La, lams[1]),
    bifdPar(bifd(matrix(0, sbasis$nbasis, tbasis$nbasis), sbasis, tbasis), Ls, Lt, lams[2], lams[3])
  )
  fit <- linmod(xfd, yfd, betaList)
  list(
    beta0estfd_coefs = unname(as.numeric(fit$beta0estfd$coefs)),
    beta1estbifd_coefs = unname(fit$beta1estbifd$coefs),
    yhatfdobj_coefs = unname(fit$yhatfdobj$coefs)
  )
}

# ---- (a) weather: log10 precipitation on temperature --------------------

weather_case <- function() {
  day <- 1:365
  harm <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  fb65 <- create.fourier.basis(c(0, 365), nbasis = 65)
  temp_fd <- smooth.basis(day, CanadianWeather$dailyAv[, , "Temperature.C"], fdPar(fb65, harm, 1e2))$fd
  prec_fd <- smooth.basis(day, CanadianWeather$dailyAv[, , "log10precip"], fdPar(fb65, harm, 1e2))$fd
  fb11 <- create.fourier.basis(c(0, 365), nbasis = 11)
  lams <- c(1e4, 1e6, 1e5)
  add_case(
    name = "linmod_weather_logprecip_on_temp",
    r_call = paste(
      "linmod(temp_fd, logprec_fd, list(fdPar(fb11, harm, 1e4),",
      "bifdPar(bifd(matrix(0,11,11), fb11, fb11), harm, harm, 1e6, 1e5)))",
      "# temp_fd, logprec_fd: smooth.basis(day, ., fdPar(fb65, harm, 1e2))$fd"
    ),
    input = list(
      x_basis = basis_info(fb65), y_basis = basis_info(fb65),
      x_coefs = unname(temp_fd$coefs), y_coefs = unname(prec_fd$coefs),
      alpha_basis = basis_info(fb11), s_basis = basis_info(fb11), t_basis = basis_info(fb11),
      penalty = "harmonic", period = 365, lambda_alpha = lams[1], lambda_s = lams[2], lambda_t = lams[3]
    ),
    output = linmod_fit(temp_fd, prec_fd, fb11, fb11, fb11, harm, harm, harm, lams)
  )
}

# ---- (b) synthetic: monomial bases, which R integrates exactly -----------
# Products of these polynomials have degree <= 4, so R's numerical inprod is
# exact to rounding.  s lives on [0, 1] and t on [0, 2]: the two continua differ.

monomial_case <- function(lams, name) {
  n <- 12
  xb <- create.monomial.basis(c(0, 1), 3)
  yb <- create.monomial.basis(c(0, 2), 3)
  ab <- create.monomial.basis(c(0, 2), 3)
  sb <- create.monomial.basis(c(0, 1), 2)
  tb <- create.monomial.basis(c(0, 2), 3)
  xc <- matrix(rnorm(3 * n), 3, n)
  yc <- matrix(rnorm(3 * n), 3, n)
  xfd <- fd(xc, xb)
  yfd <- fd(yc, yb)
  out <- linmod_fit(xfd, yfd, ab, sb, tb, 1, 1, 1, lams)
  out$r_integrals <- linmod_integrals(xfd, yfd, ab, sb, tb, 1, 1, 1)
  add_case(
    name = name,
    r_call = sprintf(
      "linmod(xfd, yfd, list(fdPar(ab, 1, %g), bifdPar(bifd(matrix(0,2,3), sb, tb), 1, 1, %g, %g)))",
      lams[1], lams[2], lams[3]
    ),
    input = list(
      x_basis = basis_info(xb), y_basis = basis_info(yb), x_coefs = unname(xc), y_coefs = unname(yc),
      alpha_basis = basis_info(ab), s_basis = basis_info(sb), t_basis = basis_info(tb),
      penalty = 1, lambda_alpha = lams[1], lambda_s = lams[2], lambda_t = lams[3]
    ),
    output = out
  )
}

# ---- (c) synthetic: cubic B-splines --------------------------------------

bspline_case <- function() {
  n <- 15
  xb <- create.bspline.basis(c(0, 1), 8, 4)
  yb <- create.bspline.basis(c(0, 1), 7, 4)
  sb <- create.bspline.basis(c(0, 1), 5, 4)
  tb <- create.bspline.basis(c(0, 1), 6, 4)
  ab <- create.bspline.basis(c(0, 1), 6, 4)
  xc <- matrix(rnorm(8 * n), 8, n)
  yc <- matrix(rnorm(7 * n), 7, n)
  xfd <- fd(xc, xb)
  yfd <- fd(yc, yb)
  lams <- c(1e-3, 1e-2, 1e-3)
  out <- linmod_fit(xfd, yfd, ab, sb, tb, 2, 2, 2, lams)
  out$r_integrals <- linmod_integrals(xfd, yfd, ab, sb, tb, 2, 2, 2)
  add_case(
    name = "linmod_synthetic_bspline",
    r_call = "linmod(xfd, yfd, list(fdPar(ab, 2, 1e-3), bifdPar(bifd(matrix(0,5,6), sb, tb), 2, 2, 1e-2, 1e-3)))",
    input = list(
      x_basis = basis_info(xb), y_basis = basis_info(yb), x_coefs = unname(xc), y_coefs = unname(yc),
      alpha_basis = basis_info(ab), s_basis = basis_info(sb), t_basis = basis_info(tb),
      penalty = 2, lambda_alpha = lams[1], lambda_s = lams[2], lambda_t = lams[3]
    ),
    output = out
  )
}

weather_case()
monomial_case(c(1e-2, 1e-1, 3e-2), "linmod_synthetic_monomial")
monomial_case(c(0, 0, 0), "linmod_synthetic_monomial_lambda0")
bspline_case()

finalize("linmod", rtol = 1e-8)
