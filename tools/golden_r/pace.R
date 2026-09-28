# tools/golden_r/pace.R
#
# Golden cases for `fabel.sparse` (sparse / longitudinal FPCA, PACE): R fda's
# smooth.sparse.mean, covPACE, pcaPACE and scoresPACE on seeded sparse
# subsamples of two real datasets (growth$hgtf and CanadianWeather daily
# temperature). Driven by tools/make_golden.py, which sources common.R first
# (defines add_case, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED and seeds
# the RNG).
#
# Clean-room note: only the public fda API is called. No fda source is read or
# copied; everything below was found by calling the functions as black boxes.
#
# fda API quirks found while writing this script (see PROGRESS notes in the
# PR / tests/parity/test_pace.py):
#   - Only the LIST data form works. The documented matrix-with-NA form gives
#     all-NA mean coefficients in smooth.sparse.mean and a "computationally
#     singular" error in covPACE. Each list element is a (n_i x 2) matrix:
#     column 1 the times, column 2 the values.
#   - smooth.sparse.mean(type = "bspline") ignores `nbasis`; the basis comes
#     from `knots` (break points, both ends included) and `norder`.
#   - covPACE needs `time` (any grid covering the data) for list data; with
#     time = NULL it fails inside eval.fd.
#   - scoresPACE fails ("evalarg contains 1 NA") unless the times are
#     integers, so it is only called on the day-numbered weather data. Many
#     curves get NA scores; see tests/parity/test_pace.py for the measured
#     defect.

#' Draw `m` sorted observation indices out of `n` for each of `n_curves`.
sparse_design <- function(n, m, n_curves) {
  lapply(seq_len(n_curves), function(i) sort(sample.int(n, m)))
}

#' Build the fda list data form from a dense (n_t x n_curves) matrix.
as_pace_list <- function(argvals, ymat, design) {
  lapply(seq_along(design), function(i) {
    k <- design[[i]]
    cbind(argvals[k], ymat[k, i])
  })
}

pace_input <- function(dl) {
  list(
    t = lapply(dl, function(m) as.numeric(m[, 1])),
    y = lapply(dl, function(m) as.numeric(m[, 2]))
  )
}

fd_output <- function(fdobj) {
  list(basis = basis_info(fdobj$basis), coefs = as.numeric(fdobj$coefs))
}

#' pcaPACE fields, plus R's inprod() cross-Gram matrices so that the parity
#' test can show which gaps come from inprod()'s Romberg quadrature alone.
pca_output <- function(pc, harm_basis, cov_basis) {
  list(
    harmonics = unname(pc$harmonics$coefs),
    values = as.numeric(pc$values),
    varprop = as.numeric(pc$varprop),
    inprod_harm_cov = unname(inprod(harm_basis, cov_basis)),
    inprod_cov_cov = unname(inprod(cov_basis, cov_basis))
  )
}

# ---- (a) growth heights of the 54 girls, 6 random ages each ----------------

growth_cases <- function() {
  age <- growth$age
  hgt <- growth$hgtf
  design <- sparse_design(length(age), 6, ncol(hgt))
  dl <- as_pace_list(age, hgt, design)
  data_in <- pace_input(dl)
  rng <- c(1, 18)
  knots <- seq(1, 18, length.out = 5)

  for (lam in c(0, 1)) {
    mfd <- smooth.sparse.mean(dl, NULL, rng = rng, type = "bspline", knots = knots,
                              norder = 4, lambda = lam)
    add_case(
      name = sprintf("smooth_sparse_mean_growth_bspline_lambda%g", lam),
      r_call = sprintf("smooth.sparse.mean(dl, NULL, rng=c(1,18), type='bspline', knots=seq(1,18,length=5), norder=4, lambda=%g)", lam),
      input = c(data_in, list(lambda = lam)),
      output = fd_output(mfd)
    )
  }

  mfd <- smooth.sparse.mean(dl, NULL, rng = rng, type = "bspline", knots = knots,
                            norder = 4, lambda = 0)
  cbasis <- create.bspline.basis(rng, 6)
  for (lam in c(0, 10)) {
    ce <- covPACE(dl, rng, age, mfd, cbasis, lam, 2)
    pc <- pcaPACE(ce, 3, fdPar(cbasis, 2, 0), FALSE)
    add_case(
      name = sprintf("pace_growth_bspline6_covlambda%g", lam),
      r_call = sprintf("pcaPACE(covPACE(dl, c(1,18), age, smooth.sparse.mean(...), create.bspline.basis(c(1,18), 6), %g, 2), 3, fdPar(cbasis, 2, 0), FALSE)", lam),
      input = c(data_in, list(
        mean_basis = basis_info(mfd$basis), mean_lambda = 0,
        cov_basis = basis_info(cbasis), cov_lambda = lam,
        harm_basis = basis_info(cbasis), harm_lambda = 0, nharm = 3
      )),
      output = c(list(mean_coefs = as.numeric(mfd$coefs),
                      cov_coefs = unname(ce$cov.estimate$coefs)), pca_output(pc, cbasis, cbasis))
    )
  }

  # Monomial bases: R's inprod() is exact for them (measured 0 against
  # eval.penalty), so every penalised quantity is comparable at 1e-8. A
  # monomial(4) covariance basis on [1, 18] makes covPACE's normal equations
  # "computationally singular" in R (t^3 up to 5832, squared by the tensor
  # product), so the covariance basis has three functions.
  mfd_m <- smooth.sparse.mean(dl, NULL, rng = rng, type = "mon", nbasis = 4, lambda = 0.5)
  mbasis <- create.monomial.basis(rng, 3)
  ce_m <- covPACE(dl, rng, age, mfd_m, mbasis, 2, 2)
  pc_m <- pcaPACE(ce_m, 2, fdPar(mbasis, 2, 5), FALSE)
  add_case(
    name = "pace_growth_monomial3_penalised",
    r_call = "pcaPACE(covPACE(dl, c(1,18), age, smooth.sparse.mean(dl, NULL, c(1,18), 'mon', nbasis=4, lambda=0.5), create.monomial.basis(c(1,18), 3), 2, 2), 2, fdPar(mbasis, 2, 5), FALSE)",
    input = c(data_in, list(
      mean_basis = basis_info(mfd_m$basis), mean_lambda = 0.5,
      cov_basis = basis_info(mbasis), cov_lambda = 2,
      harm_basis = basis_info(mbasis), harm_lambda = 5, nharm = 2
    )),
    output = c(list(mean_coefs = as.numeric(mfd_m$coefs),
                    cov_coefs = unname(ce_m$cov.estimate$coefs)), pca_output(pc_m, mbasis, mbasis))
  )
}

# ---- (b) Canadian weather temperature, 35 stations, 10 random days each ----

weather_cases <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, , "Temperature.C"]
  design <- sparse_design(length(day), 10, ncol(temp))
  dl <- as_pace_list(day, temp, design)
  data_in <- pace_input(dl)
  rng <- c(1, 365)

  mfd <- smooth.sparse.mean(dl, NULL, rng = rng, type = "fourier", nbasis = 5, lambda = 0)
  add_case(
    name = "smooth_sparse_mean_weather_fourier5",
    r_call = "smooth.sparse.mean(dl, NULL, rng=c(1,365), type='fourier', nbasis=5, lambda=0)",
    input = c(data_in, list(lambda = 0)),
    output = fd_output(mfd)
  )

  cbasis <- create.bspline.basis(rng, 7)
  ce <- covPACE(dl, rng, day, mfd, cbasis, 0, 2)
  pc <- pcaPACE(ce, 2, fdPar(cbasis, 2, 0), FALSE)
  sc <- scoresPACE(dl, day, ce, pc)
  add_case(
    name = "pace_weather_fourier5_bspline7",
    r_call = "scoresPACE(dl, 1:365, covPACE(dl, c(1,365), 1:365, smooth.sparse.mean(dl, NULL, c(1,365), 'fourier', nbasis=5, lambda=0), create.bspline.basis(c(1,365), 7), 0, 2), pcaPACE(ce, 2, fdPar(cbasis, 2, 0), FALSE))",
    input = c(data_in, list(
      mean_basis = basis_info(mfd$basis), mean_lambda = 0,
      cov_basis = basis_info(cbasis), cov_lambda = 0,
      harm_basis = basis_info(cbasis), harm_lambda = 0, nharm = 2
    )),
    output = c(list(mean_coefs = as.numeric(mfd$coefs),
                    cov_coefs = unname(ce$cov.estimate$coefs),
                    scores = unname(sc)), pca_output(pc, cbasis, cbasis))
  )
}

growth_cases()
weather_cases()

finalize("pace", rtol = 1e-8)
