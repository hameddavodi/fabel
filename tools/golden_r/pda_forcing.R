# tools/golden_r/pda_forcing.R
#
# Golden cases for principal differential analysis with forcing functions
# (pda.fd's awtlist/ufdlist) and for its stability analysis (eigen.pda).
# Driven by tools/make_golden.py, which sources common.R first (defines
# add_case, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis, fd,
# fdPar, smooth.basis, pda.fd, eigen.pda, eval.fd). The forced examples
# adapt the worked examples of pda.fd's own help page (examples 5 and 6),
# the refinery case follows the reflux -> tray 47 example of Ramsay, Hooker
# & Graves (2009, ch. 11) and the handwriting system is the example of
# eigen.pda's help page. No fda source is read or copied.
#
# Every PDA is recorded together with the exact curves it was fitted to
# (their basis and coefficients), so the Python side rebuilds R's inputs
# instead of re-smoothing them.
#
# fda API quirks discovered while writing this script (black-box behaviour):
#   - pda.fd needs every forcing function to have as many replicates as the
#     curves ("non-conformable arrays" otherwise); a single forcing curve is
#     not broadcast over several curves.
#   - pda.fd does not return ufdlist, and eigen.pda reads pdaList$ufdlist:
#     called on the pda.fd result of a forced equation it stops with
#     "rep(1, nbasis): invalid 'times' argument". Attaching the forcing
#     functions (flat list for one equation, nested per equation for a
#     system) to the result first makes it run; that is what is done below.
#   - eigen.pda's argvals are seq(lower, upper, length.out = npts); its
#     eigenvalues are those of the companion matrix of the first-order form
#     of the fitted system, ordered by decreasing modulus.
#   - pda.fd raises its 501-point grid to 5 x the largest number of basis
#     functions of the curves (120 for the handwriting basis -> 600).

pda_output <- function(result, forced = TRUE) {
  out <- list(
    bwt_coefs = lapply(result$bwtlist, function(b) unname(b$fd$coefs)),
    resfd_coefs = lapply(result$resfdlist, function(r) unname(r$coefs))
  )
  if (forced) {
    out$awt_coefs <- lapply(result$awtlist, function(a) unname(a$fd$coefs))
  }
  out
}

system_output <- function(result, nvar, order, nforce) {
  bwt <- list()
  for (i in seq_len(nvar)) {
    for (k in seq_len(nvar)) {
      for (j in seq_len(order)) {
        bwt[[length(bwt) + 1]] <- unname(result$bwtlist[[i]][[k]][[j]]$fd$coefs)
      }
    }
  }
  awt <- list()
  if (nforce > 0) {
    for (i in seq_len(nvar)) {
      for (l in seq_len(nforce)) {
        awt[[length(awt) + 1]] <- unname(result$awtlist[[i]][[l]]$fd$coefs)
      }
    }
  }
  list(
    bwt_coefs_flat = bwt,
    awt_coefs_flat = awt,
    resfd_coefs = lapply(result$resfdlist, function(r) unname(r$coefs))
  )
}

eigen_output <- function(result, npts) {
  e <- eigen.pda(result, plotresult = FALSE, npts = npts)
  ev <- as.matrix(e$eigvals)
  list(
    argvals = as.numeric(e$argvals),
    eig_re = unname(Re(ev)),
    eig_im = unname(Im(ev)),
    limvals = unname(as.matrix(e$limvals))
  )
}

fd_info <- function(f) {
  list(basis = basis_info(f$basis), coefs = unname(f$coefs))
}

# ---- (1) help example 5: Dx = -4x + 2u, u = 1 ---------------------------

forced_exp_case <- function() {
  cbasis <- create.constant.basis(c(0, 1))
  xbasis <- create.bspline.basis(c(0, 1), 24, 5)
  cfdPar <- fdPar(fd(0, cbasis))
  tvec <- seq(0, 1, len = 101)
  beta <- 4
  alpha <- 2
  xvec <- exp(-beta * tvec) * (1 + alpha * (exp(beta * tvec) - 1) / beta)
  xfd <- smooth.basis(tvec, xvec, xbasis)$fd
  ufd <- fd(1, cbasis)
  result <- pda.fd(list(xfd), list(cfdPar), list(cfdPar), list(ufd))
  add_case(
    name = "pda_forced_order1_constant_exp",
    r_call = "pda.fd(list(xfd), list(cfdPar), list(cfdPar), list(fd(1, cbasis)))  # pda.fd help example 5",
    input = list(x = fd_info(xfd), u = list(fd_info(ufd)), wbasis = basis_info(cbasis)),
    output = pda_output(result)
  )
}

# ---- (2) three curves, three forcing curves, penalised B-spline weights --

forced_multicurve_case <- function() {
  xbasis <- create.bspline.basis(c(0, 1), 13, 5)
  wbasis <- create.bspline.basis(c(0, 1), 5, 4)
  tvec <- seq(0, 1, len = 101)
  X <- cbind(
    exp(-4 * tvec) + 0.5 * (1 - exp(-4 * tvec)),
    exp(-3 * tvec) + 0.1 * sin(5 * tvec),
    cos(2 * tvec)
  )
  xfd <- smooth.basis(tvec, X, xbasis)$fd
  ufd <- fd(cbind(c(1, 1.2, 0.8, 1.1, 1.3), c(0, 1, 0, 1, 0), c(2, 1, 2, 1, 2)), wbasis)
  bwtPar <- fdPar(fd(matrix(0, 5, 1), wbasis), 2, 1e-3)
  awtPar <- fdPar(fd(matrix(0, 5, 1), wbasis), 2, 1e-2)
  result <- pda.fd(list(xfd), list(bwtPar), list(awtPar), list(ufd))
  grid <- seq(0, 1, length.out = 21)
  add_case(
    name = "pda_forced_multicurve_bspline_weights",
    r_call = "pda.fd(list(xfd), list(fdPar(wfd, 2, 1e-3)), list(fdPar(wfd, 2, 1e-2)), list(ufd))  # 3 curves, 3 forcing curves",
    input = list(
      x = fd_info(xfd), u = list(fd_info(ufd)), wbasis = basis_info(wbasis),
      bwt_lambda = 1e-3, awt_lambda = 1e-2, penalty = 2
    ),
    output = c(
      pda_output(result),
      list(
        grid = grid,
        awt_vals = unname(as.numeric(eval.fd(grid, result$awtlist[[1]]$fd))),
        bwt_vals = unname(as.numeric(eval.fd(grid, result$bwtlist[[1]]$fd)))
      )
    )
  )
}

# ---- refinery: tray 47 level forced by reflux ----------------------------

refinery_fds <- function() {
  rng <- c(0, 193)
  breaks <- c(seq(0, 60, by = 10), 67, 68, seq(70, 190, by = 10), 193)
  basis <- create.bspline.basis(rng, norder = 4, breaks = breaks)
  par <- fdPar(basis, 2, 1e-2)
  yfd <- smooth.basis(refinery$Time, refinery$Tray47, par)$fd
  ufd <- smooth.basis(refinery$Time, refinery$Reflux, par)$fd
  list(yfd = yfd, ufd = ufd, rng = rng)
}

refinery_constant_case <- function() {
  f <- refinery_fds()
  cbasis <- create.constant.basis(f$rng)
  cfdPar <- fdPar(fd(0, cbasis))
  result <- pda.fd(list(f$yfd), list(cfdPar), list(cfdPar), list(f$ufd))
  result$ufdlist <- list(f$ufd)
  add_case(
    name = "pda_forced_refinery_constant",
    r_call = "pda.fd(list(trayfd), list(cfdPar), list(cfdPar), list(refluxfd)); eigen.pda(<result + ufdlist>, npts = 51)",
    input = list(x = fd_info(f$yfd), u = list(fd_info(f$ufd)), wbasis = basis_info(cbasis)),
    output = c(pda_output(result), eigen_output(result, 51))
  )
}

refinery_two_forcings_case <- function() {
  f <- refinery_fds()
  cbasis <- create.constant.basis(f$rng)
  wbasis <- create.bspline.basis(f$rng, nbasis = 6, norder = 4)
  bwtPar <- fdPar(fd(matrix(0, 6, 1), wbasis), 2, 1e4)
  awtPar <- fdPar(fd(matrix(0, 6, 1), wbasis), 2, 1e5)
  cfdPar <- fdPar(fd(0, cbasis))
  onefd <- fd(1, cbasis)
  result <- pda.fd(list(f$yfd), list(bwtPar), list(awtPar, cfdPar), list(f$ufd, onefd))
  result$ufdlist <- list(f$ufd, onefd)
  add_case(
    name = "pda_forced_refinery_two_forcings",
    r_call = "pda.fd(list(trayfd), list(fdPar(wfd,2,1e4)), list(fdPar(wfd,2,1e5), cfdPar), list(refluxfd, fd(1,cbasis))); eigen.pda(..., npts = 41)",
    input = list(
      x = fd_info(f$yfd), u = list(fd_info(f$ufd), fd_info(onefd)),
      wbasis = basis_info(wbasis), cbasis = basis_info(cbasis),
      bwt_lambda = 1e4, awt_lambda = 1e5, penalty = 2
    ),
    output = c(pda_output(result), eigen_output(result, 41))
  )
}

# ---- (5) second order, constant forcing: D2x = -b0 x - b1 Dx + a --------

forced_order2_case <- function() {
  cbasis <- create.constant.basis(c(0, 1))
  xbasis <- create.bspline.basis(c(0, 1), 24, 5)
  cfdPar <- fdPar(fd(0, cbasis))
  tvec <- seq(0, 1, len = 101)
  X <- cbind(
    0.5 + sin(2 * pi * tvec) + 0.05 * tvec^2,
    0.5 + 0.7 * cos(2 * pi * tvec) - 0.03 * tvec,
    0.5 + 0.4 * sin(2 * pi * tvec + 1)
  )
  xfd <- smooth.basis(tvec, X, xbasis)$fd
  ufd <- fd(matrix(1, 1, 3), cbasis)
  result <- pda.fd(list(xfd), list(cfdPar, cfdPar), list(cfdPar), list(ufd))
  result$ufdlist <- list(fd(1, cbasis))
  add_case(
    name = "pda_forced_order2_constant",
    r_call = "pda.fd(list(xfd), list(cfdPar, cfdPar), list(cfdPar), list(fd(matrix(1,1,3), cbasis))); eigen.pda(..., npts = 11)",
    input = list(x = fd_info(xfd), u = list(fd_info(ufd)), wbasis = basis_info(cbasis)),
    output = c(pda_output(result), eigen_output(result, 11))
  )
}

# ---- (6) help example 6: two forced first-order equations ---------------
# The help page gives beta_22 a monomial (linear) basis; here every weight is
# constant, because Fabel's PDA takes one weight basis per derivative order,
# shared by all variables and equations.

forced_system_case <- function() {
  cbasis <- create.constant.basis(c(0, 1))
  xbasis <- create.bspline.basis(c(0, 1), 24, 5)
  cfdPar <- fdPar(fd(0, cbasis))
  tvec <- seq(0, 1, len = 101)
  beta <- 4
  xvec1 <- exp(-beta * tvec) + 2 * (1 - exp(-beta * tvec)) / beta
  vvec <- exp(beta * tvec^2 / 2)
  intv <- 0.01 * (cumsum(vvec) - 0.5 * vvec)
  xvec2 <- exp(-beta * tvec^2 / 2) * (1 - intv)
  xfd1 <- smooth.basis(tvec, xvec1, xbasis)$fd
  xfd2 <- smooth.basis(tvec, xvec2, xbasis)$fd
  bwtlist <- list(list(list(cfdPar), list(cfdPar)), list(list(cfdPar), list(cfdPar)))
  awtlist <- list(list(cfdPar), list(cfdPar))
  ufdlist <- list(list(fd(1, cbasis)), list(fd(1, cbasis)))
  result <- pda.fd(list(xfd1, xfd2), bwtlist, awtlist, ufdlist)
  result$ufdlist <- ufdlist
  add_case(
    name = "pda_forced_system_order1",
    r_call = "pda.fd(list(xfd1, xfd2), bwtlist, awtlist, ufdlist)  # pda.fd help example 6, all weights constant; eigen.pda(..., npts = 11)",
    input = list(
      x = list(fd_info(xfd1), fd_info(xfd2)),
      u = list(fd_info(fd(1, cbasis))),
      cbasis = basis_info(cbasis)
    ),
    output = c(system_output(result, 2, 1, 1), eigen_output(result, 11))
  )
}

# ---- (7) lip, unforced, B-spline weights: stability only ----------------

lip_eigen_case <- function() {
  liprange <- range(liptime)
  lipbasis <- create.bspline.basis(liprange, nbasis = 31, norder = 6)
  lipfd <- smooth.basis(liptime, lip, fdPar(lipbasis, 4, 1e-8))$fd
  wbasis <- create.bspline.basis(liprange, nbasis = 11, norder = 4)
  wfdPar <- fdPar(fd(matrix(0, 11, 1), wbasis), 2, 1e-6)
  result <- pda.fd(list(lipfd), list(wfdPar, wfdPar))
  add_case(
    name = "eigen_pda_lip_bspline11",
    r_call = "eigen.pda(pda.fd(list(lipfd), list(wfdPar, wfdPar)), npts = 51)",
    input = list(x = fd_info(lipfd), wbasis = basis_info(wbasis), lambda = 1e-6, penalty = 2),
    output = c(pda_output(result, forced = FALSE), eigen_output(result, 51))
  )
}

# ---- (8) handwriting: coupled second-order system (eigen.pda example) ---

handwriting_eigen_case <- function() {
  ni <- 281
  indx <- seq(1, 1401, length = ni)
  fdaarray <- handwrit[indx, , ]
  fdatime <- seq(0, 2.3, len = ni)
  fdarange <- c(0, 2.3)
  breaks <- seq(0, 2.3, length.out = 116)
  fdabasis <- create.bspline.basis(fdarange, norder = 6, breaks = breaks)
  nbasis <- fdabasis$nbasis
  fdaPar <- fdPar(fd(matrix(0, nbasis, 1), fdabasis), int2Lfd(4), 1e-8)
  Xfd <- smooth.basis(fdatime, fdaarray[, , 1], fdaPar)$fd
  Yfd <- smooth.basis(fdatime, fdaarray[, , 2], fdaPar)$fd
  wbasis <- create.bspline.basis(fdarange, norder = 6, nbasis = 31)
  pdaPar <- fdPar(fd(matrix(0, 31, 1), wbasis), 1, 1e-8)
  pdaParlist <- list(pdaPar, pdaPar)
  bwtlist <- list(list(pdaParlist, pdaParlist), list(pdaParlist, pdaParlist))
  result <- pda.fd(list(Xfd, Yfd), bwtlist)
  add_case(
    name = "eigen_pda_handwriting_system",
    r_call = "eigen.pda(pda.fd(list(Xfd, Yfd), bwtlist), npts = 51)  # eigen.pda help example",
    input = list(
      x = list(fd_info(Xfd), fd_info(Yfd)),
      wbasis = basis_info(wbasis), lambda = 1e-8, penalty = 1, nfine = 5 * nbasis
    ),
    output = c(system_output(result, 2, 2, 0), eigen_output(result, 51))
  )
}

forced_exp_case()
forced_multicurve_case()
refinery_constant_case()
refinery_two_forcings_case()
forced_order2_case()
forced_system_case()
lip_eigen_case()
handwriting_eigen_case()

finalize("pda_forcing", rtol = 1e-8)
