# tools/golden_r/dynamics.R
#
# Golden cases for the `dynamics` module: pda.fd (Principal Differential
# Analysis). Driven by tools/make_golden.py, which sources common.R first
# (defines add_case, eval_points, basis_info, finalize, GOLDEN_OUT,
# GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# fdPar, smooth.basis, pda.fd, fd, eval.fd). The order-1 analytic case
# below adapts the FIRST worked example from pda.fd's own installed help
# page (tools::Rd_db("fda")[["pda.fd.Rd"]], examples section: "a single
# first order constant coefficient unforced equation Dx = -4*x for
# x(t) = exp(-4t)") -- this is fda's own official documented usage
# pattern, not internal source. No fda source is read or copied.
#
# fda API quirks / notes discovered while writing this script:
#   - bwtlist nesting convention (confirmed via pda.fd's own Rd examples,
#     not memorized): for a SINGLE equation of order m, bwtlist is a flat
#     list of m fdPar objects (e.g. order-2: list(cfdPar, cfdPar) -- one
#     for the variable itself (deriv 0), one for its first derivative).
#     A multi-equation system nests three levels deep: outer list per
#     equation, middle list per variable, inner list per derivative
#     order -- not needed here since every case below is single-equation.
#   - pda.fd()'s returned bwtlist / resfdlist preserve the coefficient
#     dimensions of the WfdPar/xfdlist basis used, not any fixed shape:
#     a constant-basis bwt has coefs dim [1,1]; an 11-basis B-spline bwt
#     has coefs dim [11,1]; resfdlist coefs match the xfdlist basis's
#     nbasis x N (curves).
#   - eval.fd() of the pda-derived weight function on a grid gives
#     directly comparable values to a Python re-fit for cross-checking,
#     and for the constant-coefficient analytic case we additionally
#     record the true analytic solution x(t)=exp(-beta*t) on the same
#     grid so Python's solve_ivp can be checked against both fda's
#     estimate and ground truth.

# ---- order-1 analytic case: Dx = -4x, x(t) = exp(-4t) (official Rd ex) --

order1_analytic_case <- function() {
  cbasis <- create.constant.basis(c(0, 1))
  xbasis <- create.bspline.basis(c(0, 1), 24, 5)
  cfd0 <- fd(0, cbasis)
  cfdPar <- fdPar(cfd0)
  tvec <- seq(0, 1, len = 101)
  beta_true <- 4
  xvec <- exp(-beta_true * tvec)
  xfd <- smooth.basis(tvec, xvec, xbasis)$fd
  xfdlist <- list(xfd)
  bwtlist <- list(cfdPar)

  result <- pda.fd(xfdlist, bwtlist)

  grid <- seq(0, 1, length.out = 51)
  bwt_vals <- eval.fd(grid, result$bwtlist[[1]]$fd)
  analytic <- exp(-beta_true * grid)

  add_case(
    name = "pda_fd_order1_analytic_exp_decay",
    r_call = "pda.fd(list(xfd), list(cfdPar))  # xfd: bspline5(24basis) fit of x(t)=exp(-4t), cbasis weight",
    input = list(
      beta_true = beta_true, tvec = tvec, xvec = xvec,
      xbasis = basis_info(xbasis), cbasis = basis_info(cbasis)
    ),
    output = list(
      bwt_coefs = unname(result$bwtlist[[1]]$fd$coefs),
      resfd_coefs = unname(result$resfdlist[[1]]$coefs),
      bwt_grid = grid, bwt_vals = unname(as.numeric(bwt_vals)),
      analytic_solution = unname(analytic)
    )
  )
}

# ---- order-2 case on real lip data: constant-basis bwtlist -------------

order2_lip_constant_case <- function() {
  liprange <- range(liptime)
  lipbasis <- create.bspline.basis(liprange, nbasis = 31, norder = 6)
  lipfd <- smooth.basis(liptime, lip, fdPar(lipbasis, 4, 1e-8))$fd

  cbasis <- create.constant.basis(liprange)
  cfd0 <- fd(0, cbasis)
  cfdPar <- fdPar(cfd0)
  xfdlist <- list(lipfd)
  bwtlist <- list(cfdPar, cfdPar)

  result <- pda.fd(xfdlist, bwtlist)
  grid <- seq(liprange[1], liprange[2], length.out = 51)
  add_case(
    name = "pda_fd_order2_lip_constant",
    r_call = "pda.fd(list(lipfd), list(cfdPar, cfdPar))  # lipfd: bspline6(31basis) fit of lip data, cbasis weights",
    input = list(
      liptime = as.numeric(liptime), lip = unname(lip),
      xbasis = basis_info(lipbasis), cbasis = basis_info(cbasis)
    ),
    output = list(
      bwt0_coefs = unname(result$bwtlist[[1]]$fd$coefs),
      bwt1_coefs = unname(result$bwtlist[[2]]$fd$coefs),
      resfd_coefs = unname(result$resfdlist[[1]]$coefs),
      bwt0_vals = unname(as.numeric(eval.fd(grid, result$bwtlist[[1]]$fd))),
      bwt1_vals = unname(as.numeric(eval.fd(grid, result$bwtlist[[2]]$fd))),
      bwt_grid = grid
    )
  )
}

# ---- order-2 case on real lip data: B-spline(11) bwtlist ----------------

order2_lip_bspline_case <- function() {
  liprange <- range(liptime)
  lipbasis <- create.bspline.basis(liprange, nbasis = 31, norder = 6)
  lipfd <- smooth.basis(liptime, lip, fdPar(lipbasis, 4, 1e-8))$fd

  wbasis11 <- create.bspline.basis(liprange, nbasis = 11, norder = 4)
  wfd0 <- fd(matrix(0, 11, 1), wbasis11)
  wfdPar <- fdPar(wfd0, 2, 1e-6)
  xfdlist <- list(lipfd)
  bwtlist <- list(wfdPar, wfdPar)

  result <- pda.fd(xfdlist, bwtlist)
  grid <- seq(liprange[1], liprange[2], length.out = 51)
  add_case(
    name = "pda_fd_order2_lip_bspline11",
    r_call = "pda.fd(list(lipfd), list(wfdPar, wfdPar))  # lipfd: bspline6(31basis) fit; weight fns: bspline4(11basis), lambda=1e-6",
    input = list(
      liptime = as.numeric(liptime), lip = unname(lip),
      xbasis = basis_info(lipbasis), wbasis = basis_info(wbasis11)
    ),
    output = list(
      bwt0_coefs = unname(result$bwtlist[[1]]$fd$coefs),
      bwt1_coefs = unname(result$bwtlist[[2]]$fd$coefs),
      resfd_coefs = unname(result$resfdlist[[1]]$coefs),
      bwt0_vals = unname(as.numeric(eval.fd(grid, result$bwtlist[[1]]$fd))),
      bwt1_vals = unname(as.numeric(eval.fd(grid, result$bwtlist[[2]]$fd))),
      bwt_grid = grid
    )
  )
}

order1_analytic_case()
order2_lip_constant_case()
order2_lip_bspline_case()

finalize("dynamics", rtol = 1e-8)
