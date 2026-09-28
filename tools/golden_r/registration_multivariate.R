# tools/golden_r/registration_multivariate.R
#
# Golden cases for multivariate registration: register.fd on curves with two
# variables (gait hip + knee angles, handwriting x + y), register.newfd, and
# landmarkreg run variable by variable. Driven by tools/make_golden.py, which
# sources common.R first (defines add_case, basis_info, finalize, GOLDEN_OUT,
# GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# vec2Lfd, fdPar, fd, smooth.basis, mean.fd, eval.fd, register.fd,
# register.newfd, landmarkreg). No fda source is read or copied. The
# behaviour notes below come from calling these functions and comparing
# their outputs (black box), and from the installed help pages.
#
# fda behaviour discovered while writing this script (see PROGRESS.md):
#   - register.fd() on multivariate curves (3-D coefs) returns warps that are
#     fitted to the FIRST variable only. Its Wfd equals, to the last digit,
#     the Wfd of a univariate register.fd() on variable 1 alone, and swapping
#     the variables gives the Wfd of variable 2 alone. The other variables
#     are only warped along (regfd keeps all variables). Every case below
#     records the univariate first-variable Wfd too (`Wfd_first_variable`),
#     so the parity test can show the identity.
#   - register.fd() rejects a target with one curve per variable given as a
#     2-D coefficient matrix ("Second dimension of coefficient matrix for
#     Y0FD is neither 1 nor equal to the number of functions to be
#     registered"); a multivariate target must have 3-D coefs.
#   - landmarkreg() rejects multivariate curves ("incorrect number of
#     dimensions"). Its warps depend only on the landmarks, so the landmark
#     case runs it on each variable with the same landmarks; both runs give
#     the same Wfd, which is what a multivariate landmark registration must
#     return.
#   - register.newfd() needs as many curves as warps ("There must be as many
#     warping function replicates as y replicates").
#   - AmpPhaseDecomp() rejects multivariate curves ("non-conformable
#     arrays"), so no multivariate decomposition is recorded.

gait_fd <- function() {
  gaittime <- as.matrix(0:19) + 0.5
  gaitrange <- c(0, 20)
  gaitbasis <- create.fourier.basis(gaitrange, nbasis = 21)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 20)^2, 0), rangeval = gaitrange)
  smooth.basis(gaittime, gait, fdPar(gaitbasis, harmaccelLfd, 1e-2))$fd
}

handwriting_fd <- function() {
  basis <- create.bspline.basis(range(handwritTime), nbasis = 21, norder = 6)
  # handwrit is in metres (values ~0.03), so the registration criterion is
  # ~1e-7 and register.fd stops at its starting point; millimetres give it
  # a criterion of order one.
  smooth.basis(handwritTime, 1000 * handwrit, fdPar(basis, 2, 1e-4))$fd
}

first_variable <- function(x) {
  coefs <- x$coefs
  fd(matrix(coefs[, , 1], nrow = dim(coefs)[1], ncol = dim(coefs)[2]), x$basis)
}

continuous_case <- function(name, yfd, wbasis, lambda, crit, periodic, conv, iterlim) {
  y0fd <- mean.fd(yfd)
  ncurve <- dim(yfd$coefs)[2]
  WfdParobj <- fdPar(fd(matrix(0, wbasis$nbasis, ncurve), wbasis), 2, lambda)
  reg <- register.fd(y0fd, yfd, WfdParobj, conv = conv, iterlim = iterlim, dbglev = 0,
                     periodic = periodic, crit = crit)
  first <- register.fd(first_variable(y0fd), first_variable(yfd), WfdParobj, conv = conv,
                       iterlim = iterlim, dbglev = 0, periodic = periodic, crit = crit)
  newfd <- register.newfd(yfd, reg$Wfd, type = "monotone")
  add_case(
    name = name,
    r_call = sprintf(
      "register.fd(mean.fd(yfd), yfd, fdPar(Wfd0, 2, %g), conv=%g, iterlim=%d, dbglev=0, periodic=%s, crit=%d); register.newfd(yfd, reg$Wfd, type='monotone')",
      lambda, conv, iterlim, periodic, crit
    ),
    input = list(
      basis = basis_info(yfd$basis), wbasis = basis_info(wbasis),
      coefs = unname(yfd$coefs), y0fd_coefs = unname(y0fd$coefs),
      lambda = lambda, penalty = 2, conv = conv, iterlim = iterlim, crit = crit,
      periodic = periodic
    ),
    output = list(
      regfd_coefs = unname(reg$regfd$coefs),
      warpfd_coefs = unname(reg$warpfd$coefs),
      Wfd_coefs = unname(reg$Wfd$coefs),
      shift = unname(as.numeric(reg$shift)),
      Wfd_first_variable = unname(first$Wfd$coefs),
      shift_first_variable = unname(as.numeric(first$shift)),
      newfd_coefs = unname(newfd$coefs)
    ),
    rtol = 1e-5
  )
}

# ---- register.fd: gait hip + knee, 5 boys, eigen criterion ----------------

gait_case <- function() {
  gaitfd <- gait_fd()
  wbasis <- create.bspline.basis(c(0, 20), norder = 4, nbasis = 5)
  continuous_case("register_fd_gait_hip_knee_crit2", gaitfd[1:5], wbasis,
                  lambda = 1e-2, crit = 2, periodic = FALSE, conv = 1e-4, iterlim = 20)
}

# ---- register.fd: gait hip + knee, periodic shift, least squares ----------

gait_periodic_case <- function() {
  gaitfd <- gait_fd()
  wbasis <- create.bspline.basis(c(0, 20), norder = 4, nbasis = 5)
  continuous_case("register_fd_gait_hip_knee_periodic_crit1", gaitfd[1:4], wbasis,
                  lambda = 1e-2, crit = 1, periodic = TRUE, conv = 1e-4, iterlim = 20)
}

# ---- register.fd: handwriting x + y, 5 samples, eigen criterion -----------

handwriting_case <- function() {
  hwfd <- handwriting_fd()
  wbasis <- create.bspline.basis(range(handwritTime), norder = 4, nbasis = 5)
  continuous_case("register_fd_handwriting_xy_crit2", hwfd[1:5], wbasis,
                  lambda = 1, crit = 2, periodic = FALSE, conv = 1e-4, iterlim = 20)
}

# ---- landmarkreg: gait hip + knee, knee-peak landmark, per variable -------

landmark_case <- function() {
  gaitfd <- gait_fd()[1:6]
  # landmark = time of the knee-angle maximum over the swing phase
  tfine <- seq(10, 20, length.out = 1001)
  knee <- eval.fd(tfine, gaitfd)[, , 2]
  ximarks <- apply(knee, 2, function(v) tfine[which.max(v)])
  x0marks <- mean(ximarks)
  wbasis <- create.bspline.basis(c(0, 20), norder = 4, breaks = c(0, x0marks, 20))
  Wfd0 <- fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar <- fdPar(Wfd0, 2, 1e-4)
  per_var <- lapply(1:2, function(v) {
    uni <- fd(matrix(gaitfd$coefs[, , v], nrow = 21), gaitfd$basis)
    landmarkreg(uni, matrix(ximarks, ncol = 1), x0marks, WfdPar = WfdPar, WfdPar0 = WfdPar)
  })
  add_case(
    name = "landmarkreg_gait_hip_knee_per_variable",
    r_call = "for v in hip, knee: landmarkreg(gaitfd[1:6] variable v, matrix(ximarks,ncol=1), x0marks, WfdPar=fdPar(Wfd0,2,1e-4), WfdPar0=same)",
    input = list(
      basis = basis_info(gaitfd$basis), wbasis = basis_info(wbasis),
      coefs = unname(gaitfd$coefs), ximarks = unname(as.numeric(ximarks)), x0marks = x0marks,
      lambda = 1e-4, penalty = 2
    ),
    output = list(
      Wfd_coefs_hip = unname(per_var[[1]]$Wfd$coefs),
      Wfd_coefs_knee = unname(per_var[[2]]$Wfd$coefs)
    ),
    rtol = 1e-5
  )
}

gait_case()
gait_periodic_case()
handwriting_case()
landmark_case()

finalize("registration_multivariate", rtol = 1e-5)
