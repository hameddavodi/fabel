# tools/golden_r/registration.R
#
# Golden cases for the `registration` module: landmarkreg, register.fd,
# AmpPhaseDecomp. Driven by tools/make_golden.py, which sources common.R
# first (defines add_case, eval_points, basis_info, finalize, GOLDEN_OUT,
# GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# fdPar, smooth.basis, landmarkreg, register.fd, AmpPhaseDecomp, fd,
# mean.fd, eval.fd). No fda source is read or copied. Field names and
# argument shapes below were confirmed by calling args() and by fetching
# the installed help pages via tools::Rd_db("fda") + tools::Rd2txt(),
# then checking actual runtime output -- not from memorized API.
#
# fda API quirks discovered while writing this script (see PROGRESS.md):
#   - landmarkreg()'s installed help page (Rd) documents the return list
#     as `fdreg, warpfd, warpinvfd, Wfd`, but the ACTUAL 6.3.0 runtime
#     return list is `regfd, warpfd, warpinvfd, Wfd` (regfd, not fdreg;
#     confirmed by calling names() on the live result -- $fdreg is NULL,
#     $regfd holds the registered-curve fd). Recorded fields use the
#     runtime name, regfd.
#   - landmarkreg() in 6.3.0 has NO monwrd argument at all (args() shows
#     unregfd, ximarks, x0marks, x0lim, WfdPar, WfdPar0, ylambda) -- the
#     monwrd parameter mentioned in older references does not exist in
#     this installed version, so no monwrd=TRUE/FALSE variants are
#     recorded.
#   - landmarkreg()'s WfdPar/WfdPar0 must each wrap a SINGLE-curve fd
#     (coefs ncol=1), used internally as a per-curve template -- passing
#     a multi-curve (ncol=N) Wfd0 crashes deep inside smooth.morph with
#     "non-conformable arguments". WfdPar0 must be supplied explicitly
#     (its default NULL crashes create.bspline.basis for small landmark
#     counts).
#   - Peak-velocity landmarks for the growth curves must be searched
#     over an INTERIOR grid (here [8, 16] years), not the full age
#     range: derivative estimates near the domain boundary of a
#     boundary-heavy B-spline fit spike sharply (an edge artifact, not
#     a real growth feature), which spuriously wins argmax over the
#     true pubertal growth spurt if the full range is searched.
#   - register.fd()'s WfdParobj basis must be a B-spline basis even for
#     a periodic=TRUE registration of Fourier-basis data -- passing a
#     Fourier basis for Wfd raises "Basis for Wfd is not a B-spline
#     basis." The curves being registered (yfd/y0fd) may still be
#     Fourier; only the warping-function basis is constrained.
#   - AmpPhaseDecomp(xfd, yfd, hfd)'s first argument must be the ORIGINAL
#     UNREGISTERED curves (same N as yfd/hfd), not the 1-replicate
#     registration target -- register.fd()'s own return object
#     conveniently echoes the original yfd input back as $yfd, which is
#     the right thing to pass as AmpPhaseDecomp's first argument.

# ---- landmarkreg: growth hgtf velocity curves, 10 girls ------------------

landmarkreg_case <- function() {
  age <- growth$age
  basis <- create.bspline.basis(range(age), breaks = age, norder = 6)
  hgtf_fd <- smooth.basis(age, growth$hgtf, fdPar(basis, 4, 0.01))$fd
  girls_idx <- 1:10
  hgtf_fd10 <- hgtf_fd[girls_idx]

  # peak-velocity landmark = argmax of height velocity over an interior
  # grid (avoids the boundary-derivative artifact noted above)
  agefine <- seq(8, 16, length.out = 501)
  vel <- eval.fd(agefine, hgtf_fd10, 1)
  ximarks <- apply(vel, 2, function(v) agefine[which.max(v)])
  x0marks <- mean(ximarks)

  wbasis <- create.bspline.basis(range(age), norder = 4, breaks = c(1, 4, 18))
  Wfd0 <- fd(matrix(0, wbasis$nbasis, 1), wbasis)
  WfdPar <- fdPar(Wfd0, 2, 1)
  WfdPar0 <- fdPar(Wfd0, 2, 1)

  res <- landmarkreg(hgtf_fd10, matrix(ximarks, ncol = 1), x0marks, WfdPar = WfdPar, WfdPar0 = WfdPar0)
  add_case(
    name = "landmarkreg_growth_hgtf_pubertal_spurt",
    r_call = "landmarkreg(hgtf_fd10, matrix(ximarks,ncol=1), x0marks, WfdPar=WfdPar, WfdPar0=WfdPar0)",
    input = list(
      ximarks = unname(as.numeric(ximarks)), x0marks = x0marks,
      basis = basis_info(basis), wbasis = basis_info(wbasis),
      coefs = unname(hgtf_fd10$coefs)
    ),
    output = list(
      regfd_coefs = unname(res$regfd$coefs),
      warpfd_coefs = unname(res$warpfd$coefs),
      warpinvfd_coefs = unname(res$warpinvfd$coefs),
      Wfd_coefs = unname(res$Wfd$coefs)
    ),
    rtol = 1e-5
  )
}

# ---- register.fd: growth hgtf curves registered to their mean ------------

register_fd_growth_case <- function() {
  age <- growth$age
  basis <- create.bspline.basis(range(age), breaks = age, norder = 6)
  hgtf_fd <- smooth.basis(age, growth$hgtf, fdPar(basis, 4, 0.01))$fd
  girls_idx <- 1:10
  hgtf_fd10 <- hgtf_fd[girls_idx]
  y0fd <- mean.fd(hgtf_fd10)

  wbasis <- create.bspline.basis(range(age), norder = 4, nbasis = 5)
  Wfd0 <- fd(matrix(0, wbasis$nbasis, 10), wbasis)
  WfdParobj <- fdPar(Wfd0, 2, 1)

  reg <- register.fd(y0fd, hgtf_fd10, WfdParobj, conv = 1e-4, iterlim = 20, dbglev = 0, periodic = FALSE, crit = 2)
  add_case(
    name = "register_fd_growth_hgtf_to_mean",
    r_call = "register.fd(y0fd, hgtf_fd10, WfdParobj, conv=1e-4, iterlim=20, dbglev=0, periodic=FALSE, crit=2)",
    input = list(
      basis = basis_info(basis), wbasis = basis_info(wbasis),
      coefs = unname(hgtf_fd10$coefs), y0fd_coefs = unname(y0fd$coefs),
      conv = 1e-4, iterlim = 20, crit = 2, periodic = FALSE
    ),
    output = list(
      regfd_coefs = unname(reg$regfd$coefs),
      warpfd_coefs = unname(reg$warpfd$coefs),
      Wfd_coefs = unname(reg$Wfd$coefs),
      shift = unname(as.numeric(reg$shift))
    ),
    rtol = 1e-5
  )

  amp <- AmpPhaseDecomp(reg$yfd, reg$regfd, reg$warpfd)
  add_case(
    name = "ampphasedecomp_growth_hgtf_to_mean",
    r_call = "AmpPhaseDecomp(reg$yfd, reg$regfd, reg$warpfd)  # reg = register_fd_growth_hgtf_to_mean",
    input = list(source_case = "register_fd_growth_hgtf_to_mean"),
    output = list(MS.amp = amp$MS.amp, MS.pha = amp$MS.pha, RSQR = amp$RSQR, C = amp$C),
    rtol = 1e-5
  )
}

# ---- register.fd: periodic case on 5 weather Fourier curves, crit=1 -----

register_fd_periodic_case <- function() {
  day <- 1:365
  temp <- CanadianWeather$dailyAv[, 1:5, "Temperature.C"]
  basis <- create.fourier.basis(c(0, 365), nbasis = 65)
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / 365)^2, 0), c(0, 365))
  temp_fd <- smooth.basis(day, temp, fdPar(basis, harmaccelLfd, 1e2))$fd
  y0fd <- mean.fd(temp_fd)

  # Wfd basis must be B-spline even for periodic registration of Fourier
  # data (see header note) -- the curves themselves stay Fourier.
  wbasis <- create.bspline.basis(c(0, 365), norder = 4, nbasis = 7)
  Wfd0 <- fd(matrix(0, wbasis$nbasis, 5), wbasis)
  WfdParobj <- fdPar(Wfd0, 2, 1)

  reg <- register.fd(y0fd, temp_fd, WfdParobj, conv = 1e-4, iterlim = 20, dbglev = 0, periodic = TRUE, crit = 1)
  add_case(
    name = "register_fd_weather_periodic_crit1",
    r_call = "register.fd(y0fd, temp_fd, WfdParobj, conv=1e-4, iterlim=20, dbglev=0, periodic=TRUE, crit=1)",
    input = list(
      basis = basis_info(basis), wbasis = basis_info(wbasis),
      coefs = unname(temp_fd$coefs), y0fd_coefs = unname(y0fd$coefs),
      conv = 1e-4, iterlim = 20, crit = 1, periodic = TRUE
    ),
    output = list(
      regfd_coefs = unname(reg$regfd$coefs),
      warpfd_coefs = unname(reg$warpfd$coefs),
      Wfd_coefs = unname(reg$Wfd$coefs),
      shift = unname(as.numeric(reg$shift))
    ),
    rtol = 1e-5
  )
}

landmarkreg_case()
register_fd_growth_case()
register_fd_periodic_case()

finalize("registration", rtol = 1e-5)
