# tools/golden_r/core.R
#
# Golden cases for the `core` module: fd, eval.fd, inprod, mean.fd, sd.fd,
# center.fd, var.fd, deriv.fd, fd arithmetic, and Lfd-based evaluation.
# Driven by tools/make_golden.py, which sources common.R first (defines
# add_case, eval_points, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (fd, eval.fd, inprod,
# mean.fd, sd.fd, center.fd, var.fd, deriv.fd, arithmetic.fd, int2Lfd,
# vec2Lfd, bifd). No fda source is read or copied.
#
# Decision (recorded also in PROGRESS.md): fda exposes sd.fd, stddev.fd and
# std.fd as equivalent entry points for the same computation; this script
# uses `sd.fd` (matches base R's `sd` generic naming) and records that choice
# in each case's `r_call`.

mk_fd <- function(basisobj, ncurve) {
  fd(matrix(rnorm(basisobj$nbasis * ncurve), basisobj$nbasis, ncurve), basisobj)
}

bsA <- create.bspline.basis(c(0, 1), nbasis = 10, norder = 4)
bsB <- create.bspline.basis(c(0, 1), nbasis = 12, norder = 6)
fourA <- create.fourier.basis(c(0, 1), nbasis = 9, period = 1)
fourW <- create.fourier.basis(c(0, 365), nbasis = 15, period = 365)

fd_bsA_1 <- mk_fd(bsA, 1)
fd_bsA_3 <- mk_fd(bsA, 3)
fd_bsA_3b <- mk_fd(bsA, 3)
fd_bsA_10 <- mk_fd(bsA, 10)
fd_bsB_3 <- mk_fd(bsB, 3)
fd_fourA_3 <- mk_fd(fourA, 3)
fd_fourW_3 <- mk_fd(fourW, 3)

t_bsA <- eval_points(c(0, 1))
t_fourW <- eval_points(c(0, 365))

# ---- eval.fd, derivs 0,1,2 -------------------------------------------------
for (d in 0:2) {
  Vb <- eval.fd(t_bsA, fd_bsA_3, d)
  add_case(
    name = sprintf("evalfd_bspline_n3curves_d%d", d),
    r_call = sprintf("eval.fd(t, fdobj, %d)  # fdobj on 10-basis order-4 bspline, 3 curves", d),
    input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_3$coefs), t = t_bsA, deriv = d),
    output = list(values = unname(Vb))
  )
  Vf <- eval.fd(t_bsA, fd_fourA_3, d)
  add_case(
    name = sprintf("evalfd_fourier_n3curves_d%d", d),
    r_call = sprintf("eval.fd(t, fdobj, %d)  # fdobj on 9-basis fourier, 3 curves", d),
    input = list(basis = basis_info(fourA), coefs = unname(fd_fourA_3$coefs), t = t_bsA, deriv = d),
    output = list(values = unname(Vf))
  )
}

# ---- inprod: fd x fd, bspline, various Lfd on each side --------------------
lfd_pairs <- list(c(0, 0), c(1, 1), c(2, 0), c(0, 2))
for (pair in lfd_pairs) {
  M <- inprod(fd_bsA_3, fd_bsA_3, int2Lfd(pair[1]), int2Lfd(pair[2]))
  add_case(
    name = sprintf("inprod_bspline_L%d_%d", pair[1], pair[2]),
    r_call = sprintf("inprod(fdobj, fdobj, int2Lfd(%d), int2Lfd(%d))", pair[1], pair[2]),
    input = list(
      basis = basis_info(bsA), coefs = unname(fd_bsA_3$coefs),
      lfd1 = pair[1], lfd2 = pair[2]
    ),
    output = list(values = unname(M))
  )
}

# ---- inprod: fd x fd, fourier -----------------------------------------------
lfd_pairs_f <- list(c(0, 0), c(1, 1), c(2, 2))
for (pair in lfd_pairs_f) {
  M <- inprod(fd_fourA_3, fd_fourA_3, int2Lfd(pair[1]), int2Lfd(pair[2]))
  add_case(
    name = sprintf("inprod_fourier_L%d_%d", pair[1], pair[2]),
    r_call = sprintf("inprod(fdobj, fdobj, int2Lfd(%d), int2Lfd(%d))", pair[1], pair[2]),
    input = list(
      basis = basis_info(fourA), coefs = unname(fd_fourA_3$coefs),
      lfd1 = pair[1], lfd2 = pair[2]
    ),
    output = list(values = unname(M))
  )
}

# ---- inprod: basis x basis, mixed types (bspline x fourier, same domain) --
Mmix <- inprod(bsA, fourA)
add_case(
  name = "inprod_basis_bspline_x_fourier",
  r_call = "inprod(bspline_basis, fourier_basis)",
  input = list(basis1 = basis_info(bsA), basis2 = basis_info(fourA)),
  output = list(values = unname(Mmix))
)

# ---- mean.fd / sd.fd / center.fd -------------------------------------------
mfd <- mean.fd(fd_bsA_10)
add_case(
  name = "mean_fd_bspline_n10curves",
  r_call = "mean.fd(fdobj)  # fdobj: 10-basis order-4 bspline, 10 curves",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_10$coefs)),
  output = list(coefs = unname(mfd$coefs), basis = basis_info(mfd$basis))
)

sdfd <- sd.fd(fd_bsA_10)
add_case(
  name = "sd_fd_bspline_n10curves",
  r_call = "sd.fd(fdobj)  # fdobj: 10-basis order-4 bspline, 10 curves",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_10$coefs)),
  output = list(coefs = unname(sdfd$coefs), basis = basis_info(sdfd$basis))
)

cfd <- center.fd(fd_bsA_10)
add_case(
  name = "center_fd_bspline_n10curves",
  r_call = "center.fd(fdobj)  # fdobj: 10-basis order-4 bspline, 10 curves",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_10$coefs)),
  output = list(coefs = unname(cfd$coefs), basis = basis_info(cfd$basis))
)

# ---- var.fd -> bifd ----------------------------------------------------------
vfd <- var.fd(fd_bsA_10)
add_case(
  name = "var_fd_bspline_n10curves",
  r_call = "var.fd(fdobj)  # fdobj: 10-basis order-4 bspline, 10 curves",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_10$coefs)),
  output = list(
    coefs = unname(vfd$coefs),
    sbasis = basis_info(vfd$sbasis),
    tbasis = basis_info(vfd$tbasis)
  )
)

# ---- deriv.fd ---------------------------------------------------------------
dfd6 <- deriv.fd(fd_bsB_3, int2Lfd(2))
add_case(
  name = "deriv_fd_bspline_order6_L2",
  r_call = "deriv.fd(fdobj, int2Lfd(2))  # fdobj: 12-basis order-6 bspline, 3 curves",
  input = list(basis = basis_info(bsB), coefs = unname(fd_bsB_3$coefs), lfd = 2),
  output = list(coefs = unname(dfd6$coefs), basis = basis_info(dfd6$basis))
)

dfdF <- deriv.fd(fd_fourA_3, int2Lfd(1))
add_case(
  name = "deriv_fd_fourier_L1",
  r_call = "deriv.fd(fdobj, int2Lfd(1))  # fdobj: 9-basis fourier, 3 curves",
  input = list(basis = basis_info(fourA), coefs = unname(fd_fourA_3$coefs), lfd = 1),
  output = list(coefs = unname(dfdF$coefs), basis = basis_info(dfdF$basis))
)

# ---- arithmetic ---------------------------------------------------------------
s <- fd_bsA_3 + fd_bsA_3b
add_case(
  name = "fd_add_bspline",
  r_call = "fd1 + fd2  # both on 10-basis order-4 bspline, 3 curves",
  input = list(
    basis = basis_info(bsA), coefs1 = unname(fd_bsA_3$coefs), coefs2 = unname(fd_bsA_3b$coefs)
  ),
  output = list(coefs = unname(s$coefs), basis = basis_info(s$basis))
)

p <- fd_bsA_3 * fd_bsA_3b
add_case(
  name = "fd_mul_bspline",
  r_call = "fd1 * fd2  # both on 10-basis order-4 bspline, 3 curves",
  input = list(
    basis = basis_info(bsA), coefs1 = unname(fd_bsA_3$coefs), coefs2 = unname(fd_bsA_3b$coefs)
  ),
  output = list(coefs = unname(p$coefs), basis = basis_info(p$basis))
)

s2 <- fd_bsA_3 * 2
add_case(
  name = "fd_scalar_mul",
  r_call = "fd1 * 2",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_3$coefs)),
  output = list(coefs = unname(s2$coefs), basis = basis_info(s2$basis))
)

sq <- fd_bsA_3^2
add_case(
  name = "fd_power2",
  r_call = "fd1 ^ 2",
  input = list(basis = basis_info(bsA), coefs = unname(fd_bsA_3$coefs)),
  output = list(coefs = unname(sq$coefs), basis = basis_info(sq$basis))
)

# ---- Lfd-based evaluation -----------------------------------------------------
V_int2 <- eval.fd(t_bsA, fd_bsB_3, int2Lfd(2))
add_case(
  name = "evalfd_lfd_int2lfd2_bspline_order6",
  r_call = "eval.fd(t, fdobj, int2Lfd(2))  # fdobj: 12-basis order-6 bspline, 3 curves",
  input = list(basis = basis_info(bsB), coefs = unname(fd_bsB_3$coefs), t = t_bsA, lfd = "int2Lfd(2)"),
  output = list(values = unname(V_int2))
)

w1 <- 0.7
Lvec <- vec2Lfd(c(0, w1), c(0, 1))
V_vec <- eval.fd(t_bsA, fd_bsA_3, Lvec)
add_case(
  name = "evalfd_lfd_vec2lfd_bspline",
  r_call = "eval.fd(t, fdobj, vec2Lfd(c(0, 0.7), c(0,1)))  # fdobj: 10-basis order-4 bspline, 3 curves",
  input = list(
    basis = basis_info(bsA), coefs = unname(fd_bsA_3$coefs), t = t_bsA,
    lfd_weights = c(0, w1), lfd_rangeval = c(0, 1)
  ),
  output = list(values = unname(V_vec))
)

harm_weights <- c(0, (2 * pi / 365)^2, 0)
harmLfd <- vec2Lfd(harm_weights, c(0, 365))
V_harm <- eval.fd(t_fourW, fd_fourW_3, harmLfd)
add_case(
  name = "evalfd_lfd_harmonic_accelerator_fourier",
  r_call = "eval.fd(t, fdobj, vec2Lfd(c(0, (2*pi/365)^2, 0), c(0,365)))  # fdobj: 15-basis fourier, 3 curves",
  input = list(
    basis = basis_info(fourW), coefs = unname(fd_fourW_3$coefs), t = t_fourW,
    lfd_weights = harm_weights, lfd_rangeval = c(0, 365)
  ),
  output = list(values = unname(V_harm))
)

# ---- bifd: explicit construction from a multi-dim coefficient array ----------
# Note: fda 6.3.0's bifd() is unconditionally broken for a 3-D coef array
# (ndim == 3, i.e. a "reps" dimension with no "vars" dimension): its body sets
# `defaultnames` only inside `if (ndim == 2)` / `if (ndim == 4)` branches (a
# copy-paste bug -- ndim == 3 is never handled), but then unconditionally runs
# `names(defaultnames) <- c(...)`, so a 3-D array throws "object 'defaultnames'
# not found" no matter what -- even when `fdnames` is supplied explicitly,
# since that line does not reference the `fdnames` argument at all. Worked
# around by using a 4-D array (dim = c(nbasis_s, nbasis_t, nrep, nvar=1))
# instead, which hits the working ndim == 4 branch and still exercises the
# reps dimension the original 3-D case intended. Recorded in PROGRESS.md.
sbf <- create.bspline.basis(c(0, 1), nbasis = 4, norder = 3)
tbf <- create.bspline.basis(c(0, 1), nbasis = 3, norder = 2)
coef4d <- array(rnorm(4 * 3 * 2 * 1), dim = c(4, 3, 2, 1))
bf <- bifd(coef = coef4d, sbasisobj = sbf, tbasisobj = tbf)
add_case(
  name = "bifd_construct_4d_coef",
  r_call = "bifd(coef=array(rnorm(24), dim=c(4,3,2,1)), sbasisobj, tbasisobj)",
  input = list(
    sbasis = basis_info(sbf), tbasis = basis_info(tbf),
    coef = coef4d, coef_dim = dim(coef4d)
  ),
  output = list(coefs = bf$coefs, sbasis = basis_info(bf$sbasis), tbasis = basis_info(bf$tbasis))
)

finalize("core", rtol = 1e-6)
