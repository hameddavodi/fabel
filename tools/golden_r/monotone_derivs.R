# tools/golden_r/monotone_derivs.R
#
# Golden cases for derivatives of monotone and positive smooths: eval.monfd
# (h = integral of exp W, derivatives 0-3) and eval.posfd (exp W, derivatives
# 0-2), plus predict() of a smooth.monotone fit (beta-scaled). Driven by
# tools/make_golden.py, which sources common.R first (defines add_case,
# basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.bspline.basis,
# fd, fdPar, smooth.monotone, eval.monfd, eval.posfd, predict). The growth
# case follows the example of eval.monfd's own help page. No fda source is
# read or copied.
#
# Every case records W's basis and coefficients, so the Python side evaluates
# the same W; the monotone fit's flat direction (W + s, beta_2 e^-s) and its
# iteration path therefore play no role here.
#
# fda API notes (black-box behaviour):
#   - eval.monfd documents derivatives 0-3 only and eval.posfd 0-2 only.
#   - eval.monfd(t, W, 0) integrates exp W numerically; its accuracy is
#     measured by the parity test against fdatools' exact quadrature.

deriv_values <- function(fun, t, Wfd, orders) {
  lapply(orders, function(k) unname(as.matrix(fun(t, Wfd, k))))
}

# ---- (1) a fixed W: two curves on an order-5 B-spline basis --------------

fixed_w_case <- function() {
  basis <- create.bspline.basis(c(0, 2), nbasis = 8, norder = 5)
  coefs <- cbind(
    c(0.3, -0.5, 0.8, 1.2, -0.2, 0.4, 0.9, -0.7),
    c(-1.0, 0.2, 0.6, -0.4, 1.1, 0.5, -0.3, 0.0)
  )
  Wfd <- fd(coefs, basis)
  t <- sort(unique(c(seq(0, 2, length.out = 41), 0.1234, 0.777, 1.618, 1.9)))
  add_case(
    name = "monfd_posfd_fixed_w_two_curves",
    r_call = "eval.monfd(t, Wfd, k), k = 0..3; eval.posfd(t, Wfd, k), k = 0..2",
    input = list(basis = basis_info(basis), coefs = unname(coefs), t = t),
    output = list(
      monfd = deriv_values(eval.monfd, t, Wfd, 0:3),
      posfd = deriv_values(eval.posfd, t, Wfd, 0:2)
    )
  )
}

# ---- (2) growth: smooth.monotone of one girl (eval.monfd help example) ---

growth_case <- function() {
  age <- growth$age
  rng <- c(1, 18)
  norder <- 6
  nbasis <- length(age) + norder - 2
  wbasis <- create.bspline.basis(rng, nbasis, norder, age)
  growfdPar <- fdPar(fd(matrix(0, nbasis, 1), wbasis), 3, 10^(-0.5))
  result <- smooth.monotone(age, growth$hgtf[, 1], growfdPar)
  Wfd <- result$Wfdobj
  beta <- as.numeric(result$beta)
  t <- seq(1, 18, length.out = 69)
  add_case(
    name = "monfd_growth_girl1",
    r_call = "smooth.monotone(age, hgtf[,1], fdPar(wbasis, 3, 10^-0.5)); eval.monfd(t, Wfd, k), predict(result, t, Lfdobj = k), k = 0..3",
    input = list(
      basis = basis_info(wbasis), coefs = unname(Wfd$coefs), beta = beta, t = t
    ),
    output = list(
      monfd = deriv_values(eval.monfd, t, Wfd, 0:3),
      predict = lapply(0:3, function(k) unname(as.numeric(predict(result, t, Lfdobj = k))))
    )
  )
}

fixed_w_case()
growth_case()

finalize("monotone_derivs", rtol = 1e-8)
