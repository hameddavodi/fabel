# tools/golden_r/density.R
#
# Golden cases for the `density` module: penalised log-intensity of a point
# process (intensity.fd) and penalised log-density estimation. Driven by
# tools/make_golden.py, which sources common.R first (defines add_case,
# eval_points, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis, fd,
# fdPar, int2Lfd, intensity.fd, eval.fd). No fda source is read or copied.
#
# fda 6.3.0 black-box facts used here (see PROGRESS.md):
#   - fda 6.3.0 no longer ships density.fd (not exported, not in the
#     namespace, no help page; `intensity.fd`'s help still points to it).
#   - intensity.fd minimises F(c) = -sum_i W(x_i) + int exp W + lambda c'Rc
#     (at c = 0 its reported f is the length of the domain).
#   - When the penalty annihilates constants and the basis spans them, the
#     stationarity condition in the constant direction forces
#     int exp W = n, and the remaining conditions are exactly those of the
#     penalised log-density criterion -sum_i W(x_i) + n log int exp W +
#     lambda c'Rc. So p(x) = exp(W(x)) / n is the penalised density estimate:
#     the density cases below are intensity.fd fits normalised by n.
#   - conv = 1e-12, iterlim = 200 make intensity.fd iterate until its own
#     line search stalls; the remaining error is its quadrature of int exp W.

intensity_case <- function(name, x, basis, lfd, lambda, r_call, density = FALSE,
                           dataset = NULL) {
  nb <- basis$nbasis
  par <- fdPar(fd(matrix(0, nb, 1), basis), lfd, lambda)
  res <- intensity.fd(x, par, conv = 1e-12, iterlim = 200, dbglev = 0)
  t <- eval_points(basis$rangeval)
  w <- as.numeric(eval.fd(t, res$Wfdobj))
  output <- list(
    coefs = as.numeric(res$Wfdobj$coefs),
    f = as.numeric(res$Flist$f),
    grad_norm = as.numeric(res$Flist$norm),
    iternum = res$iternum,
    t = t,
    log_intensity = w,
    intensity = exp(w)
  )
  if (density) {
    output$density <- exp(w) / length(x)
    output$log_density <- w - log(length(x))
  }
  input <- list(
    x = as.numeric(x),
    basis = basis_info(basis),
    lfd = lfd,
    lambda = lambda,
    kind = if (density) "density" else "intensity"
  )
  if (!is.null(dataset)) input$dataset <- dataset
  add_case(name = name, r_call = r_call, input = input, output = output, rtol = 1e-5)
}

# Homogeneous Poisson process, rate 2 (the intensity.fd help example).
N <- 101
tvec <- cumsum(-log(1 - runif(N)) / 2)
intensity_case(
  "intensity_homogeneous_bspline13_L1_lam10",
  tvec, create.bspline.basis(c(0, max(tvec)), 13), 1, 10,
  "intensity.fd(tvec, fdPar(fd(matrix(0,13,1), create.bspline.basis(c(0,max(tvec)),13)), 1, 10), conv=1e-12, iterlim=200)"
)

# Non-homogeneous process mu(t) = 5 + 4 sin(t) on [0, 20], by thinning.
cand <- cumsum(-log(1 - runif(400)) / 9)
cand <- cand[cand < 20]
keep <- runif(length(cand)) < (5 + 4 * sin(cand)) / 9
tsin <- cand[keep]
intensity_case(
  "intensity_sine_bspline23_L2_lam1",
  tsin, create.bspline.basis(c(0, 20), 23), 2, 1,
  "intensity.fd(tsin, fdPar(fd(0, create.bspline.basis(c(0,20),23)), 2, 1), conv=1e-12, iterlim=200)"
)
intensity_case(
  "intensity_sine_bspline11_order3_L1_lam0",
  tsin, create.bspline.basis(c(0, 20), 11, 3), 1, 0,
  "intensity.fd(tsin, fdPar(fd(0, create.bspline.basis(c(0,20),11,3)), 1, 0), conv=1e-12, iterlim=200)"
)
intensity_case(
  "intensity_sine_fourier7_L2_lam0.1",
  tsin, create.fourier.basis(c(0, 20), 7), 2, 0.1,
  "intensity.fd(tsin, fdPar(fd(0, create.fourier.basis(c(0,20),7)), 2, 0.1), conv=1e-12, iterlim=200)"
)

# Real data: June daily precipitation in Regina between 2 and 45 mm
# (Ramsay, Hooker & Graves 2009, section 5.4 / Figure 5.8).
rain <- sort(ReginaPrecip[ReginaPrecip > 2 & ReginaPrecip <= 45])
for (lam in c(0.1, 10)) {
  intensity_case(
    sprintf("density_regina_bspline13_L2_lam%g", lam),
    rain, create.bspline.basis(c(2, 45), 13), 2, lam,
    sprintf("exp(eval.fd(t, intensity.fd(rain, fdPar(fd(0, create.bspline.basis(c(2,45),13)), 2, %g), conv=1e-12, iterlim=200)$Wfdobj)) / length(rain)", lam),
    density = TRUE, dataset = "ReginaPrecip"
  )
}

# Synthetic samples.
z <- rnorm(200)
z <- z[abs(z) < 4]
intensity_case(
  "density_normal_bspline12_order5_L3_lam0.01",
  z, create.bspline.basis(c(-4, 4), 12, 5), 3, 0.01,
  "exp(eval.fd(t, intensity.fd(z, fdPar(fd(0, create.bspline.basis(c(-4,4),12,5)), 3, 0.01), conv=1e-12, iterlim=200)$Wfdobj)) / length(z)",
  density = TRUE
)
g <- rgamma(150, shape = 3, rate = 1)
g <- g[g < 12]
intensity_case(
  "density_gamma_bspline9_L2_lam1",
  g, create.bspline.basis(c(0, 12), 9), 2, 1,
  "exp(eval.fd(t, intensity.fd(g, fdPar(fd(0, create.bspline.basis(c(0,12),9)), 2, 1), conv=1e-12, iterlim=200)$Wfdobj)) / length(g)",
  density = TRUE
)
m <- rnorm(120, 0.5, 1)
m <- m[abs(m) < 3]
intensity_case(
  "density_truncnormal_monomial3_L1_lam0",
  m, create.monomial.basis(c(-3, 3), 3), 1, 0,
  "exp(eval.fd(t, intensity.fd(m, fdPar(fd(0, create.monomial.basis(c(-3,3),3)), 1, 0), conv=1e-12, iterlim=200)$Wfdobj)) / length(m)",
  density = TRUE
)

finalize("density", rtol = 1e-5)
