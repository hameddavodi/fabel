# tools/golden_r/profiling.R
#
# Golden cases for the `profiling` module: generalized profiling estimation
# of ODE parameters, checked on R's continuously stirred tank reactor (CSTR)
# functions. Driven by tools/make_golden.py, which sources common.R first
# (defines add_case, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (CSTR2in, CSTR2,
# CSTRfitLS, CSTRfn, CSTRsse, CSTRres, quadset, create.bspline.basis,
# eval.basis, smooth.basis, fdPar) plus deSolve::lsoda and stats::nls. No fda
# source is read or copied; every fact below was found by calling the API.
#
# fda 6.3.0 behaviour found while writing this script (see PROGRESS.md):
#   - CSTRfitLS returns res$Lres with ONE column ("LConc", nquad rows): the
#     concentration equation residual only. The temperature equation residual
#     is missing, although Dres$DLres has 2 * nquad rows (both equations).
#   - Consequently CSTRfn, CSTRsse and CSTRres stop with "non-conformable
#     arguments" for every configuration tried (fit = c(1,1), c(0,1), c(1,0);
#     coef0 as matrix, vector or data.frame; scalar or vector lambda; gradwrd
#     TRUE or FALSE; observations at the quadrature points). The golden file
#     records the error messages (cstr_fn_errors).
#   - The reference inner fits and profile estimates below are therefore built
#     from R's own pieces: the CSTRfitLS data residuals and full Jacobian, the
#     concentration equation residual it does return, and the temperature
#     equation residual computed from CSTR2 (R's right-hand side) at the
#     quadrature points. The coefficients minimise that criterion by
#     Gauss-Newton; the parameters minimise the profiled data criterion by
#     stats::nls (central-difference derivatives, tight tolerance).
#   - CSTR2 right-hand side, found by probing each constant: F enters as F/V
#     in both equations, alpha = a Fc^(b+1) / (V rho Cp (Fc + a Fc^b /
#     (2 rhoc Cpc))), and the reaction heat term is (-delH / (V rho Cp)) k(T) C.
#   - CSTR2in step scenarios change inputs on [4k, 4k + 4), left-closed.

TLIM <- 24
COND <- "all.cool.step"
PAR_NAMES <- c("kref", "EoverR", "a", "b")

base_fitstruct <- function(theta = c(0.4610, 0.83301, 1.678, 0.5)) {
  fs <- list(V = 1.0, Cp = 1.0, rho = 1.0, delH = -130.0, Cpc = 1.0, rhoc = 1.0,
             Tref = 350)
  fs[PAR_NAMES] <- as.list(theta)
  fs
}

# ---- (a) quadrature: quadset ----------------------------------------------

quadset_cases <- function() {
  for (spec in list(list(n = 5, breaks = c(0, 0.3, 1)), list(n = 7, breaks = c(0, 0.3, 1)),
                    list(n = 5, breaks = seq(0, 24, by = 0.5)))) {
    qv <- quadset(nquad = spec$n, breaks = spec$breaks)
    add_case(
      name = sprintf("quadset_n%d_%dbreaks", spec$n, length(spec$breaks)),
      r_call = sprintf("quadset(nquad = %d, breaks = breaks)", spec$n),
      input = list(n_quad = spec$n, breaks = spec$breaks),
      output = list(nodes = unname(qv[, 1]), weights = unname(qv[, 2]))
    )
  }
}

# ---- (b) inputs: CSTR2in --------------------------------------------------

cstr2in_cases <- function() {
  t <- sort(unique(c(seq(-1, 66, by = 1 / 8), 0:64, 1.999999, 2.000001,
                     11.999999, 12.000001)))
  for (cond in c("all.cool.step", "all.hot.step", "Tc.hot.step", "Tc.cool.step")) {
    add_case(
      name = paste0("cstr2in_", gsub(".", "_", cond, fixed = TRUE)),
      r_call = sprintf("CSTR2in(t, '%s')", cond),
      input = list(t = t, condition = cond),
      output = list(inputs = unname(CSTR2in(t, cond)))
    )
  }
}

# ---- (c) right-hand side: CSTR2 -------------------------------------------

cstr2_cases <- function() {
  t <- c(0.5, 1, 5, 9, 17, 21, 29, 33, 41, 45, 53, 57, 61, 63.5)
  set.seed(GOLDEN_SEED + 1)
  X <- cbind(Conc = runif(length(t), 0.3, 2.0), Temp = runif(length(t), 320, 380))
  variants <- list(
    list(name = "default", cond = "all.cool.step", theta = c(0.4610, 0.83301, 1.678, 0.5),
         const = list()),
    list(name = "hot_other_theta", cond = "all.hot.step", theta = c(0.3, 1.1, 2.2, 0.7),
         const = list()),
    list(name = "tc_step_constants", cond = "Tc.hot.step", theta = c(0.5, 0.7, 1.4, 0.3),
         const = list(V = 2, Cp = 1.3, rho = 0.8, delH = -110, Cpc = 1.1, rhoc = 0.9,
                      Tref = 340))
  )
  for (v in variants) {
    fs <- base_fitstruct(v$theta)
    fs[names(v$const)] <- v$const
    parms <- list(fitstruct = fs, condition = v$cond, Tlim = 64)
    add_case(
      name = paste0("cstr2_", v$name),
      r_call = sprintf("CSTR2(t, X, list(fitstruct = fs, condition = '%s'))", v$cond),
      input = list(t = t, x = unname(X), condition = v$cond, theta = v$theta,
                   constants = if (length(v$const)) v$const else list(V = 1)),
      output = list(dx = unname(CSTR2(t, X, parms)[[1]]))
    )
  }
}

# ---- (d) simulated experiment ---------------------------------------------

set.seed(GOLDEN_SEED)
tobs <- seq(0, TLIM, by = 1 / 4)
true_fs <- base_fitstruct()
sol <- deSolve::lsoda(y = c(Conc = 1.5965, Temp = 341.3754), times = tobs, func = CSTR2,
                      parms = list(fitstruct = true_fs, condition = COND, Tlim = TLIM),
                      rtol = 1e-10, atol = 1e-10)
ytrue <- unname(sol[, 2:3])
yobs <- ytrue + cbind(rnorm(length(tobs), 0, 0.02), rnorm(length(tobs), 0, 0.5))
cbasis <- quadset(nquad = 5, basisobj = create.bspline.basis(c(0, TLIM), norder = 4,
                                                             breaks = seq(0, TLIM, by = 0.5)))
quadpts <- cbasis$quadvals[, 1]
quadwts <- cbasis$quadvals[, 2]
qin <- CSTR2in(quadpts, COND)
nb <- cbasis$nbasis
Cwt <- as.numeric(var(yobs[, 1]))
Twt <- as.numeric(var(yobs[, 2]))
Dq <- eval.basis(quadpts, cbasis, 1)
smooth_coef <- function(j) smooth.basis(tobs, yobs[, j], fdPar(cbasis, 2, 1e-4))$fd$coefs[, 1]
coef_smooth <- cbind(smooth_coef(1), smooth_coef(2))

make_datstruct <- function(fit) {
  list(basismat = eval.basis(tobs, cbasis), Dbasismat = eval.basis(tobs, cbasis, 1),
       Cwt = Cwt, Twt = Twt, y = yobs,
       quadbasismat = eval.basis(quadpts, cbasis), Dquadbasismat = Dq,
       F. = qin[, "F."], CA0 = qin[, "CA0"], T0 = qin[, "T0"], Tc = qin[, "Tcin"],
       Fc = qin[, "Fc"], quadpts = quadpts, quadwts = quadwts)
}

make_fitstruct <- function(theta, fit, coef0) {
  fs <- base_fitstruct(theta)
  fs$Tcin <- qin[, "Tcin"]
  fs$fit <- fit
  fs$estimate <- c(1, 1, 1, 1)
  fs$coef0 <- coef0
  fs
}

add_case(
  name = "cstr_data",
  r_call = "lsoda(c(Conc=1.5965, Temp=341.3754), tobs, CSTR2, parms) + noise; quadset(5, basis)",
  input = list(condition = COND, tlim = TLIM, theta_true = unname(unlist(true_fs[PAR_NAMES])),
               sd = c(0.02, 0.5)),
  output = list(t = tobs, y = unname(yobs), ytrue = ytrue, basis = basis_info(cbasis),
                quad_nodes = unname(quadpts), quad_weights = unname(quadwts),
                Cwt = Cwt, Twt = Twt, coef_smooth = unname(coef_smooth))
)

# ---- (e) residuals and Jacobian: CSTRfitLS --------------------------------

fitls <- function(coef, theta, fit, lambda) {
  CSTRfitLS(coef, make_datstruct(fit), make_fitstruct(theta, fit, coef), lambda, gradwrd = TRUE)
}

fitls_cases <- function() {
  set.seed(GOLDEN_SEED + 2)
  bumped <- coef_smooth + cbind(rnorm(nb, 0, 0.01), rnorm(nb, 0, 0.3))
  variants <- list(
    list(name = "fit11_truth", coef = coef_smooth, theta = c(0.4610, 0.83301, 1.678, 0.5),
         fit = c(1, 1), lambda = c(100, 100)),
    list(name = "fit01_other", coef = bumped, theta = c(0.4, 0.8, 1.6, 0.45),
         fit = c(0, 1), lambda = c(3, 7)),
    list(name = "fit10_other", coef = bumped, theta = c(0.5, 0.9, 1.9, 0.6),
         fit = c(1, 0), lambda = c(10, 0.5))
  )
  for (v in variants) {
    r <- fitls(v$coef, v$theta, v$fit, v$lambda)
    add_case(
      name = paste0("cstr_fitls_", v$name),
      r_call = "CSTRfitLS(coef, datstruct, fitstruct, lambda, gradwrd = TRUE)",
      input = list(coef = unname(v$coef), theta = v$theta, fit = v$fit, lambda = v$lambda,
                   source_case = "cstr_data"),
      output = list(Sres = unname(r$res$Sres), Lres = unname(r$res$Lres),
                    DSres = unname(r$Dres$DSres), DLres = unname(r$Dres$DLres))
    )
  }
}

# ---- (f) CSTRfn / CSTRsse / CSTRres fail in fda 6.3.0 ---------------------

errors_case <- function() {
  ds <- make_datstruct(c(1, 1))
  fs <- make_fitstruct(c(0.4610, 0.83301, 1.678, 0.5), c(1, 1), coef_smooth)
  msg <- function(expr) tryCatch({ force(expr); "no error" }, error = function(e) conditionMessage(e))
  add_case(
    name = "cstr_fn_errors",
    r_call = "CSTRfn(parvec, ...); CSTRsse(parvec, ...); nls(~ CSTRres(kref, EoverR, ...))",
    input = list(parvec = c(0.4610, 0.83301, 1.678, 0.5), lambda = c(100, 100),
                 source_case = "cstr_data"),
    output = list(
      CSTRfn = msg(CSTRfn(c(0.4610, 0.83301, 1.678, 0.5), ds, fs, cbasis, c(100, 100),
                          gradwrd = TRUE)),
      CSTRfn_gradwrd_false = msg(CSTRfn(c(0.4610, 0.83301, 1.678, 0.5), ds, fs, cbasis,
                                        c(100, 100), gradwrd = FALSE)),
      CSTRsse = msg(CSTRsse(c(0.4610, 0.83301, 1.678, 0.5), ds, fs, cbasis, c(100, 100))),
      CSTRres = msg(CSTRres(kref = 0.4610, EoverR = 0.83301, datstruct = ds,
                            fitstruct = within(fs, estimate <- c(1, 1, 0, 0)),
                            CSTRbasis = cbasis, lambda = c(100, 100))),
      Lres_columns = ncol(fitls(coef_smooth, c(0.4610, 0.83301, 1.678, 0.5), c(1, 1),
                                c(100, 100))$res$Lres)
    )
  )
}

# ---- (g) reference inner fits from R's pieces -----------------------------

# Full residual vector and Jacobian of the inner criterion at (coef, theta):
# CSTRfitLS's data residuals, its concentration equation residual, and the
# temperature equation residual from CSTR2; the Jacobian is CSTRfitLS's.
full_residuals <- function(coef, theta, fit, lambda) {
  r <- fitls(coef, theta, fit, lambda)
  parms <- list(fitstruct = base_fitstruct(theta), condition = COND, Tlim = TLIM)
  X <- eval.basis(quadpts, cbasis) %*% coef
  colnames(X) <- c("Conc", "Temp")
  fT <- CSTR2(quadpts, X, parms)[[1]][, 2]
  LT <- sqrt(lambda[2] * quadwts / Twt) * (as.numeric(Dq %*% coef[, 2]) - fT)
  list(res = c(as.numeric(r$res$Sres), as.numeric(r$res$Lres[, 1]), LT),
       jac = rbind(r$Dres$DSres, r$Dres$DLres),
       nS = length(r$res$Sres))
}

inner_fit <- function(theta, fit, lambda, coef0, maxit = 200) {
  coef <- coef0
  fr <- full_residuals(coef, theta, fit, lambda)
  val <- sum(fr$res^2)
  prev_rel <- Inf
  for (it in 1:maxit) {
    J <- fr$jac
    s <- 1 / sqrt(colSums(J^2))
    Js <- sweep(J, 2, s, "*")
    step <- -s * qr.solve(Js, fr$res)
    stepm <- matrix(step, ncol = 2)
    rel <- max(apply(abs(stepm), 2, max) / apply(abs(coef), 2, max))
    alpha <- 1
    # step halving only while far from the minimum; once the relative step is
    # below 1e-6 the criterion changes at rounding level and full Gauss-Newton
    # steps are taken until every coefficient has settled
    repeat {
      trial <- coef + alpha * matrix(step, ncol = 2)
      ft <- full_residuals(trial, theta, fit, lambda)
      if (rel <= 1e-6 || sum(ft$res^2) <= val || alpha < 1e-10) break
      alpha <- alpha / 2
    }
    coef <- trial
    fr <- ft
    val <- sum(ft$res^2)
    # stop once the step is at rounding level, or has stopped shrinking there
    if (alpha * rel <= 1e-14 || (rel <= 1e-12 && rel > prev_rel / 2)) break
    prev_rel <- rel
  }
  list(coef = coef, criterion = val, sse = sum(fr$res[seq_len(fr$nS)]^2),
       data_res = fr$res[seq_len(fr$nS)], n_iter = it)
}

start_coef <- function(fit) {
  if (fit[1] == 1) coef_smooth else cbind(rep(1.5, nb), coef_smooth[, 2])
}

inner_cases <- function() {
  variants <- list(
    list(name = "fit11_truth", theta = c(0.4610, 0.83301, 1.678, 0.5), fit = c(1, 1),
         lambda = c(100, 100)),
    list(name = "fit01_other", theta = c(0.4, 0.8, 1.6, 0.45), fit = c(0, 1),
         lambda = c(100, 100)),
    list(name = "fit11_small_lambda", theta = c(0.5, 0.9, 1.9, 0.6), fit = c(1, 1),
         lambda = c(3, 7))
  )
  for (v in variants) {
    inn <- inner_fit(v$theta, v$fit, v$lambda, start_coef(v$fit))
    add_case(
      name = paste0("cstr_inner_", v$name),
      r_call = "Gauss-Newton on c(Sres, Lres, LT) with rbind(DSres, DLres) from CSTRfitLS; LT from CSTR2",
      input = list(theta = v$theta, fit = v$fit, lambda = v$lambda, source_case = "cstr_data"),
      # the Gauss-Newton fixed point is exact to rounding, so the module rtol applies
      output = list(coef = unname(inn$coef), criterion = inn$criterion, sse = inn$sse)
    )
  }
}

# ---- (h) reference profile estimates by nls -------------------------------

profile_cases <- function() {
  variants <- list(
    list(name = "fit11_kref_EoverR", fit = c(1, 1), lambda = c(100, 100),
         est = c("kref", "EoverR"), start = c(0.4, 0.8)),
    list(name = "fit01_kref_EoverR", fit = c(0, 1), lambda = c(100, 100),
         est = c("kref", "EoverR"), start = c(0.4, 0.8)),
    # a and b together are not identifiable (singular gradient at any start)
    list(name = "fit11_kref_EoverR_a", fit = c(1, 1), lambda = c(100, 100),
         est = c("kref", "EoverR", "a"), start = c(0.4, 0.8, 1.6))
  )
  for (v in variants) {
    theta_of <- function(par) {
      th <- unlist(true_fs[PAR_NAMES])
      th[v$est] <- par
      unname(th)
    }
    # `idx` only gives nls a data variable; the residuals come from the inner fit
    profiled <- function(idx, ...) {
      inn <- inner_fit(theta_of(c(...)), v$fit, v$lambda, start_coef(v$fit))
      inn$data_res
    }
    rhs <- as.call(c(list(as.name("profiled"), as.name("idx")), lapply(v$est, as.name)))
    start <- setNames(as.list(v$start), v$est)
    n_res <- sum(v$fit) * length(tobs)
    fit <- nls(as.formula(call("~", rhs)), data = list(idx = seq_len(n_res)), start = start,
               control = nls.control(maxiter = 100, tol = 1e-6, nDcentral = TRUE,
                                     scaleOffset = 0))
    est <- unname(coef(fit))
    inn <- inner_fit(theta_of(est), v$fit, v$lambda, start_coef(v$fit))
    add_case(
      name = paste0("cstr_profile_", v$name),
      r_call = sprintf("nls(~ profiled(idx, %s), start, nls.control(tol = 1e-6, nDcentral = TRUE))",
                       paste(v$est, collapse = ", ")),
      input = list(fit = v$fit, lambda = v$lambda, estimate = v$est, start = v$start,
                   source_case = "cstr_data"),
      output = list(theta = est, vcov = unname(vcov(fit)), sse = sum(resid(fit)^2),
                    coef = unname(inn$coef)),
      # nls stops at relative offset 1e-6 (its floor on this problem is ~1e-7,
      # set by the rounding noise of the profiled residuals)
      rtol = 1e-5
    )
  }
}

quadset_cases()
cstr2in_cases()
cstr2_cases()
fitls_cases()
errors_case()
inner_cases()
profile_cases()

finalize("profiling", rtol = 1e-8)
