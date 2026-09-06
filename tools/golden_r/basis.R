# tools/golden_r/basis.R
#
# Golden cases for the `basis` module: eval.basis / eval.penalty for every
# fda basis type (bspline, fourier, monomial, exponential, power, constant,
# polygonal). Driven by tools/make_golden.py, which sources common.R first
# (defines add_case, eval_points, basis_info, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API is called (create.*.basis,
# eval.basis, eval.penalty). No fda source is read or copied.

bspline_cases <- function() {
  orders <- c(1, 2, 3, 4, 6)
  nbasis_opts <- c(4, 7, 10, 15, 23)
  dom <- c(0, 1)
  t <- eval_points(dom)
  for (norder in orders) {
    for (nb in nbasis_opts) {
      if (nb < norder) next
      b <- create.bspline.basis(rangeval = dom, nbasis = nb, norder = norder)
      derivs <- Filter(function(d) d < norder, c(0, 1, 2))
      for (d in derivs) {
        V <- eval.basis(t, b, d)
        add_case(
          name = sprintf("bspline_eval_k%d_n%d_dom0_1_d%d", norder, nb, d),
          r_call = sprintf(
            "eval.basis(t, create.bspline.basis(c(0,1), %d, %d), %d)", nb, norder, d
          ),
          input = list(domain = dom, n_basis = nb, order = norder, t = t, deriv = d),
          output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
        )
      }
      # eval.penalty needs norder - Lfdobj >= 2 for B-splines (empirically verified
      # against R fda 6.3.0; the naive "Lfdobj < norder" bound errors for e.g.
      # order 3, Lfdobj 2). Recorded in PROGRESS.md.
      lfds <- Filter(function(l) l <= norder - 2, c(0, 1, 2))
      for (L in lfds) {
        P <- eval.penalty(b, int2Lfd(L))
        add_case(
          name = sprintf("bspline_penalty_k%d_n%d_dom0_1_L%d", norder, nb, L),
          r_call = sprintf(
            "eval.penalty(create.bspline.basis(c(0,1), %d, %d), int2Lfd(%d))", nb, norder, L
          ),
          input = list(domain = dom, n_basis = nb, order = norder, lfd = L),
          output = list(penalty = unname(P), nbasis = b$nbasis)
        )
      }
    }
  }
}

bspline_extra_domains <- function() {
  domains <- list(c(0, 365), c(-1, 2))
  for (dom in domains) {
    norder <- 4
    nb <- 23
    b <- create.bspline.basis(rangeval = dom, nbasis = nb, norder = norder)
    t <- eval_points(dom)
    for (d in c(0, 1, 2)) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("bspline_eval_k%d_n%d_dom%g_%g_d%d", norder, nb, dom[1], dom[2], d),
        r_call = sprintf(
          "eval.basis(t, create.bspline.basis(c(%g,%g), %d, %d), %d)",
          dom[1], dom[2], nb, norder, d
        ),
        input = list(domain = dom, n_basis = nb, order = norder, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in c(0, 1, 2)) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("bspline_penalty_k%d_n%d_dom%g_%g_L%d", norder, nb, dom[1], dom[2], L),
        r_call = sprintf(
          "eval.penalty(create.bspline.basis(c(%g,%g), %d, %d), int2Lfd(%d))",
          dom[1], dom[2], nb, norder, L
        ),
        input = list(domain = dom, n_basis = nb, order = norder, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }
}

bspline_explicit_breaks <- function() {
  breaks_list <- list(
    c(0, 0.1, 0.15, 0.4, 0.42, 1),
    c(0, 50, 120, 200, 365),
    c(-1, -0.5, 0, 0.3, 1, 2)
  )
  for (i in seq_along(breaks_list)) {
    breaks <- breaks_list[[i]]
    dom <- range(breaks)
    norder <- 4
    b <- create.bspline.basis(rangeval = dom, breaks = breaks, norder = norder)
    t <- eval_points(dom)
    for (d in c(0, 1, 2)) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("bspline_breaks%d_eval_d%d", i, d),
        r_call = sprintf(
          "eval.basis(t, create.bspline.basis(breaks=c(%s), norder=%d), %d)",
          paste(breaks, collapse = ","), norder, d
        ),
        input = list(domain = dom, breaks = breaks, order = norder, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in c(0, 1, 2)) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("bspline_breaks%d_penalty_L%d", i, L),
        r_call = sprintf(
          "eval.penalty(create.bspline.basis(breaks=c(%s), norder=%d), int2Lfd(%d))",
          paste(breaks, collapse = ","), norder, L
        ),
        input = list(domain = dom, breaks = breaks, order = norder, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }
}

fourier_cases <- function() {
  nbasis_opts <- c(3, 5, 9, 65)
  dom <- c(0, 1)
  t <- eval_points(dom)
  for (nb in nbasis_opts) {
    b <- create.fourier.basis(rangeval = dom, nbasis = nb, period = diff(dom))
    for (d in 0:3) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("fourier_eval_n%d_d%d", nb, d),
        r_call = sprintf(
          "eval.basis(t, create.fourier.basis(c(0,1), %d, 1), %d)", nb, d
        ),
        input = list(domain = dom, n_basis = nb, period = diff(dom), t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis)
      )
    }
    for (L in 0:2) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("fourier_penalty_n%d_L%d", nb, L),
        r_call = sprintf(
          "eval.penalty(create.fourier.basis(c(0,1), %d, 1), int2Lfd(%d))", nb, L
        ),
        input = list(domain = dom, n_basis = nb, period = diff(dom), lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }

  # period != domain width
  b_pm <- create.fourier.basis(rangeval = dom, nbasis = 9, period = 2)
  for (d in 0:3) {
    V <- eval.basis(t, b_pm, d)
    add_case(
      name = sprintf("fourier_periodmismatch_eval_d%d", d),
      r_call = sprintf("eval.basis(t, create.fourier.basis(c(0,1), 9, 2), %d)", d),
      input = list(domain = dom, n_basis = 9, period = 2, t = t, deriv = d),
      output = list(values = unname(V), nbasis = b_pm$nbasis)
    )
  }
  for (L in 0:2) {
    P <- eval.penalty(b_pm, int2Lfd(L))
    add_case(
      name = sprintf("fourier_periodmismatch_penalty_L%d", L),
      r_call = sprintf("eval.penalty(create.fourier.basis(c(0,1), 9, 2), int2Lfd(%d))", L),
      input = list(domain = dom, n_basis = 9, period = 2, lfd = L),
      output = list(penalty = unname(P), nbasis = b_pm$nbasis)
    )
  }

  # even nbasis request -> fda forces odd; record actual nbasis returned
  b_even <- create.fourier.basis(rangeval = dom, nbasis = 4, period = 1)
  V_even <- eval.basis(t, b_even, 0)
  add_case(
    name = "fourier_even_nbasis_request_eval_d0",
    r_call = "eval.basis(t, create.fourier.basis(c(0,1), 4, 1), 0)",
    input = list(domain = dom, n_basis = 4, period = 1, t = t, deriv = 0),
    output = list(values = unname(V_even), nbasis = b_even$nbasis)
  )

  # canonical (0, 365) domain, period 365
  dom365 <- c(0, 365)
  t365 <- eval_points(dom365)
  b365 <- create.fourier.basis(rangeval = dom365, nbasis = 65, period = 365)
  for (d in 0:3) {
    V <- eval.basis(t365, b365, d)
    add_case(
      name = sprintf("fourier_eval_dom365_n65_d%d", d),
      r_call = sprintf("eval.basis(t, create.fourier.basis(c(0,365), 65, 365), %d)", d),
      input = list(domain = dom365, n_basis = 65, period = 365, t = t365, deriv = d),
      output = list(values = unname(V), nbasis = b365$nbasis)
    )
  }
  for (L in 0:2) {
    P <- eval.penalty(b365, int2Lfd(L))
    add_case(
      name = sprintf("fourier_penalty_dom365_n65_L%d", L),
      r_call = sprintf("eval.penalty(create.fourier.basis(c(0,365), 65, 365), int2Lfd(%d))", L),
      input = list(domain = dom365, n_basis = 65, period = 365, lfd = L),
      output = list(penalty = unname(P), nbasis = b365$nbasis)
    )
  }
}

monomial_cases <- function() {
  dom <- c(0, 1)
  t <- eval_points(dom)
  for (nb in 1:6) {
    b <- create.monomial.basis(rangeval = dom, nbasis = nb)
    for (d in 0:2) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("monomial_eval_n%d_d%d", nb, d),
        r_call = sprintf("eval.basis(t, create.monomial.basis(c(0,1), %d), %d)", nb, d),
        input = list(domain = dom, n_basis = nb, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in 0:2) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("monomial_penalty_n%d_L%d", nb, L),
        r_call = sprintf("eval.penalty(create.monomial.basis(c(0,1), %d), int2Lfd(%d))", nb, L),
        input = list(domain = dom, n_basis = nb, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }

  # custom exponents
  exps <- c(0, 2, 3, 5)
  b <- create.monomial.basis(rangeval = dom, exponents = exps)
  for (d in 0:2) {
    V <- eval.basis(t, b, d)
    add_case(
      name = sprintf("monomial_customexp_eval_d%d", d),
      r_call = sprintf("eval.basis(t, create.monomial.basis(c(0,1), exponents=c(0,2,3,5)), %d)", d),
      input = list(domain = dom, exponents = exps, t = t, deriv = d),
      output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
    )
  }
  for (L in 0:2) {
    P <- eval.penalty(b, int2Lfd(L))
    add_case(
      name = sprintf("monomial_customexp_penalty_L%d", L),
      r_call = sprintf(
        "eval.penalty(create.monomial.basis(c(0,1), exponents=c(0,2,3,5)), int2Lfd(%d))", L
      ),
      input = list(domain = dom, exponents = exps, lfd = L),
      output = list(penalty = unname(P), nbasis = b$nbasis)
    )
  }
}

exponential_cases <- function() {
  configs <- list(
    list(dom = c(0, 1), ratevec = c(0, 0.5, -1.0, 2)),
    list(dom = c(0, 5), ratevec = c(0, -1, -5)),
    list(dom = c(0, 1), ratevec = c(1, 2, 3))
  )
  for (i in seq_along(configs)) {
    cfg <- configs[[i]]
    dom <- cfg$dom
    rv <- cfg$ratevec
    t <- eval_points(dom)
    b <- create.exponential.basis(rangeval = dom, nbasis = length(rv), ratevec = rv)
    for (d in 0:2) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("exponential_cfg%d_eval_d%d", i, d),
        r_call = sprintf(
          "eval.basis(t, create.exponential.basis(c(%g,%g), %d, c(%s)), %d)",
          dom[1], dom[2], length(rv), paste(rv, collapse = ","), d
        ),
        input = list(domain = dom, ratevec = rv, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in 0:2) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("exponential_cfg%d_penalty_L%d", i, L),
        r_call = sprintf(
          "eval.penalty(create.exponential.basis(c(%g,%g), %d, c(%s)), int2Lfd(%d))",
          dom[1], dom[2], length(rv), paste(rv, collapse = ","), L
        ),
        input = list(domain = dom, ratevec = rv, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }
}

power_cases <- function() {
  configs <- list(
    list(dom = c(0.1, 2), exponents = c(0, 0.5, 1, 2)),
    list(dom = c(0.5, 3), exponents = c(-1, 0, 1))
  )
  for (i in seq_along(configs)) {
    cfg <- configs[[i]]
    dom <- cfg$dom
    exps <- cfg$exponents
    t <- eval_points(dom)
    b <- create.power.basis(rangeval = dom, nbasis = length(exps), exponents = exps)
    for (d in 0:1) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("power_cfg%d_eval_d%d", i, d),
        r_call = sprintf(
          "eval.basis(t, create.power.basis(c(%g,%g), %d, c(%s)), %d)",
          dom[1], dom[2], length(exps), paste(exps, collapse = ","), d
        ),
        input = list(domain = dom, exponents = exps, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in 0:1) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("power_cfg%d_penalty_L%d", i, L),
        r_call = sprintf(
          "eval.penalty(create.power.basis(c(%g,%g), %d, c(%s)), int2Lfd(%d))",
          dom[1], dom[2], length(exps), paste(exps, collapse = ","), L
        ),
        input = list(domain = dom, exponents = exps, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }
}

constant_cases <- function() {
  for (dom in list(c(0, 1), c(0, 365))) {
    b <- create.constant.basis(rangeval = dom)
    t <- eval_points(dom)
    V <- eval.basis(t, b, 0)
    add_case(
      name = sprintf("constant_eval_dom%g_%g_d0", dom[1], dom[2]),
      r_call = sprintf("eval.basis(t, create.constant.basis(c(%g,%g)), 0)", dom[1], dom[2]),
      input = list(domain = dom, t = t, deriv = 0),
      output = list(values = unname(V), nbasis = b$nbasis)
    )
    P <- eval.penalty(b, int2Lfd(0))
    add_case(
      name = sprintf("constant_penalty_dom%g_%g_L0", dom[1], dom[2]),
      r_call = sprintf(
        "eval.penalty(create.constant.basis(c(%g,%g)), int2Lfd(0))", dom[1], dom[2]
      ),
      input = list(domain = dom, lfd = 0),
      output = list(penalty = unname(P), nbasis = b$nbasis)
    )
  }
}

polygonal_cases <- function() {
  argvals_list <- list(
    seq(0, 1, length.out = 5),
    seq(0, 1, length.out = 11),
    seq(0, 1, length.out = 20)
  )
  for (i in seq_along(argvals_list)) {
    av <- argvals_list[[i]]
    b <- create.polygonal.basis(argvals = av)
    t <- eval_points(range(av))
    for (d in 0:1) {
      V <- eval.basis(t, b, d)
      add_case(
        name = sprintf("polygonal_n%d_eval_d%d", length(av), d),
        r_call = sprintf(
          "eval.basis(t, create.polygonal.basis(argvals=c(%s)), %d)",
          paste(av, collapse = ","), d
        ),
        input = list(argvals = av, t = t, deriv = d),
        output = list(values = unname(V), nbasis = b$nbasis, params = as.numeric(b$params))
      )
    }
    for (L in 0:1) {
      P <- eval.penalty(b, int2Lfd(L))
      add_case(
        name = sprintf("polygonal_n%d_penalty_L%d", length(av), L),
        r_call = sprintf(
          "eval.penalty(create.polygonal.basis(argvals=c(%s)), int2Lfd(%d))",
          paste(av, collapse = ","), L
        ),
        input = list(argvals = av, lfd = L),
        output = list(penalty = unname(P), nbasis = b$nbasis)
      )
    }
  }
}

bspline_cases()
bspline_extra_domains()
bspline_explicit_breaks()
fourier_cases()
monomial_cases()
exponential_cases()
power_cases()
constant_cases()
polygonal_cases()

finalize("basis", rtol = 1e-8)
