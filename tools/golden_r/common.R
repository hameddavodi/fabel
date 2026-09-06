# tools/golden_r/common.R
#
# Shared helpers for golden-file generator R scripts. Sourced by
# tools/make_golden.py BEFORE the per-module script (tools/golden_r/<module>.R),
# via both the rpy2 path and the Rscript subprocess path, so module scripts
# can assume these are already in the global environment.
#
# The driver (tools/make_golden.py) sets these globals before sourcing this
# file:
#   GOLDEN_OUT   - absolute path to the output tests/golden/<module>.json
#   GOLDEN_SEED  - integer seed (20260906)
#
# Clean-room note: this file only CALLS the public fda API (create.*.basis,
# eval.basis, getbasispenalty, fd, eval.fd, inprod, ...). No fda source code
# is read or copied here.

suppressMessages(library(fda))
suppressMessages(library(jsonlite))

stopifnot(packageVersion("fda") == "6.3.0")

set.seed(GOLDEN_SEED)

CASES <- list()

#' Append one golden case to the module-level CASES accumulator.
add_case <- function(name, r_call, input, output) {
  CASES[[length(CASES) + 1]] <<- list(
    name = name,
    r_call = r_call,
    input = input,
    output = output
  )
}

#' 41 equally spaced points across rangeval, plus a few interior points,
#' always including both endpoints.
eval_points <- function(rangeval) {
  lo <- rangeval[1]
  hi <- rangeval[2]
  main <- seq(lo, hi, length.out = 41)
  extra <- lo + c(0.1234, 0.37, 0.618, 0.859) * (hi - lo)
  sort(unique(c(main, extra)))
}

#' Compact, JSON-friendly description of a basisfd object.
basis_info <- function(b) {
  info <- list(
    type = b$type,
    nbasis = b$nbasis,
    rangeval = as.numeric(b$rangeval)
  )
  if (!is.null(b$params)) {
    info$params <- as.numeric(b$params)
  }
  info
}

#' Write meta + cases to GOLDEN_OUT as JSON (jsonlite, full precision,
#' row-major matrices, auto-unboxed scalars, NA -> null).
finalize <- function(module_name, rtol) {
  meta <- list(
    module = module_name,
    r_version = R.version.string,
    fda_version = as.character(packageVersion("fda")),
    generated = format(Sys.Date(), "%Y-%m-%d"),
    seed = GOLDEN_SEED,
    rtol = rtol
  )
  obj <- list(meta = meta, cases = CASES)
  json <- jsonlite::toJSON(
    obj,
    auto_unbox = TRUE,
    digits = NA,
    matrix = "rowmajor",
    na = "null",
    null = "null",
    pretty = FALSE
  )
  writeLines(json, GOLDEN_OUT)
  invisible(NULL)
}
