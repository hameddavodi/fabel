# tools/export_dynamics_inputs.R
#
# Writes tests/parity/data/dynamics_inputs.json: the smoothed functional-data
# INPUTS of the three pda.fd golden cases in tools/golden_r/dynamics.R.
#
# tests/golden/dynamics.json records the design (raw data, bases) and pda.fd's
# OUTPUTS, but not the smoothed curves pda.fd was run on.  The lip fit
# (order-6 B-spline, D^4 penalty, lambda = 1e-8) has a normal matrix with
# condition number ~3.7e9, so two correct double-precision solvers disagree in
# the 7th significant digit of its coefficients.  Parity of pda.fd itself at
# rtol 1e-8 therefore needs R's own curves as input.
#
# Clean-room note: only the public fda API is called (create.bspline.basis,
# fdPar, smooth.basis); no fda source is read.
#
# Usage: Rscript tools/export_dynamics_inputs.R   (from the repository root)

suppressMessages(library(fda))
suppressMessages(library(jsonlite))
stopifnot(packageVersion("fda") == "6.3.0")

tvec <- seq(0, 1, len = 101)
xbasis <- create.bspline.basis(c(0, 1), 24, 5)
expfd <- smooth.basis(tvec, exp(-4 * tvec), xbasis)$fd

liprange <- range(liptime)
lipbasis <- create.bspline.basis(liprange, nbasis = 31, norder = 6)
lipfd <- smooth.basis(liptime, lip, fdPar(lipbasis, 4, 1e-8))$fd

out <- list(
  `_provenance` = paste(
    "Coefficients of the smoothed curves passed to pda.fd by",
    "tools/golden_r/dynamics.R, written by tools/export_dynamics_inputs.R",
    "with", paste0("R ", R.version$major, ".", R.version$minor), "and fda",
    as.character(packageVersion("fda")), "via the public API only",
    "(smooth.basis); jsonlite digits = 17, row-major."
  ),
  exp_decay_fd_coefs = unname(expfd$coefs),
  lip_fd_coefs = unname(lipfd$coefs)
)
writeLines(
  jsonlite::toJSON(out, digits = 17, matrix = "rowmajor", auto_unbox = TRUE, pretty = TRUE),
  "tests/parity/data/dynamics_inputs.json"
)
