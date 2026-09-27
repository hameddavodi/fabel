# tests/fixtures/regression_scalar_design.R
#
# One-off R reference for tests/parity/test_regression.py, written to
# tests/fixtures/regression_scalar_design.json. Run from the repository root:
#
#     Rscript tests/fixtures/regression_scalar_design.R
#
# Why this exists: every case in tests/golden/regression.json puts a
# functional term into the design, and R's fRegress evaluates the integrals
# of such a term by an approximate quadrature (measured 1e-8 .. 3e-3 relative
# error, see tests/parity/test_regression.py). None of those cases can
# therefore check the rest of the regression pipeline -- the penalised solve,
# df, GCV, OCV, the coefficient covariance, the standard errors,
# leave-one-out cross-validation and prediction -- at the 1e-8 parity
# tolerance. A design made only of scalar covariates has no integral to
# approximate, so R computes it exactly and every one of those quantities can
# be compared at 1e-8.
#
# Clean-room note: only the public fda API is called (create.constant.basis,
# fdPar, fd, fRegress, fRegress.CV, predict). No fda source is read or copied.

suppressMessages(library(fda))
suppressMessages(library(jsonlite))
stopifnot(packageVersion("fda") == "6.3.0")

set.seed(20260927)
n <- 30
z1 <- rnorm(n)
z2 <- runif(n)
y <- 1.5 - 0.7 * z1 + 2.0 * z2 + rnorm(n, sd = 0.3)
cbasis <- create.constant.basis(c(0, 1))
xfdlist <- list(rep(1, n), z1, z2)
betalist <- list(fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0), fdPar(cbasis, 0, 0))

fr <- fRegress(y, xfdlist, betalist)
yhat <- as.numeric(fr$yhatfdobj)
sigma2 <- sum((y - yhat)^2) / (n - fr$df)
fr_se <- fRegress(y, xfdlist, betalist, y2cMap = diag(n), SigmaE = diag(rep(sigma2, n)))

# fRegress.CV needs a constant covariate as an fd object (see PROGRESS.md).
cv <- fRegress.CV(y, list(fd(matrix(1, 1, n), cbasis), z1, z2), betalist)

z1_new <- c(-1, 0, 0.5, 2)
z2_new <- c(0.1, 0.9, 0.4, 0.0)
# predict.fRegress wants every newdata entry as an fd object.
as_fd <- function(v) fd(matrix(v, 1, length(v)), cbasis)
pred <- predict(fr, newdata = list(as_fd(rep(1, 4)), as_fd(z1_new), as_fd(z2_new)))

out <- list(
  meta = list(
    r_version = R.version.string,
    fda_version = as.character(packageVersion("fda")),
    seed = 20260927,
    rtol = 1e-8
  ),
  input = list(y = y, z1 = z1, z2 = z2, z1_new = z1_new, z2_new = z2_new),
  output = list(
    beta = sapply(fr$betaestlist, function(b) as.numeric(b$fd$coefs)),
    yhat = yhat,
    df = fr$df,
    gcv = fr$GCV,
    ocv = fr$OCV,
    Cmat = unname(fr$Cmat),
    Dmat = as.numeric(fr$Dmat),
    sigma2 = sigma2,
    bvar = unname(fr_se$bvar),
    betastderr = sapply(fr_se$betastderrlist, function(b) as.numeric(b$coefs)),
    sse_cv = as.numeric(cv$SSE.CV),
    errfd_cv = as.numeric(cv$errfd.cv),
    predicted = as.numeric(pred)
  )
)
writeLines(
  jsonlite::toJSON(out, auto_unbox = TRUE, digits = 17, matrix = "rowmajor", pretty = FALSE),
  "tests/fixtures/regression_scalar_design.json"
)
