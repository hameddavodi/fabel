# tools/export_datasets.R
#
# Dumps each real fda dataset referenced by tools/golden_r/datasets.R to
# data_export/<name>.json in full precision, for consumption by the
# datasets-loader agent/tests. `data_export/` is gitignored -- these dumps
# are regenerable from fda itself and are not committed.
#
# Usage:
#   Rscript --vanilla tools/export_datasets.R
#
# Output format: jsonlite::toJSON(x, auto_unbox=TRUE, digits=17,
# matrix="rowmajor", na="null", null="null", pretty=FALSE) applied directly
# to the dataset object -- a list-of-components object (e.g.
# CanadianWeather, growth) serializes as a JSON object of its named
# components; a bare array/vector (e.g. handwrit, ReginaPrecip) serializes
# as nested/flat arrays; a data.frame (refinery, seabird) serializes as a
# JSON object of columns (jsonlite's default data.frame representation),
# not an array of row objects.
#
# Clean-room note: only the public fda API (data loading) and jsonlite are
# used. No fda source is read or copied.

suppressMessages(library(fda))
suppressMessages(library(jsonlite))

stopifnot(packageVersion("fda") == "6.3.0")

OUT_DIR <- file.path("data_export")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dataset_names <- c(
  "CanadianWeather", "growth", "gait", "handwrit", "handwritTime",
  "pinch", "melanoma", "refinery", "seabird", "ReginaPrecip",
  "MontrealTemp", "daily", "infantGrowth", "nondurables",
  "lip", "liptime", "lipmarks", "pinchtime", "pinchraw"
)

for (nm in dataset_names) {
  x <- get(nm)
  json <- jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    digits = 17,
    matrix = "rowmajor",
    na = "null",
    null = "null",
    pretty = FALSE
  )
  out_path <- file.path(OUT_DIR, paste0(nm, ".json"))
  writeLines(json, out_path)
  cat(sprintf("wrote %s\n", out_path))
}
