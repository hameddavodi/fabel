# tools/golden_r/datasets.R
#
# Golden cases for the `datasets` module: NOT numerical parity cases like
# the other modules -- this records structure + checksums (spot values and
# sums) for each real fda dataset, so the Python-side dataset loaders can
# verify they parsed the exported data (see tools/export_datasets.R)
# correctly, without committing the datasets themselves to git.
#
# Driven by tools/make_golden.py, which sources common.R first (defines
# add_case, finalize, GOLDEN_OUT, GOLDEN_SEED).
#
# Clean-room note: only the public fda API / base R introspection is used
# (class, dim, names, head, tail, sum on the dataset objects themselves).
# No fda source is read or copied.
#
# fda 6.3.0 dataset catalog note: the task spec asked for a "CSTR" dataset
# and a separate "Chinese script" dataset. Neither exists in fda 6.3.0 --
# the full catalog (data(package="fda")$results[, "Item"]) is recorded in
# meta.catalog below. The handwriting-related dataset that does exist is
# `handwrit` (+ its companion time vector `handwritTime`), which is covered
# here instead.

#' First/last 3 values (in storage order) and the sum of a numeric object.
numeric_spot <- function(x) {
  v <- as.numeric(x)
  list(
    head3 = utils::head(v, 3),
    tail3 = utils::tail(v, 3),
    sum = sum(v, na.rm = TRUE)
  )
}

#' Compact, JSON-friendly description of one dataset component (a column,
#' a list element, or the object itself for a bare array/vector).
component_info <- function(x) {
  info <- list(class = class(x)[1])
  if (is.numeric(x)) {
    info$dim <- if (!is.null(dim(x))) as.integer(dim(x)) else length(x)
    info <- c(info, numeric_spot(x))
  } else if (is.character(x) || is.factor(x)) {
    v <- as.character(x)
    info$length <- length(v)
    info$head3 <- utils::head(v, 3)
    info$tail3 <- utils::tail(v, 3)
  } else {
    info$length <- length(x)
  }
  info
}

#' Structure + spot-value summary for one named dataset in the fda package.
dataset_info <- function(name) {
  x <- get(name)
  if (is.data.frame(x)) {
    components <- lapply(as.list(x), component_info)
    names(components) <- names(x)
    list(class = "data.frame", dim = as.integer(dim(x)), columns = names(x), components = components)
  } else if (is.list(x) && is.null(dim(x))) {
    components <- lapply(x, component_info)
    names(components) <- names(x)
    list(class = "list", names = names(x), components = components)
  } else {
    list(class = class(x)[1], components = list(value = component_info(x)))
  }
}

dataset_names <- c(
  "CanadianWeather", "growth", "gait", "handwrit", "handwritTime",
  "pinch", "melanoma", "refinery", "seabird", "ReginaPrecip"
)

for (nm in dataset_names) {
  add_case(
    name = paste0("dataset_", nm),
    r_call = sprintf("str(%s)  # structure/spot-values only, see component_info()", nm),
    input = list(dataset = nm),
    output = dataset_info(nm)
  )
}

catalog <- data(package = "fda")$results[, "Item"]

meta_note <- list(
  catalog = as.character(catalog),
  note = paste(
    "No 'CSTR' dataset and no separate 'Chinese script' dataset exist in",
    "fda 6.3.0 (see catalog above). The handwriting dataset present is",
    "`handwrit` (+ `handwritTime`), covered by the dataset_handwrit and",
    "dataset_handwritTime cases."
  )
)
add_case(
  name = "dataset_catalog",
  r_call = "data(package='fda')$results[, 'Item']",
  input = list(query = "data(package='fda')$results[, 'Item']"),
  output = meta_note
)

finalize("datasets", rtol = 1e-8)
