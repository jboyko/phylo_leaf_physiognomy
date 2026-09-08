# Integration regression test for occurrence-specific fossil placement.
#
# This test intentionally exercises the real script because placement was
# historically nested inside 04_fossil_predictions.R. It runs in an isolated
# temporary project root so generated tables never overwrite working outputs.

project_root <- normalizePath(".", mustWork = TRUE)
required <- c(
  "README.md",
  "code/04_fossil_predictions.R",
  "data/fossil_traits.csv",
  "models/pip_components.rds",
  "models/nophy_models.rds",
  "models/site_models.rds"
)

if (!all(file.exists(file.path(project_root, required)))) {
  cat("Skipping fossil occurrence integration test: generated models unavailable.\n")
  quit(status = 0L)
}

test_root <- tempfile("fossil-occurrence-placement-")
dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

for (path in c("README.md", "code", "data", "models")) {
  linked <- file.symlink(file.path(project_root, path), file.path(test_root, path))
  if (!linked) stop("Could not create test symlink for ", path)
}
dir.create(file.path(test_root, "tables"))

old_wd <- setwd(test_root)
on.exit(setwd(old_wd), add = TRUE)

output <- system2(
  file.path(R.home("bin"), "Rscript"),
  "code/04_fossil_predictions.R",
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(output, "status")
if (is.null(status)) status <- 0L
if (status != 0L) {
  stop("04_fossil_predictions.R failed:\n", paste(output, collapse = "\n"))
}

placement <- read.csv(
  "tables/fossil_placement_log_formal_only.csv",
  stringsAsFactors = FALSE
)
occurrence <- placement[placement$tip_scope == "site_occurrence", ]
input <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)

if (nrow(occurrence) != nrow(input)) {
  stop("Expected ", nrow(input), " occurrence placements; found ", nrow(occurrence))
}
if (!all(occurrence$placed)) {
  failed <- occurrence[!occurrence$placed, c("fossil_name", "site", "error")]
  stop(
    "Every fossil occurrence must be placed; failures:\n",
    paste(capture.output(print(failed, row.names = FALSE)), collapse = "\n")
  )
}

predictions <- read.csv(
  "tables/fossil_predictions_formal_only.csv",
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(predictions) == nrow(input),
  setequal(predictions$fossil_name, input$fossil_name)
)

cat("Fossil occurrence placement integration checks passed.\n")
