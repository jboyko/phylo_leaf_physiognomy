# Integration regression test for occurrence-specific fossil LMA placement.

project_root <- normalizePath(".", mustWork = TRUE)
required <- c(
  "README.md",
  "code/04b_lma_fossil_predictions.R",
  "data/fossil_traits.csv",
  "models/lma_pip_components.rds",
  "models/lma_nophy_models.rds"
)

if (!all(file.exists(file.path(project_root, required)))) {
  cat("Skipping LMA fossil occurrence integration test: generated models unavailable.\n")
  quit(status = 0L)
}

test_root <- tempfile("lma-fossil-occurrence-placement-")
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
  "code/04b_lma_fossil_predictions.R",
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(output, "status")
if (is.null(status)) status <- 0L
if (status != 0L) {
  stop("04b_lma_fossil_predictions.R failed:\n", paste(output, collapse = "\n"))
}

input <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)
eligible <- input[is.finite(log10(input$pw2.a.ratio)), ]
placement <- read.csv(
  "tables/lma_fossil_placement_log.csv",
  stringsAsFactors = FALSE
)
predictions <- read.csv(
  "tables/lma_fossil_species_predictions.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(placement) == nrow(eligible),
  all(placement$placed),
  setequal(placement$fossil_name, eligible$fossil_name),
  nrow(predictions) == nrow(eligible),
  setequal(predictions$fossil_name, eligible$fossil_name),
  all(placement$age_ma == eligible$age_ma[match(
    placement$fossil_name, eligible$fossil_name
  )])
)

cat("LMA fossil occurrence placement integration checks passed.\n")

