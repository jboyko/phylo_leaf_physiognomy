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

pipeline <- new.env(parent = globalenv())
invisible(capture.output(sys.source("code/04_fossil_predictions.R", envir = pipeline)))

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

# Removed comparators must not survive in tables or placement logs.
stopifnot(all(placement$tip_scope == "site_occurrence"))
comparison <- pipeline$scenario_results$formal_only$site
stopifnot(!any(c("lm_sp", "pip_sp") %in% names(comparison)))
stopifnot(!any(grepl("^(lm_sp|pip_sp)_", names(read.csv(
  "tables/fossil_taxonomy_sensitivity.csv"
)))))

# Change one site's measured traits for a species shared with another site.
# Every other site's prediction must remain unchanged.
shared <- names(which(vapply(split(input$site, input$species),
  function(x) length(unique(x)) > 1L, logical(1))))
stopifnot(length(shared) > 0L)
changed_site <- input$site[match(shared[1], input$species)]
rows <- pipeline$foss_base$site == changed_site
pipeline$foss_base$ln.leaf.area.mm2[rows] <-
  pipeline$foss_base$ln.leaf.area.mm2[rows] + 2
invisible(capture.output(changed <- pipeline$run_fossil_scenario("formal_only")))
other <- predictions$site != changed_site
stopifnot(isTRUE(all.equal(
  predictions[other, c("mat_pip_site", "map_pip_site")],
  changed$species[other, c("mat_pip_site", "map_pip_site")],
  tolerance = 1e-10, check.attributes = FALSE
)))
other_sites <- comparison$site != changed_site
stopifnot(isTRUE(all.equal(comparison[other_sites, ],
                          changed$site[other_sites, ], tolerance = 1e-10)))
stopifnot(any(abs(predictions$mat_pip_site[!other] -
                 changed$species$mat_pip_site[!other]) > 1e-8))

cat("Fossil occurrence placement and site-local prediction checks passed.\n")
