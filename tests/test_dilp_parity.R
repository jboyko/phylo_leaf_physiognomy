# From the analysis root, after script 04 has generated the scenario outputs:
# R_LIBS_USER=/path/to/package-library Rscript tests/test_dilp_parity.R
# DILP_SOURCE=/path/to/dilp can also be used during package development.
package_source <- Sys.getenv("DILP_SOURCE")
if (nzchar(package_source)) {
  pkgload::load_all(package_source, quiet = TRUE)
} else {
  library(dilp)
}
source("code/fossil_taxonomy.R")
source("release/code/dilp_preprocessing.R")

align_rows <- function(actual, expected, columns, label) {
  key <- function(data) do.call(paste, c(data[columns], sep = "\r"))
  actual_key <- key(actual)
  expected_key <- key(expected)
  if (anyDuplicated(actual_key) || anyDuplicated(expected_key) ||
      !setequal(actual_key, expected_key)) {
    stop(label, " identifiers differ.")
  }
  actual[match(expected_key, actual_key), , drop = FALSE]
}

compare_values <- function(actual, expected, label, tolerance = 0) {
  if (length(actual) != length(expected) || !length(expected)) {
    stop(label, " lengths differ or reference is empty.")
  }
  difference <- max(abs(actual - expected))
  if (!is.finite(difference) || difference > tolerance) {
    stop(label, " maximum absolute difference ", difference,
         " exceeds ", tolerance, ".")
  }
  difference
}

verify_lma <- function(raw, figure, submission) {
  actual <- lma_pgls(raw)
  expected <- figure[figure$target == "LMA" & figure$record == "fossil_occurrence", ]
  names(expected)[names(expected) == "id"] <- "fossil_name"
  tolerance <- 1e-7 # Absolute g/m^2
  occurrences <- align_rows(actual$species_predictions, expected,
                           c("site", "fossil_name"), "LMA occurrences")
  compare_values(occurrences$age_ma, expected$age_ma, "LMA occurrence ages")
  compare_values(log10(occurrences$petiole_metric), expected$x,
                 "Submitted log10 petiole metric", 1e-12)
  occurrence_diff <- compare_values(occurrences$LMA.PIP, expected$y,
                                    "Submitted LMA occurrences", tolerance)
  placement <- read.csv("tables/lma_fossil_placement_log.csv")
  placed <- align_rows(actual$placement_log, placement, c("site", "fossil_name"),
                       "LMA placements")
  stopifnot(all(placed$placed),
            identical(placed$placement_level, placement$placement_level),
            identical(placed$placement_target, placement$placement_target),
            identical(placed$age_fallback, placement$age_fallback))

  reference <- read.csv(file.path(submission, "tables/Table_7.csv"))
  sites <- align_rows(actual$results, reference, "site", "LMA sites")
  compare_values(sites$age_ma, reference$age_ma, "LMA site ages")
  compare_values(sites$n_occurrences, reference$n_occurrences, "LMA site counts")
  eligible_sites <- reference$n_occurrences > 0
  stopifnot(identical(is.na(sites$LMA.PIP), !eligible_sites),
            identical(actual$missing_sites$site, "Cerrejon"))
  site_diff <- compare_values(sites$LMA.PIP[eligible_sites],
                              reference$lma_phylodilp_g_m2[eligible_sites],
                              "Submitted LMA site estimates", tolerance)

  climate <- figure[figure$target == "MAT" & figure$record == "fossil_occurrence", ]
  expected_missing <- climate[!climate$id %in% expected$fossil_name, ]
  names(expected_missing)[names(expected_missing) == "id"] <- "fossil_name"
  missing <- align_rows(actual$missing_occurrences, expected_missing,
                        c("site", "fossil_name"), "LMA missing occurrences")
  compare_values(missing$age_ma, expected_missing$age_ma, "Missing occurrence ages")
  compare_values(sites$n_missing_petiole,
                 vapply(sites$site, function(s) sum(missing$site == s), integer(1)),
                 "LMA missing site counts")
  cat("LMA matched", nrow(expected), "occurrences and", sum(eligible_sites),
      "site estimates. Maximum absolute differences (occurrence, site):",
      format(c(occurrence_diff, site_diff), digits = 3), "g/m^2\n")
  cat(nrow(missing), "occurrences lacked usable petiole measurements.",
      "No site estimate for", actual$missing_sites$site, "\n")
}

verify_package <- function() {
  submission <- "submission/submission_1"
  scratch <- tempfile("dilp-parity-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  leaves <- "data/Peppe_2011_fossil_data_June_2026_leaf_level_clean.csv"
  site_ages <- "data/fossil_site_ages.csv"
  unzip(file.path(submission, "Data_S1.zip"), c(leaves, site_ages), exdir = scratch)
  raw <- read_mixed_utf8_csv(file.path(scratch, leaves))
  ages <- read.csv(file.path(scratch, site_ages))
  # Apply the study's site and morphotype labels before passing in raw traits.
  raw$site <- normalise_fossil_site(raw$site)
  raw$morphotype <- trimws(raw$morphotype)
  raw$age_ma <- ages$age_ma[match(raw$site, ages$site)]
  stopifnot(!anyNA(raw$age_ma))
  standard <- suppressWarnings(dilp(raw))
  bundled <- pip_preprocess_dilp(raw)
  stopifnot(identical(bundled, standard[names(bundled)]))

  figure <- read.csv(file.path(submission, "figures/Figure_2_data.csv"))
  table <- read.csv(file.path(submission, "tables/Table_6.csv"))
  models <- c(mat = "MAT.PIP", log_map = "MAP.PIP")
  reference_columns <- c(mat = "mat_pip_site", log_map = "map_pip_site")
  figure_targets <- c(mat = "MAT", log_map = "MAP")
  table_targets <- c(mat = "mat", log_map = "map")
  # Absolute tolerances in response units, across platforms.
  tolerances <- c(mat = 1e-7, log_map = 1e-5)

  cat("Testing DiLP", as.character(packageVersion("dilp")),
      "from", getNamespaceInfo(asNamespace("dilp"), "path"), "\n")
  for (scenario in c("formal_only", "include_informal")) {
    actual <- dilp_pgls(raw, taxonomy_scenario = scenario)
    stopifnot(identical(actual$processed_leaf_data, standard$processed_leaf_data),
              identical(actual$processed_site_data,
                        as.data.frame(standard$processed_site_data)))
    expected <- read.csv(paste0("tables/fossil_predictions_", scenario, ".csv"))
    expected_sites <- read.csv(paste0("tables/fossil_site_uncertainty_", scenario, ".csv"))
    expected_placement <- read.csv(paste0("tables/fossil_placement_log_", scenario, ".csv"))
    stopifnot(nrow(expected) == 360L, nrow(actual$results) == 10L,
              all(actual$placement_log$placed))
    occurrences <- align_rows(actual$species_predictions, expected,
                             c("fossil_name", "site", "species"), scenario)
    compare_values(occurrences$age_ma, expected$age_ma, "Occurrence ages")
    placement <- align_rows(actual$placement_log, expected_placement,
                            c("fossil_name", "site", "species"), scenario)
    stopifnot(identical(placement$placement_level, expected_placement$placement_level),
              identical(placement$placement_target, expected_placement$placement_target),
              identical(placement$age_fallback, expected_placement$age_fallback),
              setequal(actual$unmatched_taxa$fossil_name,
                       placement$fossil_name[placement$placement_level == "root"]))
    sites <- align_rows(actual$results, table, "site", scenario)
    compare_values(sites$n_species, table$mat_n_occurrences, "Site counts")

    for (target in names(models)) {
      label <- paste(scenario, figure_targets[[target]])
      column <- models[[target]]
      lower <- paste0(column, ".lower_95")
      upper <- paste0(column, ".upper_95")
      tolerance <- tolerances[[target]]
      reference <- expected_sites[expected_sites$target == target, ]
      by_site <- align_rows(sites, reference, "site", label)
      compare_values(by_site$n_species, reference$n_species, "Scenario site counts")
      diffs <- c(
        occurrence = compare_values(occurrences[[column]],
                                    expected[[reference_columns[[target]]]],
                                    paste(label, "occurrences"), tolerance),
        site = compare_values(by_site[[column]], reference$estimate_response,
                              paste(label, "sites"), tolerance),
        lower = compare_values(by_site[[lower]], reference$calibrated_lower_response,
                               paste(label, "lower endpoints"), tolerance),
        upper = compare_values(by_site[[upper]], reference$calibrated_upper_response,
                               paste(label, "upper endpoints"), tolerance)
      )

      if (scenario == "formal_only") {
        from_figure <- figure[figure$target == figure_targets[[target]] &
                                figure$record == "fossil_occurrence", ]
        names(from_figure)[names(from_figure) == "id"] <- "fossil_name"
        plotted <- align_rows(occurrences, from_figure,
                             c("fossil_name", "site"), label)
        compare_values(plotted$age_ma, from_figure$age_ma, "Submitted occurrence ages")
        compare_values(plotted[[column]], from_figure$y,
                       paste(label, "submitted figure"), tolerance)
        figure_sites <- figure[figure$target == figure_targets[[target]] &
                                  figure$record == "fossil_site", ]
        plotted_sites <- align_rows(sites, figure_sites, "site", label)
        compare_values(plotted_sites$n_species, figure_sites$n_occurrences,
                       "Submitted figure site counts")
        compare_values(plotted_sites[[column]], figure_sites$y,
                       paste(label, "submitted figure sites"), tolerance)
        compare_values(plotted_sites[[lower]], figure_sites$lower_95,
                       paste(label, "submitted figure lower endpoints"), tolerance)
        compare_values(plotted_sites[[upper]], figure_sites$upper_95,
                       paste(label, "submitted figure upper endpoints"), tolerance)
        prefix <- table_targets[[target]]
        compare_values(sites[[column]], table[[paste0(prefix, "_phylodilp")]],
                       paste(label, "submitted site estimates"), tolerance)
        compare_values(sites[[lower]], table[[paste0(prefix, "_lower_95")]],
                       paste(label, "submitted lower endpoints"), tolerance)
        compare_values(sites[[upper]], table[[paste0(prefix, "_upper_95")]],
                       paste(label, "submitted upper endpoints"), tolerance)
      }

      cat(label, "matched 360 occurrences and 10 sites.",
          "Maximum absolute differences (occurrence, site, lower, upper):",
          format(diffs, digits = 3), "\n")
    }
  }
  verify_lma(raw, figure, submission)
}

verify_package()
