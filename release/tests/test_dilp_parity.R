# Run after 04 using either an installed updated dilp or a source checkout:
# DILP_SOURCE=/path/to/dilp Rscript tests/test_dilp_parity.R
package_source <- Sys.getenv("DILP_SOURCE")
if (nzchar(package_source)) {
  pkgload::load_all(package_source, quiet = TRUE)
} else {
  library(dilp)
}
source("code/fossil_taxonomy.R")

raw <- read_mixed_utf8_csv("data/Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv")
raw$site <- trimws(raw$site)
raw$morphotype <- trimws(raw$morphotype)
is_palacio <- raw$site == "Palacio de los Loros"
raw$site[is_palacio] <- palacio_analysis_site(raw$morphotype[is_palacio])
raw$age_ma <- unname(fossil_site_ages()[raw$site])
# dilp_pgls's species field is a complete operational taxon identifier.
genus_label <- strip_taxon_quotes(raw$genus)
species_epithet <- strip_taxon_quotes(raw$species)
genus_label[is.na(genus_label)] <- ""
species_epithet[is.na(species_epithet)] <- ""
raw$species <- ifelse(nzchar(genus_label) & nzchar(species_epithet),
  trimws(gsub("[[:space:]]+", " ", paste(genus_label, species_epithet))), raw$morphotype)

for (scenario in c("formal_only", "include_informal")) {
  actual <- dilp_pgls(raw, taxonomy_scenario = scenario)
  expected <- read.csv(paste0("tables/fossil_predictions_", scenario, ".csv"))
  actual_sp <- actual$species_predictions
  key <- function(d) paste(d$site, d$species, sep = "\r")
  stopifnot(nrow(actual_sp) == 361L, nrow(expected) == 361L,
            all(actual$placement_log$placed),
            !anyDuplicated(key(actual_sp)), setequal(key(actual_sp), key(expected)))
  actual_sp <- actual_sp[match(key(expected), key(actual_sp)), ]
  stopifnot(isTRUE(all.equal(actual_sp$MAT.PIP, expected$mat_pip_site, tolerance = 1e-8)),
            isTRUE(all.equal(actual_sp$MAP.PIP, expected$map_pip_site, tolerance = 1e-8)),
            identical(actual_sp$age_ma, expected$age_ma),
            all(is.finite(actual_sp$MAT.PIP)), all(is.finite(actual_sp$MAP.PIP)))
  expected_sites <- aggregate(cbind(mat_pip_site, map_pip_site) ~ site,
    transform(expected, map_pip_site = log(map_pip_site)), mean)
  observed_sites <- actual$results[match(expected_sites$site, actual$results$site), ]
  stopifnot(isTRUE(all.equal(observed_sites$MAT.PIP, expected_sites$mat_pip_site, tolerance = 1e-8)),
            isTRUE(all.equal(observed_sites$MAP.PIP, exp(expected_sites$map_pip_site), tolerance = 1e-8)))
  cat("Raw fossil pipeline/package parity passed:", scenario, "361 occurrences\n")
}
