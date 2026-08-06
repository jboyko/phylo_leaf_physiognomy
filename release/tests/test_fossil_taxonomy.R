source("code/fossil_taxonomy.R")

fixture <- data.frame(
  genus = c('"Vitis"', "Acer", ""),
  family = c("Platanaceae", '"Salicaceae"', ""),
  order = c("Proteales", "Malpighiales", '"Ericales"'),
  stringsAsFactors = FALSE
)

formal <- taxonomy_for_scenario(fixture, "formal_only")
provisional <- taxonomy_for_scenario(fixture, "include_informal")

stopifnot(
  !is_known_taxon("unknown"),
  !is_known_taxon(""),
  is_known_taxon("Acer"),
  identical(formal$genus, c("unknown", "Acer", "unknown")),
  identical(formal$family, c("Platanaceae", "unknown", "unknown")),
  identical(formal$order, c("Proteales", "Malpighiales", "unknown")),
  identical(provisional$genus, c("Vitis", "Acer", "unknown")),
  identical(provisional$family, c("Platanaceae", "Salicaceae", "unknown")),
  identical(provisional$order, c("Proteales", "Malpighiales", "Ericales"))
)

ages <- fossil_site_ages()
stopifnot(
  identical(unname(ages["Fox Hills"]), 66.5),
  identical(unname(ages["Williston Basin I"]), 64.75),
  identical(unname(ages["Williston Basin II"]), 63.5),
  identical(unname(ages["Williston Basin III"]), 59.75)
)

stopifnot(
  identical(palacio_analysis_site("SA049"), "Palacio de los Loros PL1"),
  identical(palacio_analysis_site("SA050"), "Palacio de los Loros PL2"),
  identical(palacio_analysis_site("SA060"), "Palacio de los Loros PL1")
)

raw <- read_mixed_utf8_csv(
  "data/Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv"
)
stopifnot(
  nrow(raw) == 1413L,
  length(unique(paste(trimws(raw$site), trimws(raw$morphotype)))) == 365L
)

fossil_traits <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)
stopifnot(
  nrow(fossil_traits) == 361L,
  !anyDuplicated(paste(fossil_traits$site, fossil_traits$species)),
  sum(fossil_traits$genus_informal) == 20L,
  sum(fossil_traits$family_informal) == 12L,
  sum(fossil_traits$order_informal) == 6L,
  all(fossil_traits$genus[fossil_traits$genus_informal] == "unknown"),
  all(fossil_traits$family[fossil_traits$family_informal] == "unknown"),
  all(fossil_traits$order[fossil_traits$order_informal] == "unknown"),
  all(fossil_traits$age_ma == unname(ages[fossil_traits$site]))
)

cat("Fossil taxonomy regression checks passed.\n")
