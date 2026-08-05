# ==============================================================================
# 00c_fossil_data_cleaning.R
# Build species-by-site fossil traits from Dana Royer's April 2026 leaf-level
# data while preserving formal vs informal taxonomy for sensitivity analysis.
# ==============================================================================

if (!file.exists("README.md")) {
  stop("Run this script from the phylo_leaf_physiognomy repository root.")
}

source("code/fossil_taxonomy.R")

if (!requireNamespace("dilp", quietly = TRUE)) {
  dilp_source <- path.expand("~/dilp")
  if (!dir.exists(dilp_source)) {
    stop("Install the dilp package or place its source at ~/dilp.")
  }
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("The devtools package is required to load the local dilp source.")
  }
  devtools::load_all(dilp_source, quiet = TRUE)
}
dilp_fn <- get("dilp", envir = asNamespace("dilp"))

mean_or_na <- function(x) {
  value <- mean(x, na.rm = TRUE)
  if (is.nan(value)) NA_real_ else value
}

input_path <- "data/Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv"
input <- read_mixed_utf8_csv(input_path)

input$site <- trimws(input$site)
input$morphotype <- trimws(input$morphotype)

# Restore the PL1/PL2 distinction omitted from the April leaf-level CSV.
is_palacio <- input$site == "Palacio de los Loros"
input$site[is_palacio] <- palacio_analysis_site(input$morphotype[is_palacio])

ages <- fossil_site_ages()
unknown_sites <- setdiff(unique(input$site), names(ages))
if (length(unknown_sites) > 0) {
  stop("No age configured for site(s): ", paste(unknown_sites, collapse = ", "))
}
input$age_ma <- unname(ages[input$site])

# Taxonomy must be constant within each site-morphotype before aggregation.
taxonomy_cols <- c("site", "morphotype", "genus", "family", "order", "species")
taxonomy <- unique(input[, taxonomy_cols])
taxonomy_key <- paste(taxonomy$site, taxonomy$morphotype, sep = "\r")
if (anyDuplicated(taxonomy_key)) {
  stop("Conflicting taxonomy found within at least one site-morphotype.")
}

processed <- dilp_fn(input)$processed_morphotype_data

# `dilp()` returns NaN for tooth variables when every leaf in a morphotype is
# untoothed. Match the extant training-data convention explicitly.
untoothed <- !is.na(processed$margin) & processed$margin == 1
tooth_zero <- c("tc_p", "tc_ip", "avg_ta", "ta_ba", "ta_p", "ta_ip", "tc_ba")
for (column in tooth_zero) {
  processed[[column]][untoothed & !is.finite(processed[[column]])] <- 0
}
processed$perimeter_ratio[
  untoothed & !is.finite(processed$perimeter_ratio)
] <- 1

processed_key <- paste(processed$site, processed$morphotype, sep = "\r")
taxonomy <- taxonomy[match(processed_key, taxonomy_key), ]
if (anyNA(taxonomy$site)) {
  stop("Failed to reattach taxonomy to one or more processed morphotypes.")
}

reported <- data.frame(
  genus_reported = normalise_missing_taxon(taxonomy$genus),
  family_reported = normalise_missing_taxon(taxonomy$family),
  order_reported = normalise_missing_taxon(taxonomy$order),
  stringsAsFactors = FALSE
)
reported$genus_informal <- is_informal_taxon(reported$genus_reported)
reported$family_informal <- is_informal_taxon(reported$family_reported)
reported$order_informal <- is_informal_taxon(reported$order_reported)

genus_label <- strip_taxon_quotes(taxonomy$genus)
species_epithet <- strip_taxon_quotes(taxonomy$species)
genus_label[is.na(genus_label)] <- ""
species_epithet[is.na(species_epithet)] <- ""
named_taxon <- trimws(paste(genus_label, species_epithet))
named_taxon <- gsub("[[:space:]]+", " ", named_taxon)
species_label <- ifelse(
  nzchar(genus_label) & nzchar(species_epithet),
  named_taxon,
  taxonomy$morphotype
)

morphotype_traits <- data.frame(
  species = species_label,
  site = processed$site,
  morphotype = processed$morphotype,
  reported,
  age_ma = unname(ages[processed$site]),
  ln.leaf.area.mm2 = processed$ln_leaf_area,
  feret.diam.ratio = processed$fdr,
  perim.ratio = processed$perimeter_ratio,
  teeth.perimeter.percm = processed$tc_p,
  teeth.interior.percm = processed$tc_ip,
  teeth.blade.area.ratio = processed$tc_ba,
  tooth.area.perimeter = processed$ta_p,
  tooth.area.interior = processed$ta_ip,
  tooth.area.blade.area.ratio = processed$ta_ba,
  avt.tooth.area = processed$avg_ta,
  pw2.a.ratio = processed$petiole_metric,
  margin.score = processed$margin,
  stringsAsFactors = FALSE
)

trait_cols <- c(
  "ln.leaf.area.mm2", "feret.diam.ratio", "perim.ratio",
  "teeth.perimeter.percm", "teeth.interior.percm",
  "teeth.blade.area.ratio", "tooth.area.perimeter",
  "tooth.area.interior", "tooth.area.blade.area.ratio",
  "avt.tooth.area", "pw2.a.ratio", "margin.score"
)

species_site <- aggregate(
  morphotype_traits[, trait_cols],
  by = morphotype_traits[, c("species", "site", "age_ma")],
  FUN = mean_or_na
)

morphotypes <- aggregate(
  morphotype ~ species + site,
  data = morphotype_traits,
  FUN = function(x) paste(sort(unique(x)), collapse = "|")
)

resolve_reported_taxon <- function(x) {
  x <- normalise_missing_taxon(x)
  known <- unique(x[x != "unknown"])
  if (length(known) > 1) {
    stop("Conflicting non-missing taxonomy within a species-site: ",
         paste(known, collapse = ", "))
  }
  if (length(known) == 1) known else "unknown"
}

species_taxonomy <- aggregate(
  cbind(genus_reported, family_reported, order_reported) ~ species + site,
  data = morphotype_traits,
  FUN = resolve_reported_taxon
)
species_taxonomy$genus_informal <- is_informal_taxon(
  species_taxonomy$genus_reported
)
species_taxonomy$family_informal <- is_informal_taxon(
  species_taxonomy$family_reported
)
species_taxonomy$order_informal <- is_informal_taxon(
  species_taxonomy$order_reported
)

species_site <- merge(
  species_site,
  morphotypes,
  by = c("species", "site"),
  all.x = TRUE,
  sort = FALSE
)
species_site <- merge(
  species_site,
  species_taxonomy,
  by = c("species", "site"),
  all.x = TRUE,
  sort = FALSE
)

species_site <- taxonomy_for_scenario(species_site, "formal_only")
species_site$fossil_name <- make.unique(
  gsub("[^[:alnum:]]+", "_", paste(species_site$species, species_site$site))
)

output_cols <- c(
  "fossil_name", "species", "site", "morphotype",
  "genus", "family", "order",
  "genus_reported", "family_reported", "order_reported",
  "genus_informal", "family_informal", "order_informal",
  "age_ma", trait_cols
)
species_site <- species_site[, output_cols]
species_site <- species_site[order(-species_site$age_ma,
                                   species_site$site,
                                   species_site$species), ]
row.names(species_site) <- NULL

if (anyDuplicated(paste(species_site$site, species_site$species))) {
  stop("Output contains duplicate species-site rows.")
}
if (nrow(species_site) != 361L) {
  stop("Expected 361 species-site rows; produced ", nrow(species_site), ".")
}

write.csv(species_site, "data/fossil_traits.csv", row.names = FALSE)

cat("Wrote data/fossil_traits.csv (", nrow(species_site),
    " species-site rows)\n", sep = "")
cat("Informal ranks: genus=", sum(species_site$genus_informal),
    ", family=", sum(species_site$family_informal),
    ", order=", sum(species_site$order_informal), "\n", sep = "")
