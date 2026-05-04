devtools::load_all("~/dilp")
library(tidyverse)
library(ape)
library(phytools)

setwd("/Users/jboyko/phylo_leaf_physiognomy")

BUILD_PHYLOGENY <- FALSE

input_1 <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv", fileEncoding = "latin1") #or write in your file path
input_2 <- read.csv("data/extra_calibration_data_for_LMA.csv", fileEncoding = "latin1")

# Combine leaf-level and petiole-width-only datasets (same columns; input_2 is already species-site level)
lma_cols <- c("specimen_number", "site", "morphotype", "order", "family", "genus", "species",
              "margin", "lma_g_m2",
              "petiole_width", "petiole_area", "blade_area", "blade_perimeter",
              "feret", "minimum_feret", "raw_blade_area", "raw_blade_perimeter",
              "internal_raw_blade_area", "internal_raw_blade_perimeter",
              "length_of_cut_perimeter", "no_of_primary_teeth", "no_of_subsidiary_teeth",
              "MAT_C", "MAP_mm", "latitude_N", "longitude_E")
input <- bind_rows(
  input_1 %>% select(any_of(lma_cols)),
  input_2 %>% select(any_of(lma_cols))
)

input <- dilp(input)  #runs the dilp package
input.leaf <- input$processed_leaf_data #includes calculated variables at the leaf level
input.morphotype <- input$processed_morphotype_data #includes calculated variables at the species-site level

# ==============================================================================
# 2. ADD TAXONOMY AND FILL TOOTH TRAITS
# ==============================================================================

# Join taxonomy back via site+morphotype (safer than positional indexing)
taxonomy_lookup <- input.leaf %>%
  select(site, morphotype, order, family, genus, species) %>%
  distinct()

input.morphotype <- input.morphotype %>%
  select(-any_of("measurer_comments")) %>%
  left_join(taxonomy_lookup, by = c("site", "morphotype")) %>%
  relocate(order, family, genus, species, .after = morphotype)

# Fill tooth traits for confirmed untoothed species (margin == 1)
input.morphotype <- input.morphotype %>%
  mutate(
    tc_p            = ifelse(margin == 1, 0, tc_p),
    tc_ip           = ifelse(margin == 1, 0, tc_ip),
    tc_ba           = ifelse(margin == 1, 0, tc_ba),
    ta_p            = ifelse(margin == 1, 0, ta_p),
    ta_ip           = ifelse(margin == 1, 0, ta_ip),
    ta_ba           = ifelse(margin == 1, 0, ta_ba),
    avg_ta          = ifelse(margin == 1, 0, avg_ta),
    perimeter_ratio = ifelse(margin == 1, 1, perimeter_ratio)
  )

# ==============================================================================
# 3. AGGREGATE TO SPECIES LEVEL
# ==============================================================================

# Create genus_species key; drop rows with no genus/species identification
input.morphotype <- input.morphotype %>%
  filter(!is.na(genus), genus != "", !is.na(species), species != "") %>%
  mutate(genusSpecies = paste(genus, species, sep = "_"))

numeric_cols <- names(input.morphotype)[sapply(input.morphotype, is.numeric)]

dat_sp_lma <- input.morphotype %>%
  group_by(genusSpecies, genus, family, order) %>%
  summarise(across(all_of(numeric_cols), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

# Log-transform LMA and petiole_metric for modelling
dat_sp_lma <- dat_sp_lma %>%
  mutate(
    log10_lma           = log10(lma_g_m2),
    log10_petiole_metric = log10(petiole_metric)
  )

# ==============================================================================
# 4. MATCH TO PHYLOGENY
# ==============================================================================

if (BUILD_PHYLOGENY) {

# Reuse the scaffold tree built in 00_data_cleaning.R (family backbone +
# all climate training species). Graft any LMA species not already present.

tree_scaffold <- read.tree("data/tre_scaffold.tre")
name_table    <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)

# Only use backbone tips (those with valid name_table entries) for MRCA lookups.
# Grafted training species (905 tips) have no name_table rows and would produce
# NAs in genus/family/order comparisons, corrupting the tip indexing.
all_labels   <- tree_scaffold$tip.label
bb_matches   <- match(all_labels, paste(name_table$genus, name_table$species, sep = "_"))
bb_valid     <- !is.na(bb_matches)
orig_labels_bb      <- all_labels[bb_valid]
name_table_backbone <- name_table[bb_matches[bb_valid], ]

h <- max(nodeHeights(tree_scaffold))

cat("Grafting LMA species onto scaffold...\n")
for (i in seq_len(nrow(dat_sp_lma))) {
  cat("\r  ", round(i / nrow(dat_sp_lma) * 100, 1), "%")
  current_sp <- dat_sp_lma$genusSpecies[i]
  if (current_sp %in% tree_scaffold$tip.label) next

  target_node <- NULL

  g_tips <- orig_labels_bb[name_table_backbone$genus == dat_sp_lma$genus[i]]
  g_tips <- g_tips[g_tips %in% tree_scaffold$tip.label]
  if (length(g_tips) >= 2) {
    target_node <- getMRCA(tree_scaffold, g_tips)
  } else if (length(g_tips) == 1) {
    target_node <- which(tree_scaffold$tip.label == g_tips[1])
  }

  if (is.null(target_node)) {
    f_tips <- orig_labels_bb[name_table_backbone$family == dat_sp_lma$family[i]]
    f_tips <- f_tips[f_tips %in% tree_scaffold$tip.label]
    if (length(f_tips) >= 2) {
      target_node <- getMRCA(tree_scaffold, f_tips)
    } else if (length(f_tips) == 1) {
      target_node <- which(tree_scaffold$tip.label == f_tips[1])
    }
  }

  if (is.null(target_node)) {
    o_tips <- orig_labels_bb[name_table_backbone$order == dat_sp_lma$order[i]]
    o_tips <- o_tips[o_tips %in% tree_scaffold$tip.label]
    if (length(o_tips) >= 2) {
      target_node <- getMRCA(tree_scaffold, o_tips)
    } else if (length(o_tips) == 1) {
      target_node <- which(tree_scaffold$tip.label == o_tips[1])
    }
  }

  if (!is.null(target_node)) {
    edge_len      <- h - nodeheight(tree_scaffold, target_node)
    tree_scaffold <- bind.tip(tree_scaffold, tip.label = current_sp,
                              where = target_node, edge.length = edge_len)
  } else {
    warning("Could not place LMA species: ", current_sp)
  }
}
cat("\n")

# Pruned tree: LMA species only
tree_lma_pruned <- keep.tip(tree_scaffold,
                             dat_sp_lma$genusSpecies[dat_sp_lma$genusSpecies %in%
                                                       tree_scaffold$tip.label])

# Align species data to tree tip order
dat_sp_lma <- dat_sp_lma[match(tree_lma_pruned$tip.label, dat_sp_lma$genusSpecies), ]

write.tree(tree_lma_pruned, file = "data/tre_lma_pruned.tre")
cat("Wrote data/tre_lma_pruned.tre (", length(tree_lma_pruned$tip.label), "tips)\n")

} # end BUILD_PHYLOGENY

# ==============================================================================
# 5. WRITE OUTPUT
# ==============================================================================

write.csv(dat_sp_lma, file = "data/lma_species.csv", row.names = FALSE)
cat("Wrote data/lma_species.csv (", nrow(dat_sp_lma), "species)\n")

