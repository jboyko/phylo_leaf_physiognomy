setwd("/Users/jboyko/phylo_leaf_physiognomy")

# Set to TRUE to rebuild phylogenies (slow; only needed when tree or species
# list changes). Set to FALSE to skip sections 4-8 and only regenerate data.
BUILD_PHYLOGENY <- TRUE

library(ape)
library(phytools)
library(dilp)
library(dplyr)

# ==============================================================================
# 1. LOAD AND PROCESS SPECIMEN-LEVEL DATA VIA DILP
# ==============================================================================

raw      <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv",
                     fileEncoding = "latin1")
dilp_out <- dilp(raw)
dat_leaf <- dilp_out$processed_leaf_data
dat      <- dilp_out$processed_morphotype_data  # species×site means

# Re-attach taxonomy columns (dilp drops them) — Dana Royer, pers. comm.
dat <- subset(dat, select = -measurer_comments)

taxonomy_lookup <- dat_leaf %>%
  dplyr::select(site, morphotype, order, family, genus, species) %>%
  distinct()
dat <- dat %>%
  left_join(taxonomy_lookup, by = c("site", "morphotype"))

dat$genusSpecies <- gsub("[() ]", "_", paste(dat$genus, dat$species, sep = "_"))
dat$map          <- dat$map_mm / 10   # mm -> cm (pipeline convention; CLAUDE.md §8)

# ==============================================================================
# 2. FILL TOOTH TRAITS FOR UNTOOTHED LEAVES (Dana Royer, pers. comm.)
# ==============================================================================

# Save pre-fill copy for Peppe-style aggregation (untoothed excluded from tooth
# trait averages via na.rm = TRUE rather than pulled toward zero).
dat_prefill <- dat

# margin == 1 indicates untoothed; tooth trait absence is biologically real
# ONLY when the cell is also blank (NA) — a margin-scored-untoothed leaf that
# still carries a nonzero tooth measurement must not be overwritten with 0.
# Set to 0 (or 1 for perimeter_ratio) before aggregating (Dana Royer, pers. comm.).
dat <- dat %>%
  mutate(
    tc_p            = ifelse(margin == 1 & is.na(tc_p), 0, tc_p),
    tc_ip           = ifelse(margin == 1 & is.na(tc_ip), 0, tc_ip),
    tc_ba           = ifelse(margin == 1 & is.na(tc_ba), 0, tc_ba),
    ta_p            = ifelse(margin == 1 & is.na(ta_p), 0, ta_p),
    ta_ip           = ifelse(margin == 1 & is.na(ta_ip), 0, ta_ip),
    ta_ba           = ifelse(margin == 1 & is.na(ta_ba), 0, ta_ba),
    avg_ta          = ifelse(margin == 1 & is.na(avg_ta), 0, avg_ta),
    perimeter_ratio = ifelse(margin == 1 & is.na(perimeter_ratio), 1, perimeter_ratio)
  )

# Rename dilp output columns to pipeline convention
rename_pipeline <- function(df) {
  df %>% rename(
    pw2.a.ratio                 = petiole_metric,
    ln.leaf.area.mm2            = ln_leaf_area,
    feret.diam.ratio            = fdr,
    margin.score                = margin,
    perim.ratio                 = perimeter_ratio,
    teeth.perimeter.percm       = tc_p,
    teeth.interior.percm        = tc_ip,
    avt.tooth.area              = avg_ta,
    tooth.area.blade.area.ratio = ta_ba,
    tooth.area.perimeter        = ta_p,
    tooth.area.interior         = ta_ip,
    teeth.blade.area.ratio      = tc_ba,
    mat                         = mat_c,
    Site                        = site,
    Family                      = family,
    Order                       = order
  )
}
dat         <- rename_pipeline(dat)
dat_prefill <- rename_pipeline(dat_prefill)

# ==============================================================================
# 3. AGGREGATE TO SPECIES LEVEL
# ==============================================================================

agg_cols <- c("pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
              "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
              "avt.tooth.area", "tooth.area.blade.area.ratio", "tooth.area.perimeter",
              "tooth.area.interior", "teeth.blade.area.ratio",
              "mat", "map", "latitude_n", "longitude_e")

dat_sp <- aggregate(dat[, agg_cols],
  by = list(genusSpecies = dat$genusSpecies,
            genus        = dat$genus,
            Family       = dat$Family,
            Order        = dat$Order),
  FUN = mean, na.rm = TRUE)
dat_sp <- dat_sp[!is.na(dat_sp$Order) & dat_sp$Order != "unknown", ]
dat_sp <- dat_sp[grepl("[A-Za-z]", dat_sp$genusSpecies), ]

# ==============================================================================
# 3b. AGGREGATE TO SITE LEVEL (three variants)
# ==============================================================================

# Variant 1: morphotype → site (tooth-filled; each morphotype equally weighted)
dat_site <- aggregate(dat[, agg_cols],
                      by  = list(Site = dat$Site),
                      FUN = mean, na.rm = TRUE)
write.csv(dat_site, file = "data/dat_site.csv", row.names = FALSE)
cat("Wrote data/dat_site.csv (", nrow(dat_site), "sites)\n")

# Variant 2: with dilp data already at morphotype (species×site) level, the
# two-step specimen→species→site path collapses to the same result as variant 1.
write.csv(dat_site, file = "data/dat_site_sp_zero.csv", row.names = FALSE)
cat("Wrote data/dat_site_sp_zero.csv (", nrow(dat_site), "sites)\n")

# Variant 3: morphotype → site without tooth-fill; untoothed species are
# excluded from tooth trait site means via na.rm = TRUE. Matches Peppe et al. (2011).
dat_site_peppe <- aggregate(dat_prefill[, agg_cols],
                            by  = list(Site = dat_prefill$Site),
                            FUN = mean, na.rm = TRUE)
write.csv(dat_site_peppe, file = "data/dat_site_peppe.csv", row.names = FALSE)
cat("Wrote data/dat_site_peppe.csv (", nrow(dat_site_peppe), "sites)\n")

if (BUILD_PHYLOGENY) {

# ==============================================================================
# 4. LOAD FULL WCVP TREE AND BUILD NAME TABLE
# ==============================================================================

tree <- ladderize(read.tree("data/best_wcvp.tre_dated"))

split_names <- sapply(tree$tip.label, function(x) strsplit(x, "_"))
for (i in seq_along(split_names)) {
  if (length(split_names[[i]]) > 4)
    split_names[[i]] <- c(split_names[[i]][1:3],
                          paste(split_names[[i]][4:length(split_names[[i]])], collapse = "_"))
}
name_table <- as.data.frame(do.call(rbind, split_names), stringsAsFactors = FALSE)
colnames(name_table) <- c("order", "family", "genus", "species")

# Relabel to genus_species; orig_labels is a fixed-length reference used for
# all taxonomy lookups — never modified as the tree grows
tree$tip.label <- paste(name_table$genus, name_table$species, sep = "_")
orig_labels    <- tree$tip.label

# Save name_table so downstream scripts don't need to reload the full tree
write.csv(name_table, file = "data/name_table_full.csv", row.names = FALSE)

# ==============================================================================
# 5. HELPER: CROWN TIPS OF A CLADE
# ==============================================================================

# Given a set of tip labels that form a clade, returns one tip from each
# daughter of their MRCA — the pair that anchors the correct crown age.
crown_tips <- function(tr, tip_labels) {
  tip_idx <- which(tr$tip.label %in% tip_labels)
  if (length(tip_idx) == 0) return(character(0))
  if (length(tip_idx) == 1) return(tr$tip.label[tip_idx[1]])

  mrca_node <- getMRCA(tr, tip_idx)
  children  <- tr$edge[tr$edge[, 1] == mrca_node, 2]

  get_one_tip <- function(node) {
    while (node > length(tr$tip.label))
      node <- tr$edge[tr$edge[, 1] == node, 2][1]
    tr$tip.label[node]
  }
  unique(sapply(children, get_one_tip))
}

# ==============================================================================
# 6. BUILD FAMILY-LEVEL BACKBONE
# ==============================================================================

# For each angiosperm family in the WCVP tree, keep the 2 tips from the basal
# split. This ~2 × n_families tip tree is small enough that all subsequent
# getMRCA / bind.tip calls are fast, while still spanning the full angiosperm
# topology so any fossil order/family can be placed in 03_.

cat("Building family-level backbone...\n")
all_families  <- unique(name_table$family)
backbone_tips <- character(0)
for (fam in all_families) {
  backbone_tips <- c(backbone_tips,
                     crown_tips(tree, orig_labels[name_table$family == fam]))
}
backbone_tips <- unique(backbone_tips)
cat("  ", length(backbone_tips), "backbone tips across",
    length(all_families), "families\n")

tree_scaffold <- keep.tip(tree, backbone_tips)

# Taxonomy restricted to backbone tips — used for lookup during grafting
bb_idx              <- which(orig_labels %in% backbone_tips)
name_table_backbone <- name_table[bb_idx, ]
orig_labels_bb      <- orig_labels[bb_idx]

# ==============================================================================
# 7. GRAFT TRAINING SPECIES ONTO BACKBONE
# ==============================================================================

# All getMRCA / bind.tip calls now run on the small scaffold

h <- max(nodeHeights(tree_scaffold))

cat("Grafting", nrow(dat_sp), "training species onto backbone...\n")
for (i in seq_along(dat_sp$genusSpecies)) {
  cat("\r  ", round(i / nrow(dat_sp) * 100, 1), "%")
  current_sp <- dat_sp$genusSpecies[i]
  if (current_sp %in% tree_scaffold$tip.label) next

  target_node <- NULL

  # Genus: find backbone tips sharing the genus
  g_orig <- orig_labels_bb[name_table_backbone$genus == dat_sp$genus[i]]
  g_tips <- g_orig[g_orig %in% tree_scaffold$tip.label]
  if (length(g_tips) >= 2) {
    target_node <- getMRCA(tree_scaffold, g_tips)
  } else if (length(g_tips) == 1) {
    target_node <- which(tree_scaffold$tip.label == g_tips[1])
  }

  # Family fallback
  if (is.null(target_node)) {
    f_orig <- orig_labels_bb[name_table_backbone$family == dat_sp$Family[i]]
    f_tips <- f_orig[f_orig %in% tree_scaffold$tip.label]
    if (length(f_tips) >= 2) {
      target_node <- getMRCA(tree_scaffold, f_tips)
    } else if (length(f_tips) == 1) {
      target_node <- which(tree_scaffold$tip.label == f_tips[1])
    }
  }

  # Order fallback
  if (is.null(target_node)) {
    o_orig <- orig_labels_bb[name_table_backbone$order == dat_sp$Order[i]]
    o_tips <- o_orig[o_orig %in% tree_scaffold$tip.label]
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
    warning("Could not place training species: ", current_sp)
  }
}
cat("\n")

# ==============================================================================
# 8. WRITE OUTPUTS
# ==============================================================================

# tre_scaffold.tre — family backbone + all training species; used by
#   04_fossil_predictions.R to graft fossils into the full angiosperm topology
write.tree(tree_scaffold, file = "data/tre_scaffold.tre")
cat("Wrote data/tre_scaffold.tre (", length(tree_scaffold$tip.label), "tips)\n")

# tre_pruned.tre — training species only; used for VCV in PGLS / PIP
tree_pruned <- keep.tip(tree_scaffold, dat_sp$genusSpecies)
write.tree(tree_pruned, file = "data/tre_pruned.tre")
cat("Wrote data/tre_pruned.tre (", length(tree_pruned$tip.label), "tips)\n")

dat_pruned <- dat_sp[match(tree_pruned$tip.label, dat_sp$genusSpecies), ]
write.csv(dat_pruned, file = "data/data_species.csv", row.names = FALSE)
cat("Wrote data/data_species.csv (", nrow(dat_pruned), "species)\n")

} # end BUILD_PHYLOGENY
