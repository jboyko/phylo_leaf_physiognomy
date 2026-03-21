setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(phytools)

# ==============================================================================
# 1. LOAD AND CLEAN SPECIMEN-LEVEL DATA
# ==============================================================================

dat <- read.csv("data/RoyerLeafShapeClimateDataFixedNames_June2012.csv")
dat$genusSpecies[dat$genusSpecies == " Dialyanthera sp."] <- "Dialyanthera sp."
dat$genusSpecies <- gsub(" ", "_", dat$genusSpecies)
genera_names <- sapply(strsplit(dat$genusSpecies, "_"), `[`, 1)

# ==============================================================================
# 2. FILL TOOTH TRAITS FOR UNTOOTHED LEAVES (Dana Royer, pers. comm.)
# ==============================================================================

# For confirmed untoothed leaves (blank AND margin.score == 1), tooth trait
# absence is real — set to 0 (or 1 for perim.ratio) BEFORE aggregating so
# species means are biologically correct rather than NaN.

vars_tooth_0 <- c("primary.teeth.number", "secondary.teeth.number",
                   "teeth.number", "teeth.perimeter.percm",
                   "teeth.interior.percm", "tooth.area.cm2",
                   "avt.tooth.area", "tooth.area.perimeter",
                   "tooth.area.interior", "tooth.area.blade.area.ratio",
                   "teeth.blade.area.ratio")

untoothed <- !is.na(dat$margin.score) & dat$margin.score == 1
for (v in vars_tooth_0) dat[[v]][untoothed & is.na(dat[[v]])] <- 0
dat[["perim.ratio"]][untoothed & is.na(dat[["perim.ratio"]])] <- 1

# ==============================================================================
# 3. AGGREGATE TO SPECIES LEVEL
# ==============================================================================

dat_sp <- aggregate(dat[, 9:ncol(dat)],
  by = list(dat$genusSpecies, genera_names, dat$Family, dat$Order),
  FUN = mean, na.rm = TRUE)
dat_sp <- dat_sp[dat_sp$Group.4 != "unknown", ]

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

# All getMRCA / bind.tip calls now run on the small scaffold rather than the
# full 123k-tip tree. h is constant because bind.tip always reaches crown age.

h <- max(nodeHeights(tree_scaffold))

cat("Grafting", nrow(dat_sp), "training species onto backbone...\n")
for (i in seq_along(dat_sp$Group.1)) {
  cat("\r  ", round(i / nrow(dat_sp) * 100, 1), "%")
  current_sp <- dat_sp$Group.1[i]
  if (current_sp %in% tree_scaffold$tip.label) next

  target_node <- NULL

  # Genus: find backbone tips sharing the genus
  g_orig <- orig_labels_bb[name_table_backbone$genus == dat_sp$Group.2[i]]
  g_tips <- g_orig[g_orig %in% tree_scaffold$tip.label]
  if (length(g_tips) >= 2) {
    target_node <- getMRCA(tree_scaffold, g_tips)
  } else if (length(g_tips) == 1) {
    target_node <- which(tree_scaffold$tip.label == g_tips[1])
  }

  # Family fallback
  if (is.null(target_node)) {
    f_orig <- orig_labels_bb[name_table_backbone$family == dat_sp$Group.3[i]]
    f_tips <- f_orig[f_orig %in% tree_scaffold$tip.label]
    if (length(f_tips) >= 2) {
      target_node <- getMRCA(tree_scaffold, f_tips)
    } else if (length(f_tips) == 1) {
      target_node <- which(tree_scaffold$tip.label == f_tips[1])
    }
  }

  # Order fallback
  if (is.null(target_node)) {
    o_orig <- orig_labels_bb[name_table_backbone$order == dat_sp$Group.4[i]]
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
#   03_fossil_predictions.R to graft fossils into the full angiosperm topology
write.tree(tree_scaffold, file = "data/tre_scaffold.tre")
cat("Wrote data/tre_scaffold.tre (", length(tree_scaffold$tip.label), "tips)\n")

# tre_pruned.tre — training species only; used for VCV in PGLS / PIP
tree_pruned <- keep.tip(tree_scaffold, dat_sp$Group.1)
write.tree(tree_pruned, file = "data/tre_pruned.tre")
cat("Wrote data/tre_pruned.tre (", length(tree_pruned$tip.label), "tips)\n")

dat_pruned <- dat_sp[match(tree_pruned$tip.label, dat_sp$Group.1), ]
write.csv(dat_pruned, file = "data/data_species.csv", row.names = FALSE)
cat("Wrote data/data_species.csv (", nrow(dat_pruned), "species)\n")
