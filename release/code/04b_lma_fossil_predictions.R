source("code/_setup.R")

library(ape)
library(phytools)
source("code/Phylogenetically-Informed_Predictions_Source.R")

PLACEMENT_FALLBACK <- "ancestral_branch"

# 1. LOAD COMPONENTS

pip    <- readRDS("models/lma_pip_components.rds")
nophy  <- readRDS("models/lma_nophy_models.rds")
foss   <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)

cat("Loaded", nrow(foss), "fossil species-site rows\n")

# pw2.a.ratio is the fossil equivalent of dilp's petiole_metric (PW^2/A, same units)
foss$log10_petiole_metric <- log10(foss$pw2.a.ratio)
foss_lma <- foss[!is.na(foss$log10_petiole_metric), ]
cat("Fossils with petiole metric:", nrow(foss_lma), "of", nrow(foss), "\n")

# 2. AGGREGATE: grand means for grafting, within-site means for prediction

# Grand means — one row per species; used only for tree grafting
foss_sp <- aggregate(
  foss_lma[, c("log10_petiole_metric", "age_ma")],
  by  = list(species = foss_lma$species,
             genus   = foss_lma$genus,
             family  = foss_lma$family,
             order   = foss_lma$order),
  FUN = mean, na.rm = TRUE
)

# Within-site means — one row per species × site; used to build X for prediction
foss_site <- aggregate(
  foss_lma[, "log10_petiole_metric", drop = FALSE],
  by  = list(species = foss_lma$species,
             site    = foss_lma$site,
             age_ma  = foss_lma$age_ma),
  FUN = mean, na.rm = TRUE
)

cat("Unique fossil species with petiole metric:", nrow(foss_sp), "\n")
cat("Fossil species-site rows:", nrow(foss_site), "\n")

# 3. LM PREDICTIONS

lm_preds <- predict(nophy$model,
                    newdata = foss_site[, "log10_petiole_metric", drop = FALSE])
foss_site$lma_lm <- 10^as.numeric(lm_preds)

site_lm <- aggregate(lma_lm ~ site + age_ma, data = foss_site, FUN = mean)

cat("LM predictions done\n")

# 4. GRAFT FOSSIL SPECIES ONTO LMA TRAINING TREE

# Graft fossils onto tre_lma_pruned.tre. The tree already contains all LMA
# training species; family/order placements use name_table_full for lookups.

tree     <- read.tree("data/tre_lma_pruned.tre")
h        <- max(nodeHeights(tree))
tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
name_tbl   <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)

for (j in seq_len(nrow(foss_sp))) {
  fossil_name <- foss_sp$species[j]
  age_ma      <- foss_sp$age_ma[j]

  if (fossil_name %in% tree$tip.label) next

  target_node <- NULL

  # Genus match within LMA tree
  genus_tips <- tree$tip.label[tip_genera == foss_sp$genus[j]]
  if (length(genus_tips) >= 2) {
    target_node <- getMRCA(tree, genus_tips)
  } else if (length(genus_tips) == 1) {
    target_node <- which(tree$tip.label == genus_tips[1])
  }

  # Family fallback
  if (is.null(target_node) && nzchar(foss_sp$family[j]) &&
      foss_sp$family[j] != "unknown") {
    fam_genera <- unique(name_tbl$genus[name_tbl$family == foss_sp$family[j]])
    fam_tips   <- tree$tip.label[tip_genera %in% fam_genera]
    if (length(fam_tips) >= 2) target_node <- getMRCA(tree, fam_tips)
    else if (length(fam_tips) == 1) target_node <- which(tree$tip.label == fam_tips[1])
  }

  # Order fallback
  if (is.null(target_node) && nzchar(foss_sp$order[j]) &&
      foss_sp$order[j] != "unknown") {
    ord_genera <- unique(name_tbl$genus[name_tbl$order == foss_sp$order[j]])
    ord_tips   <- tree$tip.label[tip_genera %in% ord_genera]
    if (length(ord_tips) >= 2) target_node <- getMRCA(tree, ord_tips)
    else if (length(ord_tips) == 1) target_node <- which(tree$tip.label == ord_tips[1])
  }

  # Root fallback
  if (is.null(target_node)) {
    warning("No taxonomy match for '", fossil_name, "'. Placing at root.")
    target_node <- length(tree$tip.label) + 1L
  }

  node_ht  <- nodeheight(tree, target_node)
  edge_len <- (h - age_ma) - node_ht

  if (edge_len >= 0) {
    tree <- bind.tip(tree, tip.label = fossil_name,
                     where = target_node, edge.length = edge_len)
  } else if (PLACEMENT_FALLBACK == "ancestral_branch") {
    fossil_ht <- h - age_ma
    node      <- target_node
    found     <- FALSE
    repeat {
      parent_idx <- which(tree$edge[, 2] == node)
      if (length(parent_idx) == 0) break
      parent_node <- tree$edge[parent_idx, 1]
      parent_ht   <- nodeheight(tree, parent_node)
      if (parent_ht <= fossil_ht) {
        edge_idx <- which(tree$edge[, 2] == node)
        max_pos  <- if (length(edge_idx) > 0) tree$edge.length[edge_idx] else 0
        position <- min(nodeheight(tree, node) - fossil_ht, max_pos)
        tree  <- bind.tip(tree, tip.label = fossil_name,
                          where = node, position = position, edge.length = 0)
        found <- TRUE
        break
      }
      node <- parent_node
    }
    if (!found) {
      warning("No spanning branch for '", fossil_name, "'. Placing at root.")
      tree <- bind.tip(tree, tip.label = fossil_name,
                       where = length(tree$tip.label) + 1L, edge.length = 0.001)
    }
  } else {
    tree <- bind.tip(tree, tip.label = fossil_name,
                     where = target_node, edge.length = 0.001)
  }

  tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
}

placed   <- foss_sp$species[foss_sp$species %in% tree$tip.label]
cat("Placed", length(placed), "/", nrow(foss_sp), "fossil species on tree\n")
foss_sp  <- foss_sp[foss_sp$species %in% placed, ]

# 5. COMPUTE VCV CROSS-COVARIANCES

idx_extant <- rownames(pip$dat_fit)   # LMA training species used in PGLS
idx_fossil  <- foss_sp$species

tree_small     <- keep.tip(tree, c(idx_extant, idx_fossil))
phylomat_small <- vcv(tree_small)

V_cross <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda
V_inv   <- solve(pip$V_lam)
resid   <- pip$resid[idx_extant]

phylo_adj <- as.numeric(t(V_cross) %*% V_inv %*% resid)
names(phylo_adj) <- idx_fossil

# 6. PIP PREDICTIONS

# Build design matrix from species-within-site means; phylo_adj is per species
placed_site <- foss_site[foss_site$species %in% placed, ]

X_fossil <- cbind(1, placed_site$log10_petiole_metric)
colnames(X_fossil) <- colnames(pip$X)
rownames(X_fossil) <- placed_site$species

# GLS prediction + phylogenetic correction (adjustment indexed by species)
yhat_log10 <- as.numeric(X_fossil %*% t(pip$beta)) + phylo_adj[placed_site$species]
placed_site$lma_pip <- 10^yhat_log10

site_pip <- aggregate(lma_pip ~ site + age_ma, data = placed_site, FUN = mean)

cat("PIP predictions done\n")

# 7. COMPILE AND WRITE OUTPUT

site_out <- Reduce(function(x, y) merge(x, y, by = c("site", "age_ma")),
                   list(site_lm, site_pip))
site_out <- site_out[order(-site_out$age_ma), ]
site_out[, c("lma_lm", "lma_pip")] <- round(site_out[, c("lma_lm", "lma_pip")], 1)

cat("\nSite-level LMA predictions (g/m²):\n")
print(site_out, row.names = FALSE)

# Per-species-site predictions
sp_out <- placed_site[, c("species", "site", "age_ma", "lma_pip")]
sp_out <- merge(sp_out,
                foss_site[foss_site$species %in% placed, c("species", "site", "lma_lm")],
                by = c("species", "site"))
sp_out <- sp_out[order(-sp_out$age_ma, sp_out$site, sp_out$species), ]

write.csv(site_out, "tables/lma_fossil_site_predictions.csv", row.names = FALSE)
write.csv(sp_out,   "tables/lma_fossil_species_predictions.csv", row.names = FALSE)
cat("\nSaved tables/lma_fossil_site_predictions.csv\n")
cat("Saved tables/lma_fossil_species_predictions.csv\n")
