# ==============================================================================
# 04_fossil_predictions.R
# Predict MAT and log(MAP) for fossil sites using four methods:
#
#  LM sp    : LM trained on extant species means; fossil species grand means
#             predicted per species, averaged to site
#  LM site  : LM trained on extant site means; fossil site-mean traits
#             predicted directly per site
#  PIP sp   : PIP using fossil species grand means, averaged to site
#  PIP+site : PIP using site-specific fossil species means, averaged to site
#
# PIP sp and PIP+site share the same phylogenetic adjustment (same tip
# placement); they differ only in the trait values used for the GLS prediction.
#
# Requires:
#   models/pip_components.rds   -- output of 02_phy_regression.R
#   models/nophy_models.rds     -- output of 01_nophy_regression.R
#   models/site_models.rds      -- output of 01_nophy_regression.R
#   data/tre_scaffold.tre       -- output of 00_data_cleaning.R
#   data/fossil_traits.csv      -- species-site rows with trait + age_ma columns
#
# Output:
#   tables/fossil_site_comparison.csv   -- site-level 4-model comparison
#   tables/fossil_predictions.csv       -- per-species PIP predictions
# ==============================================================================

setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(phytools)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# USER SETTINGS
# ==============================================================================

PLACEMENT_FALLBACK <- "ancestral_branch"

# ==============================================================================
# 1. LOAD COMPONENTS
# ==============================================================================

pip       <- readRDS("models/pip_components.rds")
nophy     <- readRDS("models/nophy_models.rds")
site_mods <- readRDS("models/site_models.rds")
foss      <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)

cat("Loaded", nrow(foss), "fossil species-site rows\n")

trait_cols <- nophy$pred_names   # 12 fossil-measurable traits

# ==============================================================================
# 2. PREPARE AGGREGATIONS
# ==============================================================================

# 2a. Species grand means (average across all sites a species appears in)
foss_sp <- aggregate(foss[, c(trait_cols, "age_ma")],
                     by  = list(species = foss$species,
                                genus   = foss$genus,
                                family  = foss$family,
                                order   = foss$order),
                     FUN = mean, na.rm = TRUE)

# 2b. Site means (average all species-site rows within each site)
foss_site_means <- aggregate(foss[, trait_cols],
                             by  = list(site   = foss$site,
                                        age_ma = foss$age_ma),
                             FUN = mean, na.rm = TRUE)

# Species → sites lookup (needed to aggregate species predictions back to sites)
sp_site_map <- unique(foss[, c("species", "site", "age_ma")])

# ==============================================================================
# 3. LM SPECIES-LEVEL
# ==============================================================================

# Impute using species-level training imputation model
foss_sp_imp <- predict(nophy$impute_model, newdata = foss_sp[, trait_cols])

lm_sp_mat <- predict(nophy$mat$LM,      newdata = foss_sp_imp)
lm_sp_map <- predict(nophy$log_map$LM,  newdata = foss_sp_imp)

foss_sp$mat_lm_sp  <- as.numeric(lm_sp_mat)
foss_sp$map_lm_sp  <- exp(as.numeric(lm_sp_map))

# Propagate species predictions back to site rows, then average per site
sp_lm <- merge(sp_site_map,
               foss_sp[, c("species", "mat_lm_sp", "map_lm_sp")],
               by = "species")
site_lm_sp <- aggregate(cbind(mat_lm_sp, map_lm_sp) ~ site + age_ma,
                         data = sp_lm, FUN = mean)

cat("LM species-level predictions done\n")

# ==============================================================================
# 4. LM SITE-LEVEL
# ==============================================================================

pred_names_site <- site_mods$pred_names
foss_site_imp   <- predict(site_mods$impute_model,
                           newdata = foss_site_means[, pred_names_site])

site_lm_mat <- predict(site_mods$mat$LM,     newdata = foss_site_imp)
site_lm_map <- predict(site_mods$log_map$LM, newdata = foss_site_imp)

site_lm_site <- data.frame(
  site        = foss_site_means$site,
  age_ma      = foss_site_means$age_ma,
  mat_lm_site = as.numeric(site_lm_mat),
  map_lm_site = exp(as.numeric(site_lm_map))
)

cat("LM site-level predictions done\n")

# ==============================================================================
# 5. GRAFT UNIQUE FOSSIL SPECIES ONTO SCAFFOLD TREE
# ==============================================================================

tree     <- read.tree("data/tre_scaffold.tre")
h        <- max(nodeHeights(tree))
tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
name_tbl   <- pip$name_table_full

for (j in seq_len(nrow(foss_sp))) {
  fossil_name <- foss_sp$species[j]
  age_ma      <- foss_sp$age_ma[j]

  if (fossil_name %in% tree$tip.label) next   # already placed

  target_node <- NULL

  # Genus match
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
        position <- nodeheight(tree, node) - fossil_ht
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

placed <- foss_sp$species[foss_sp$species %in% tree$tip.label]
cat("Placed", length(placed), "/", nrow(foss_sp), "fossil species on tree\n")
foss_sp <- foss_sp[foss_sp$species %in% placed, ]

# ==============================================================================
# 6. COMPUTE VCV CROSS-COVARIANCES
# ==============================================================================

idx_extant <- rownames(pip$dat_imputed_mat)
idx_fossil  <- foss_sp$species

tree_small     <- keep.tip(tree, c(idx_extant, idx_fossil))
phylomat_small <- vcv(tree_small)

V_cross_mat <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_mat
V_cross_map <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_map

V_inv_mat <- solve(pip$V_lam_mat)
V_inv_map <- solve(pip$V_lam_map)

resid_ord_mat <- pip$resid_mat[idx_extant]
resid_ord_map <- pip$resid_map[idx_extant]

# Phylogenetic adjustment — same for PIP sp and PIP+site
phylo_adj_mat <- as.numeric(t(V_cross_mat) %*% V_inv_mat %*% resid_ord_mat)
phylo_adj_map <- as.numeric(t(V_cross_map) %*% V_inv_map %*% resid_ord_map)
names(phylo_adj_mat) <- idx_fossil
names(phylo_adj_map) <- idx_fossil

# ==============================================================================
# 7. IMPUTE FOSSIL TRAITS FOR PIP
# ==============================================================================

pad_cols <- function(df, cols) {
  out <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(cols),
                               dimnames = list(NULL, cols)))
  shared <- intersect(cols, names(df))
  out[, shared] <- df[, shared]
  out
}

trait_cols_mat <- all.vars(pip$formula_mat)[-1]
trait_cols_map <- all.vars(pip$formula_map)[-1]

# Species grand means for PIP sp
sp_mat_raw <- pad_cols(foss_sp, trait_cols_mat)
sp_map_raw <- pad_cols(foss_sp, trait_cols_map)
sp_mat_imp <- predict(pip$impute_model_mat_traits, newdata = sp_mat_raw)
sp_map_imp <- predict(pip$impute_model_map_traits, newdata = sp_map_raw)
rownames(sp_mat_imp) <- foss_sp$species
rownames(sp_map_imp) <- foss_sp$species

# ==============================================================================
# 8. BUILD DESIGN MATRICES FOR PIP SP
# ==============================================================================

sp_mat_imp$mat     <- 0
sp_map_imp$log_map <- 0
X_sp_mat <- model.matrix(pip$formula_mat, sp_mat_imp)[, colnames(pip$X_mat), drop = FALSE]
X_sp_map <- model.matrix(pip$formula_map, sp_map_imp)[, colnames(pip$X_map), drop = FALSE]

# ==============================================================================
# 9. PIP SPECIES-LEVEL PREDICTIONS
# ==============================================================================

yhat_sp_mat <- as.numeric(X_sp_mat %*% pip$beta_mat) + phylo_adj_mat[foss_sp$species]
yhat_sp_map <- as.numeric(X_sp_map %*% pip$beta_map) + phylo_adj_map[foss_sp$species]

foss_sp$mat_pip_sp  <- yhat_sp_mat
foss_sp$map_pip_sp  <- exp(yhat_sp_map)

sp_pip <- merge(sp_site_map,
                foss_sp[, c("species", "mat_pip_sp", "map_pip_sp")],
                by = "species")
site_pip_sp <- aggregate(cbind(mat_pip_sp, map_pip_sp) ~ site + age_ma,
                          data = sp_pip, FUN = mean)

cat("PIP species-level predictions done\n")

# ==============================================================================
# 10. PIP+SITE PREDICTIONS
# ==============================================================================

# For each species-site row, use site-specific trait values with the
# species-level phylogenetic adjustment

site_sp_mat_raw <- pad_cols(foss, trait_cols_mat)
site_sp_map_raw <- pad_cols(foss, trait_cols_map)
site_sp_mat_imp <- predict(pip$impute_model_mat_traits, newdata = site_sp_mat_raw)
site_sp_map_imp <- predict(pip$impute_model_map_traits, newdata = site_sp_map_raw)
rownames(site_sp_mat_imp) <- foss$fossil_name
rownames(site_sp_map_imp) <- foss$fossil_name

site_sp_mat_imp$mat     <- 0
site_sp_map_imp$log_map <- 0
X_site_mat <- model.matrix(pip$formula_mat, site_sp_mat_imp)[, colnames(pip$X_mat), drop = FALSE]
X_site_map <- model.matrix(pip$formula_map, site_sp_map_imp)[, colnames(pip$X_map), drop = FALSE]

# GLS prediction uses site-specific traits; phylo adjustment from species placement
yhat_site_mat <- as.numeric(X_site_mat %*% pip$beta_mat) +
                   phylo_adj_mat[foss$species]
yhat_site_map <- as.numeric(X_site_map %*% pip$beta_map) +
                   phylo_adj_map[foss$species]

foss$mat_pip_site <- yhat_site_mat
foss$map_pip_site <- exp(yhat_site_map)

site_pip_site <- aggregate(cbind(mat_pip_site, map_pip_site) ~ site + age_ma,
                            data = foss, FUN = mean)

# Also save per-species PIP+site predictions
results_per_species <- foss[, c("fossil_name", "species", "site", "age_ma",
                                 "mat_pip_site", "map_pip_site")]
write.csv(results_per_species, "tables/fossil_predictions.csv", row.names = FALSE)
cat("PIP+site predictions done\n")

# ==============================================================================
# 11. OUTPUT COMPARISON TABLE (LONG FORMAT)
# ==============================================================================

site_comparison <- Reduce(
  function(x, y) merge(x, y, by = c("site", "age_ma")),
  list(site_lm_sp, site_lm_site, site_pip_sp, site_pip_site)
)
site_comparison <- site_comparison[order(-site_comparison$age_ma), ]

# Round for readability
mat_cols <- grep("^mat_", names(site_comparison), value = TRUE)
map_cols <- grep("^map_", names(site_comparison), value = TRUE)
site_comparison[, mat_cols] <- round(site_comparison[, mat_cols], 1)
site_comparison[, map_cols] <- round(site_comparison[, map_cols], 0)

mat_df <- site_comparison[, c("site", "age_ma", mat_cols)]
names(mat_df) <- sub("^mat_", "", names(mat_df))
mat_df$variable <- "MAT"
map_df <- site_comparison[, c("site", "age_ma", map_cols)]
names(map_df) <- sub("^map_", "", names(map_df))
map_df$variable <- "MAP"
site_long <- rbind(mat_df, map_df)
site_long <- site_long[, c("site", "age_ma", "variable",
                          "lm_sp", "lm_site", "pip_sp", "pip_site")]
cat("\nLong-form site comparison:\n")
print(site_long, row.names = FALSE)
write.csv(site_long, "tables/fossil_site_comparison.csv", row.names = FALSE)
cat("\nResults written to tables/fossil_site_comparison.csv\n")
