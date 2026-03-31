# ==============================================================================
# 03_fossil_predictions.R
# Predict MAT and log(MAP) for fossil specimens using PIP
# (Phylogenetically-Informed Predictions)
#
# Scaffold: data/tre_scaffold.tre  -- family-level angiosperm backbone with all
#           training species grafted in. Fossils placed here so any angiosperm
#           family resolves to the correct phylogenetic neighbourhood.
#
# Key property: vcv(tree_pruned)[i,j] == vcv(tree_scaffold)[i,j] for any two
# training species, so pip$V_lam_mat / pip$V_lam_map are reused unchanged.
# Only the cross-covariance V[training, fossil] needs to be computed.
#
# Requires:
#   models/pip_components.rds   -- output of 02_phy_regression.R
#   data/tre_scaffold.tre       -- output of 00_data_cleaning.R
#   data/fossil_traits.csv      -- user-provided fossil data
#
# Output:
#   tables/fossil_predictions.csv
# ==============================================================================

setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(phytools)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# USER SETTINGS
# ==============================================================================

# What to do when a fossil predates the MRCA of its genus/family/order in the
# scaffold tree (i.e., the fossil is older than the crown age of its clade):
#
#   "ancestral_branch" — walk up the tree to find the branch that spans the
#                        fossil's age and split it there; tip terminates at the
#                        correct time depth (recommended: trusts the fossil age)
#
#   "node"             — attach directly at the MRCA node with a minimal edge
#                        (0.001 Ma); simpler, but ignores the fossil's true age

PLACEMENT_FALLBACK <- "ancestral_branch"

# ==============================================================================
# 1. LOAD COMPONENTS AND FOSSIL DATA
# ==============================================================================

pip  <- readRDS("models/pip_components.rds")
foss <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)

cat("Loaded", nrow(foss), "fossil specimens\n")

# ==============================================================================
# 2. GRAFT FOSSILS ONTO THE FULL BACKBONE TREE
# ==============================================================================

# Load full backbone: all tips in genus_species format, 97 angiosperm orders
tree <- read.tree("data/tre_scaffold.tre")
h    <- max(nodeHeights(tree))

# Genus lookup from tip labels (format: genus_species)
tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)

# Full-tree taxonomy: covers all WCVP genera (all 97 orders)
name_tbl <- pip$name_table_full   # columns: order, family, genus, species, tip_label

for (j in 1:nrow(foss)) {
  fossil_name <- foss$fossil_name[j]
  age_ma      <- foss$age_ma[j]

  target_node <- NULL

  # ── Genus match ───────────────────────────────────────────────────────────
  genus_tips <- tree$tip.label[tip_genera == foss$genus[j]]
  if (length(genus_tips) >= 2) {
    target_node <- getMRCA(tree, genus_tips)
  } else if (length(genus_tips) == 1) {
    target_node <- which(tree$tip.label == genus_tips[1])
  }

  # ── Family fallback ───────────────────────────────────────────────────────
  if (is.null(target_node)) {
    fam_genera <- unique(name_tbl$genus[name_tbl$family == foss$family[j]])
    fam_tips   <- tree$tip.label[tip_genera %in% fam_genera]
    if (length(fam_tips) >= 2) target_node <- getMRCA(tree, fam_tips)
    else if (length(fam_tips) == 1) target_node <- which(tree$tip.label == fam_tips[1])
  }

  # ── Order fallback ────────────────────────────────────────────────────────
  if (is.null(target_node)) {
    ord_genera <- unique(name_tbl$genus[name_tbl$order == foss$order[j]])
    ord_tips   <- tree$tip.label[tip_genera %in% ord_genera]
    if (length(ord_tips) >= 2) target_node <- getMRCA(tree, ord_tips)
    else if (length(ord_tips) == 1) target_node <- which(tree$tip.label == ord_tips[1])
  }

  # ── Fallback: root (phylo adjustment → 0, most conservative) ─────────────
  if (is.null(target_node)) {
    warning("Order '", foss$order[j], "' not found in full WCVP tree. ",
            "Placing fossil '", fossil_name, "' at tree root (conservative).")
    target_node <- length(tree$tip.label) + 1L   # root node in ape convention
  }

  # Edge length: tip terminates at fossil age before present
  node_ht  <- nodeheight(tree, target_node)
  edge_len <- (h - age_ma) - node_ht

  if (edge_len >= 0) {
    # Normal case: fossil is younger than its placement node
    tree <- bind.tip(tree, tip.label = fossil_name,
                     where = target_node, edge.length = edge_len)
    cat("Grafted '", fossil_name, "' onto ", foss$genus[j],
        " clade at ", age_ma, " Ma\n", sep = "")
  } else if (PLACEMENT_FALLBACK == "ancestral_branch") {
    # Fossil is older than its placement node — walk up toward the root to find
    # the branch that spans the fossil's age, then split it at the correct time.
    fossil_ht <- h - age_ma
    node      <- target_node
    found     <- FALSE
    repeat {
      parent_idx <- which(tree$edge[, 2] == node)
      if (length(parent_idx) == 0) break   # reached root without finding a span
      parent_node <- tree$edge[parent_idx, 1]
      parent_ht   <- nodeheight(tree, parent_node)
      if (parent_ht <= fossil_ht) {
        # fossil_ht falls on the edge from parent_node → node; split it here
        position <- nodeheight(tree, node) - fossil_ht   # distance back from child
        tree <- bind.tip(tree, tip.label = fossil_name,
                         where = node, position = position, edge.length = 0)
        message("Note: '", fossil_name, "' (", age_ma, " Ma) predates its ",
                foss$genus[j], " MRCA — attached on ancestral branch at correct age.")
        found <- TRUE
        break
      }
      node <- parent_node
    }
    if (!found) {
      warning("Could not find spanning branch for '", fossil_name,
              "'. Placing at root with minimal edge.")
      tree <- bind.tip(tree, tip.label = fossil_name,
                       where = length(tree$tip.label) + 1L, edge.length = 0.001)
    }
  } else {
    # PLACEMENT_FALLBACK == "node": attach at the MRCA with a minimal edge
    warning("Fossil '", fossil_name, "' (", age_ma, " Ma) predates its ",
            foss$genus[j], " MRCA — attaching at node with minimal edge (0.001 Ma).")
    tree <- bind.tip(tree, tip.label = fossil_name,
                     where = target_node, edge.length = 0.001)
  }

  # Update tip_genera for subsequent fossils
  tip_genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
}

# Fossils successfully placed
placed_fossils <- foss$fossil_name[foss$fossil_name %in% tree$tip.label]
if (length(placed_fossils) == 0) stop("No fossils could be placed on the tree.")
foss <- foss[foss$fossil_name %in% placed_fossils, ]

# ==============================================================================
# 3. EXTRACT CROSS-COVARIANCE FROM PRUNED SUBTREE
# ==============================================================================

idx_extant <- rownames(pip$dat_imputed_mat)   # training species names
idx_fossil  <- foss$fossil_name               # placed fossil names

# Prune expanded tree to training + fossils only — makes vcv() fast
tree_small     <- keep.tip(tree, c(idx_extant, idx_fossil))
phylomat_small <- vcv(tree_small)

# Cross-covariance V[training, fossil]: off-diagonal → scale by lambda only
# (the +1-lambda diagonal correction applies only to same-tip entries)
V_cross_mat <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_mat
V_cross_map <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_map

# Reuse pre-stored training VCV inverse (unchanged by what fossils are added)
V_inv_mat <- solve(pip$V_lam_mat)
V_inv_map <- solve(pip$V_lam_map)

# ==============================================================================
# 4. IMPUTE MISSING FOSSIL TRAITS
# ==============================================================================

# Helper: build a data frame with all required trait columns (NAs for missing)
pad_trait_cols <- function(fossil_df, required_cols) {
  out <- as.data.frame(matrix(NA, nrow = nrow(fossil_df),
                               ncol = length(required_cols),
                               dimnames = list(NULL, required_cols)))
  shared <- intersect(required_cols, names(fossil_df))
  out[, shared] <- fossil_df[, shared]
  return(out)
}

# MAT traits (use traits-only imputer — response is unknown for fossils)
trait_cols_mat  <- all.vars(pip$formula_mat)[-1]
fossil_mat_raw  <- pad_trait_cols(foss, trait_cols_mat)
fossil_mat_imp  <- predict(pip$impute_model_mat_traits, newdata = fossil_mat_raw)
rownames(fossil_mat_imp) <- foss$fossil_name

# MAP traits
trait_cols_map  <- all.vars(pip$formula_map)[-1]
fossil_map_raw  <- pad_trait_cols(foss, trait_cols_map)
fossil_map_imp  <- predict(pip$impute_model_map_traits, newdata = fossil_map_raw)
rownames(fossil_map_imp) <- foss$fossil_name

# ==============================================================================
# 5. BUILD FOSSIL DESIGN MATRICES
# ==============================================================================

# Add dummy response columns so model.matrix can parse the formula
fossil_mat_imp$mat     <- 0
fossil_map_imp$log_map <- 0

X_fossil_mat <- model.matrix(pip$formula_mat, fossil_mat_imp)
X_fossil_mat <- X_fossil_mat[, colnames(pip$X_mat), drop = FALSE]

X_fossil_map <- model.matrix(pip$formula_map, fossil_map_imp)
X_fossil_map <- X_fossil_map[, colnames(pip$X_map), drop = FALSE]

# ==============================================================================
# 6. APPLY PIP FORMULA
# ==============================================================================

# Reorder residuals to match VCV row order
resid_ord_mat <- pip$resid_mat[idx_extant]
resid_ord_map <- pip$resid_map[idx_extant]

phylo_adj_mat <- t(V_cross_mat) %*% V_inv_mat %*% resid_ord_mat   # m × 1
phylo_adj_map <- t(V_cross_map) %*% V_inv_map %*% resid_ord_map   # m × 1

yhat_mat <- as.numeric(X_fossil_mat %*% pip$beta_mat + phylo_adj_mat)
yhat_map <- as.numeric(X_fossil_map %*% pip$beta_map + phylo_adj_map)

# Prediction SE from coefficient uncertainty (design-matrix propagation)
se_mat <- sqrt(diag(X_fossil_mat %*% pip$vcv_mat %*% t(X_fossil_mat)))
se_map <- sqrt(diag(X_fossil_map %*% pip$vcv_map %*% t(X_fossil_map)))

# ==============================================================================
# 7. OUTPUT RESULTS
# ==============================================================================

results <- data.frame(
  fossil_name  = foss$fossil_name,
  mat_pred     = yhat_mat,
  mat_95lo     = yhat_mat - 1.96 * se_mat,
  mat_95hi     = yhat_mat + 1.96 * se_mat,
  log_map_pred = yhat_map,
  map_pred     = exp(yhat_map),
  map_95lo     = exp(yhat_map - 1.96 * se_map),
  map_95hi     = exp(yhat_map + 1.96 * se_map)
)

print(results)
write.csv(results, "tables/fossil_predictions.csv", row.names = FALSE)
cat("\nResults written to tables/fossil_predictions.csv\n")
