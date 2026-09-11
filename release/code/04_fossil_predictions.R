# ==============================================================================
# 04_fossil_predictions.R
# Predict fossil site climate using local traits only:
#  LM site : average species-within-site traits, then predict site climate.
#  PIP     : predict each species-site occurrence using its local traits and
#            site age, then average predictions within the site.
# Fossil traits and ages must never be pooled across sites.
#
# Requires:
#   models/pip_components.rds   -- output of 02_phy_regression.R
#   models/site_models.rds      -- output of 01_nophy_regression.R
#   data/tre_scaffold.tre       -- output of 00_data_cleaning.R
#   data/fossil_traits.csv      -- species-site rows with trait + age_ma columns
#
# Output:
#   tables/fossil_site_comparison.csv       -- formal-only site comparison
#   tables/fossil_predictions.csv           -- formal-only species predictions
#   tables/fossil_taxonomy_sensitivity.csv  -- scenario deltas by site
#   tables/*_<scenario>.csv                  -- scenario predictions + placements
# ==============================================================================

source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(ape)
library(phytools)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")
source("code/fossil_taxonomy.R")
source("code/fossil_placement.R")

# ==============================================================================
# USER SETTINGS
# ==============================================================================

PLACEMENT_FALLBACK <- "ancestral_branch"

# ==============================================================================
# 1. LOAD COMPONENTS
# ==============================================================================

pip       <- readRDS("models/pip_components.rds")
site_mods <- readRDS("models/site_models.rds")
foss_base <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)

# Dana's requested non-phylogenetic comparison is the species-to-site,
# zero-fill, bag-imputed LM.  Select it explicitly instead of relying on the
# backwards-compatible specimen-impute alias.
site_lm_config <- site_mods$configs$sp_zero_impute
if (is.null(site_lm_config)) {
  stop("models/site_models.rds lacks the required sp_zero_impute site LM")
}

cat("Loaded", nrow(foss_base), "fossil species-site rows\n")

trait_cols <- site_lm_config$pred_names   # fossil-measurable site LM traits

run_fossil_scenario <- function(taxonomy_scenario) {
foss <- taxonomy_for_scenario(foss_base, taxonomy_scenario)
cat("\n=== Taxonomy scenario:", taxonomy_scenario, "===\n")

# ==============================================================================
# 2. PREPARE AGGREGATIONS
# ==============================================================================

# Site means (average all species-site rows within each site)
foss_site_means <- aggregate(foss[, trait_cols],
                             by  = list(site   = foss$site,
                                        age_ma = foss$age_ma),
                             FUN = mean, na.rm = TRUE)

# ==============================================================================
# 3. LM SITE-LEVEL
# ==============================================================================

pred_names_site <- site_lm_config$pred_names
foss_site_imp   <- predict(site_lm_config$impute_model,
                           newdata = foss_site_means[, pred_names_site])

site_lm_mat <- predict(site_lm_config$mat$LM,     newdata = foss_site_imp)
site_lm_map <- predict(site_lm_config$log_map$LM, newdata = foss_site_imp)

site_lm_site <- data.frame(
  site        = foss_site_means$site,
  age_ma      = foss_site_means$age_ma,
  mat_lm_site = as.numeric(site_lm_mat),
  map_lm_site = exp(as.numeric(site_lm_map))
)

cat("LM site-level predictions done\n")

# ==============================================================================
# 4. GRAFT FOSSIL TIPS ONTO SCAFFOLD TREE
# ==============================================================================
#
# One tip per species-site occurrence, placed at that site's age.
# Placement failures are logged and excluded downstream.

tree <- read.tree("data/tre_scaffold.tre")
# This snapshot is the only permitted taxonomic anchor set.  `tree` grows
# during placement, but an earlier fossil must not affect a later placement.
scaffold_tip_labels <- tree$tip.label
name_tbl <- pip$name_table_full

placement_rows_site <- vector("list", nrow(foss))
for (j in seq_len(nrow(foss))) {
  res <- graft_fossil_tip(
    tree, scaffold_tip_labels, foss$fossil_name[j], foss$age_ma[j],
    foss$genus[j], foss$family[j], foss$order[j],
    placement_fallback = PLACEMENT_FALLBACK, name_table = name_tbl
  )
  tree <- res$tree
  placement_rows_site[[j]] <- data.frame(
    taxonomy_scenario = taxonomy_scenario,
    tip_scope = "site_occurrence",
    species = foss$species[j],
    fossil_name = foss$fossil_name[j],
    site = foss$site[j],
    placed = res$placed,
    placement_level = res$placement_level,
    placement_target = res$placement_target,
    age_fallback = res$age_fallback,
    error = res$error,
    stringsAsFactors = FALSE
  )
}

placed_occ  <- foss$fossil_name[vapply(placement_rows_site, function(r) r$placed, logical(1))]
dropped_occ <- setdiff(foss$fossil_name, placed_occ)
cat("PIP: placed", length(placed_occ), "/", nrow(foss), "fossil species-site rows on tree\n")
if (length(dropped_occ) > 0) {
  dropped_rows <- foss[foss$fossil_name %in% dropped_occ, c("fossil_name", "species", "site")]
  warning("Dropping ", length(dropped_occ), " fossil species-site occurrences ",
          "that failed phylogenetic placement (PIP): ",
          paste(dropped_rows$fossil_name, collapse = ", "))
  cat("  Dropped species-site rows (sites affected):",
      paste(unique(dropped_rows$site), collapse = ", "), "\n")
}
foss <- foss[foss$fossil_name %in% placed_occ, ]

placement_log <- do.call(rbind, placement_rows_site)

# ==============================================================================
# 5. COMPUTE VCV CROSS-COVARIANCES
# ==============================================================================

idx_extant  <- rownames(pip$dat_imputed_mat)
idx_fossil_occ <- foss$fossil_name    # site-occurrence tips (PIP, issue #11)

tree_small     <- keep.tip(tree, c(idx_extant, idx_fossil_occ))
phylomat_small <- vcv(tree_small)

V_inv_mat <- solve(pip$V_lam_mat)
V_inv_map <- solve(pip$V_lam_map)

resid_ord_mat <- pip$resid_mat[idx_extant]
resid_ord_map <- pip$resid_map[idx_extant]

compute_phylo_adj <- function(idx_fossil) {
  V_cross_mat <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_mat
  V_cross_map <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda_map
  adj_mat <- as.numeric(t(V_cross_mat) %*% V_inv_mat %*% resid_ord_mat)
  adj_map <- as.numeric(t(V_cross_map) %*% V_inv_map %*% resid_ord_map)
  names(adj_mat) <- idx_fossil
  names(adj_map) <- idx_fossil
  list(mat = adj_mat, map = adj_map)
}

# Each occurrence receives its own phylogenetic adjustment.
phylo_adj_occ     <- compute_phylo_adj(idx_fossil_occ)
phylo_adj_site_mat <- phylo_adj_occ$mat
phylo_adj_site_map <- phylo_adj_occ$map

# ==============================================================================
# 6. IMPUTE FOSSIL TRAITS FOR PIP
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

# ==============================================================================
# 7. PIP OCCURRENCE PREDICTIONS
# ==============================================================================

# For each species-site row, use site-specific traits and the phylogenetic
# adjustment from that occurrence's own placement.

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

# GLS prediction uses site-specific traits; phylo adjustment from this
# occurrence's own site-specific tip placement.
yhat_site_mat <- as.numeric(X_site_mat %*% pip$beta_mat) +
                   phylo_adj_site_mat[foss$fossil_name]
yhat_site_map <- as.numeric(X_site_map %*% pip$beta_map) +
                   phylo_adj_site_map[foss$fossil_name]

foss$mat_pip_site <- yhat_site_mat
foss$map_pip_site <- exp(yhat_site_map)

site_pip_site <- aggregate(cbind(mat_pip_site, map_pip_site) ~ site + age_ma,
                            data = foss, FUN = mean)
# MAP is the geometric mean of occurrence predictions; MAT remains arithmetic.
site_pip_site$map_pip_site <- exp(aggregate(log(map_pip_site) ~ site + age_ma,
  data = foss, FUN = mean)[[3]])

# Also save per-species PIP predictions
results_per_species <- foss[, c(
  "fossil_name", "species", "site", "age_ma",
  "genus_reported", "family_reported", "order_reported",
  "genus_informal", "family_informal", "order_informal",
  "mat_pip_site", "map_pip_site"
)]
results_per_species$taxonomy_scenario <- taxonomy_scenario
scenario_predictions_path <- file.path(
  "tables", paste0("fossil_predictions_", taxonomy_scenario, ".csv")
)
write.csv(results_per_species, scenario_predictions_path, row.names = FALSE)
cat("PIP predictions done\n")

# ==============================================================================
# 8. OUTPUT COMPARISON TABLE (LONG FORMAT)
# ==============================================================================

site_comparison <- Reduce(
  function(x, y) merge(x, y, by = c("site", "age_ma")),
  list(site_lm_site, site_pip_site)
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
                          "lm_site", "pip_site")]
cat("\nLong-form site comparison:\n")
print(site_long, row.names = FALSE)
site_long$taxonomy_scenario <- taxonomy_scenario
scenario_site_path <- file.path(
  "tables", paste0("fossil_site_comparison_", taxonomy_scenario, ".csv")
)
scenario_placement_path <- file.path(
  "tables", paste0("fossil_placement_log_", taxonomy_scenario, ".csv")
)
write.csv(site_long, scenario_site_path, row.names = FALSE)
write.csv(placement_log, scenario_placement_path, row.names = FALSE)
cat("\nResults written to", scenario_site_path, "\n")

list(
  site = site_long,
  species = results_per_species,
  placement = placement_log
)
}

# ==============================================================================
# 9. RUN PRIMARY + SENSITIVITY SCENARIOS
# ==============================================================================

scenario_names <- c("formal_only", "include_informal")
scenario_results <- setNames(
  lapply(scenario_names, run_fossil_scenario),
  scenario_names
)

# Backward-compatible outputs use Dana's requested formal-only interpretation.
write.csv(
  scenario_results$formal_only$species,
  "tables/fossil_predictions.csv",
  row.names = FALSE
)
write.csv(
  scenario_results$formal_only$site,
  "tables/fossil_site_comparison.csv",
  row.names = FALSE
)

formal_site <- scenario_results$formal_only$site
informal_site <- scenario_results$include_informal$site
sensitivity <- merge(
  formal_site,
  informal_site,
  by = c("site", "age_ma", "variable"),
  suffixes = c("_formal_only", "_include_informal")
)

metric_cols <- c("lm_site", "pip_site")
for (metric in metric_cols) {
  sensitivity[[paste0(metric, "_delta")]] <-
    sensitivity[[paste0(metric, "_include_informal")]] -
    sensitivity[[paste0(metric, "_formal_only")]]
  sensitivity[[paste0(metric, "_delta")]] <- ifelse(
    sensitivity$variable == "MAT",
    round(sensitivity[[paste0(metric, "_delta")]], 1),
    round(sensitivity[[paste0(metric, "_delta")]], 0)
  )
}
sensitivity <- sensitivity[, c(
  "site", "age_ma", "variable",
  as.vector(rbind(
    paste0(metric_cols, "_formal_only"),
    paste0(metric_cols, "_include_informal"),
    paste0(metric_cols, "_delta")
  ))
)]
write.csv(
  sensitivity,
  "tables/fossil_taxonomy_sensitivity.csv",
  row.names = FALSE
)

cat("\nPrimary formal-only outputs also written to the legacy filenames.\n")
cat("Sensitivity comparison written to",
    "tables/fossil_taxonomy_sensitivity.csv\n")
