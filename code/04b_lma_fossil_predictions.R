source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(ape)
library(phytools)
source("code/Phylogenetically-Informed_Predictions_Source.R")
source("code/fossil_taxonomy.R")
source("code/fossil_placement.R")

PLACEMENT_FALLBACK <- "ancestral_branch"

# ==============================================================================
# 1. LOAD COMPONENTS
# ==============================================================================

pip_all   <- readRDS("models/lma_pip_components.rds")
nophy_all <- readRDS("models/lma_nophy_models.rds")
pip       <- pip_all$univariate
nophy     <- nophy_all$univariate
foss      <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)
foss      <- taxonomy_for_scenario(foss, "formal_only")

cat("Loaded", nrow(foss), "fossil species-site rows\n")

# pw2.a.ratio is the fossil equivalent of the LMA calibration predictor
# log10_pw2a_ratio (PW^2/A, same units).
foss$log10_pw2a_ratio <- log10(foss$pw2.a.ratio)
foss_lma <- foss[is.finite(foss$log10_pw2a_ratio), ]
cat("Fossils with petiole metric:", nrow(foss_lma), "of", nrow(foss), "\n")

# ==============================================================================
# 2. PREPARE OCCURRENCE-SPECIFIC PREDICTION ROWS
# ==============================================================================

# One row per species × site. The input data already satisfy this contract, but
# the aggregation makes the intended prediction unit explicit and robust to
# duplicated specimen rows.
foss_site <- aggregate(
  foss_lma[, "log10_pw2a_ratio", drop = FALSE],
  by  = list(fossil_name = foss_lma$fossil_name,
             species     = foss_lma$species,
             site        = foss_lma$site,
             age_ma      = foss_lma$age_ma,
             genus       = foss_lma$genus,
             family      = foss_lma$family,
             order       = foss_lma$order),
  FUN = mean, na.rm = TRUE
)

cat("Unique fossil species with petiole metric:",
    length(unique(foss_site$species)), "\n")
cat("Fossil species-site rows:", nrow(foss_site), "\n")

# ==============================================================================
# 3. LM PREDICTIONS
# ==============================================================================

lm_preds <- predict(nophy$model,
                    newdata = foss_site[, "log10_pw2a_ratio", drop = FALSE])
foss_site$lma_lm <- 10^as.numeric(lm_preds)

site_lm <- aggregate(lma_lm ~ site + age_ma, data = foss_site, FUN = mean)

cat("LM predictions done\n")

# ==============================================================================
# 4. GRAFT FOSSIL OCCURRENCES ONTO LMA TRAINING TREE
# ==============================================================================

tree                 <- read.tree("data/tre_lma_pruned.tre")
h                    <- max(nodeHeights(tree))
reference_tip_labels <- tree$tip.label
name_tbl             <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)
placement_rows       <- vector("list", nrow(foss_site))

for (j in seq_len(nrow(foss_site))) {
  result <- graft_fossil_tip(
    tree,
    tip_label = foss_site$fossil_name[j],
    age_ma = foss_site$age_ma[j],
    genus = foss_site$genus[j],
    family = foss_site$family[j],
    order = foss_site$order[j],
    name_table = name_tbl,
    reference_tip_labels = reference_tip_labels,
    tree_height = h,
    placement_fallback = PLACEMENT_FALLBACK
  )
  tree <- result$tree
  placement_rows[[j]] <- data.frame(
    fossil_name = foss_site$fossil_name[j],
    species = foss_site$species[j],
    site = foss_site$site[j],
    age_ma = foss_site$age_ma[j],
    placed = result$placed,
    placement_level = result$placement_level,
    placement_target = result$placement_target,
    age_fallback = result$age_fallback,
    error = result$error,
    stringsAsFactors = FALSE
  )
}

placement_log <- do.call(rbind, placement_rows)
write.csv(
  placement_log,
  "tables/lma_fossil_placement_log.csv",
  row.names = FALSE
)
if (!all(placement_log$placed)) {
  failed <- placement_log$fossil_name[!placement_log$placed]
  stop("LMA occurrence placement failed for: ", paste(failed, collapse = ", "))
}
cat("Placed", nrow(placement_log), "/", nrow(foss_site),
    "fossil species-site occurrences on tree\n")

# ==============================================================================
# 5. COMPUTE VCV CROSS-COVARIANCES
# ==============================================================================

idx_extant <- rownames(pip$dat_fit)   # LMA training species used in PGLS
idx_fossil <- foss_site$fossil_name

tree_small     <- keep.tip(tree, c(idx_extant, idx_fossil))
phylomat_small <- vcv(tree_small)

V_cross <- phylomat_small[idx_extant, idx_fossil, drop = FALSE] * pip$lambda
V_inv   <- solve(pip$V_lam)
resid   <- pip$resid[idx_extant]

phylo_adj <- as.numeric(t(V_cross) %*% V_inv %*% resid)
names(phylo_adj) <- idx_fossil

# ==============================================================================
# 6. PIP PREDICTIONS
# ==============================================================================

# Build the design matrix from species-within-site means; the adjustment is
# indexed by the occurrence-specific fossil tip.
placed_site <- foss_site

X_fossil <- cbind(1, placed_site$log10_pw2a_ratio)
colnames(X_fossil) <- colnames(pip$X)
rownames(X_fossil) <- placed_site$fossil_name

# GLS prediction + occurrence-specific phylogenetic correction
yhat_log10 <- as.numeric(X_fossil %*% pip$beta) +
  phylo_adj[placed_site$fossil_name]
placed_site$lma_pip <- 10^yhat_log10

site_pip <- aggregate(lma_pip ~ site + age_ma, data = placed_site, FUN = mean)

cat("PIP predictions done\n")

# ==============================================================================
# 7. COMPILE AND WRITE OUTPUT
# ==============================================================================

site_out <- Reduce(function(x, y) merge(x, y, by = c("site", "age_ma")),
                   list(site_lm, site_pip))
site_out <- site_out[order(-site_out$age_ma), ]
site_out[, c("lma_lm", "lma_pip")] <- round(site_out[, c("lma_lm", "lma_pip")], 1)

cat("\nSite-level LMA predictions (g/m²):\n")
print(site_out, row.names = FALSE)

# Per-species-site predictions
sp_out <- placed_site[, c("fossil_name", "species", "site", "age_ma", "lma_pip")]
sp_out <- merge(sp_out,
                foss_site[, c("fossil_name", "species", "site", "lma_lm")],
                by = c("fossil_name", "species", "site"))
sp_out <- sp_out[order(-sp_out$age_ma, sp_out$site, sp_out$species), ]

write.csv(site_out, "tables/lma_fossil_site_predictions.csv", row.names = FALSE)
write.csv(sp_out,   "tables/lma_fossil_species_predictions.csv", row.names = FALSE)
cat("\nSaved tables/lma_fossil_site_predictions.csv\n")
cat("Saved tables/lma_fossil_species_predictions.csv\n")
