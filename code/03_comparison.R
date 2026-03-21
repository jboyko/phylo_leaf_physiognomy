# ==============================================================================
# 04_comparison.R
# All model reporting: nophy CV summaries, PIP LOOCV, comparison tables, plots
#
# Requires:
#   models/pip_components.rds   -- output of 02_phy_regression.R
#   models/nophy_models.rds     -- output of 01_nophy_regression.R
# ==============================================================================

setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(caret)
library(ggplot2)
library(gridExtra)
library(pdp)

pip          <- readRDS("models/pip_components.rds")
all_results  <- readRDS("models/nophy_models.rds")
site_results <- readRDS("models/site_models.rds")

# ==============================================================================
# 1. NOPHY MODEL SUMMARIES AND VARIABLE IMPORTANCE
# ==============================================================================

for (target_name in c("mat", "log_map")) {
  comparison <- resamples(all_results[[target_name]])
  all_stats  <- summary(comparison)$statistics
  for (metric in names(all_stats)) {
    write.csv(all_stats[[metric]],
              file = paste0("tables/", target_name, "_", metric, "_summary.csv"))
  }
  imp_data <- varImp(all_results[[target_name]][["RF"]], scale = FALSE)$importance
  imp_data <- imp_data[order(imp_data[, 1], decreasing = TRUE), , drop = FALSE]
  write.csv(imp_data,
            file = paste0("tables/", target_name, "_RF_variable_importance.csv"))
  cat(target_name, "CV summary (RMSE):\n")
  print(summary(comparison)$statistics$RMSE)
}

# ==============================================================================
# 2. PARTIAL DEPENDENCE PLOTS — TOP 4 RF PREDICTORS
# ==============================================================================

export_pdp_grid <- function(target_name) {
  imp      <- varImp(all_results[[target_name]][["RF"]], scale = FALSE)
  top_vars <- rownames(imp$importance)[order(imp$importance$Overall, decreasing = TRUE)][1:4]
  p_list   <- lapply(top_vars, function(v) {
    pd_data <- partial(all_results[[target_name]][["RF"]], pred.var = v)
    ggplot(pd_data, aes(x = .data[[v]], y = .data$yhat)) +
      geom_line(linewidth = 0.8) +
      theme_minimal() +
      labs(title = v, y = "Predicted Response", x = v)
  })
  g <- grid.arrange(grobs = p_list, ncol = 2, nrow = 2,
                    top = paste("Top 4 Drivers for", target_name))
  ggsave(paste0("plots/", target_name, "_top4_PDP_grid.png"), g, width = 10, height = 8)
  cat("Saved plots/", target_name, "_top4_PDP_grid.png\n", sep = "")
}

export_pdp_grid("mat")
export_pdp_grid("log_map")

# ==============================================================================
# 3. VECTORISED PIP LOOCV
# ==============================================================================

# Standard PIP LOOCV requires solving an (N-1)x(N-1) matrix N times — O(N^4).
# Using the Schur complement identity:
#
#   V[i,-i] · (V[-i,-i])^{-1}  =  -K[i,-i] / K[ii]     where K = V^{-1}
#
# the LOO prediction simplifies to a single matrix solve + one matrix-vector
# product:
#
#   ŷ_{-i}  =  y_i  -  (K ε)_i / K_{ii}
#
# This is O(N^3) total, identical in result to the naive loop.

pip_loocv <- function(y, X, beta, V_lam) {
  K     <- solve(V_lam)
  resid <- y - as.numeric(X %*% beta)
  r     <- as.numeric(K %*% resid)
  y - r / diag(K)
}

cat("\nRunning MAT PIP LOOCV...\n")
mat_loo <- pip_loocv(
  y     = pip$dat_imputed_mat$mat,
  X     = pip$X_mat,
  beta  = pip$beta_mat,
  V_lam = pip$V_lam_mat
)

cat("Running MAP PIP LOOCV...\n")
map_loo <- pip_loocv(
  y     = pip$dat_imputed_map$log_map,
  X     = pip$X_map,
  beta  = pip$beta_map,
  V_lam = pip$V_lam_map
)

# ==============================================================================
# 3b. SITE-AGGREGATED PIP: AVERAGE SPECIES LOOCV PREDICTIONS WITHIN SITES
# ==============================================================================

# For each training site, average the species-level LOOCV predictions across
# all training species present at that site. This combines phylogenetic
# information (via PIP) with site-level averaging (noise reduction), and
# mirrors how a paleobotanist would apply PIP to a fossil assemblage.

raw_dat <- read.csv("data/RoyerLeafShapeClimateDataFixedNames_June2012.csv",
                    stringsAsFactors = FALSE)
raw_dat$genusSpecies[raw_dat$genusSpecies == " Dialyanthera sp."] <- "Dialyanthera sp."
raw_dat$genusSpecies <- gsub(" ", "_", raw_dat$genusSpecies)

train_spp <- names(mat_loo)
site_sp   <- tapply(raw_dat$genusSpecies, raw_dat$Site,
                    function(x) intersect(unique(x), train_spp))
site_sp   <- site_sp[sapply(site_sp, length) > 0]

mat_pip_site <- sapply(site_sp, function(spp) mean(mat_loo[spp]))
map_pip_site <- sapply(site_sp, function(spp) mean(map_loo[spp]))

dat_site         <- read.csv("data/dat_site.csv")
dat_site$log_map <- log(dat_site$map)
rownames(dat_site) <- dat_site$Site
dat_site_aln     <- dat_site[names(site_sp), ]

rmse_pip_site_mat <- sqrt(mean((dat_site_aln$mat     - mat_pip_site)^2))
rmse_pip_site_map <- sqrt(mean((dat_site_aln$log_map - map_pip_site)^2))
cat("Site-aggregated PIP RMSE — MAT:", round(rmse_pip_site_mat, 3),
    " log(MAP):", round(rmse_pip_site_map, 3), "\n")

# Save per-site predictions for all models with site-level output
# RF: extract CV predictions, align to dat_site row order
rf_site_mat <- site_results$mat$RF$pred[
  order(site_results$mat$RF$pred$rowIndex), ]
rf_site_map <- site_results$log_map$RF$pred[
  order(site_results$log_map$RF$pred$rowIndex), ]

site_preds <- data.frame(
  site        = dat_site$Site,
  obs_mat     = dat_site$mat,
  obs_log_map = dat_site$log_map,
  obs_map     = dat_site$map,
  rf_mat      = rf_site_mat$pred,
  rf_log_map  = rf_site_map$pred,
  pip_mat     = mat_pip_site[dat_site$Site],
  pip_log_map = map_pip_site[dat_site$Site]
)
write.csv(site_preds, "tables/site_predictions.csv", row.names = FALSE)
cat("Saved tables/site_predictions.csv\n")

# ==============================================================================
# 4. RMSE COMPARISON TABLE
# ==============================================================================

get_best_rmse <- function(model) min(model$results$RMSE)
rmse          <- function(obs, pred) sqrt(mean((obs - pred)^2))

comparison_mat <- data.frame(
  Model = c("LM (species)", "Elastic Net (species)", "Random Forest (species)",
            "PIP (species)",
            "LM (site)", "Elastic Net (site)", "Random Forest (site)",
            "PIP (site-aggregated)"),
  RMSE  = c(
    get_best_rmse(all_results$mat$LM),
    get_best_rmse(all_results$mat$ENet),
    get_best_rmse(all_results$mat$RF),
    rmse(pip$dat_imputed_mat$mat, mat_loo),
    get_best_rmse(site_results$mat$LM),
    get_best_rmse(site_results$mat$ENet),
    get_best_rmse(site_results$mat$RF),
    rmse_pip_site_mat
  ),
  Level = c(rep("Species", 3), "Phylogenetic", rep("Site", 3), "Site + Phylogenetic")
)
comparison_mat <- comparison_mat[order(comparison_mat$RMSE), ]

comparison_map <- data.frame(
  Model = c("LM (species)", "Elastic Net (species)", "Random Forest (species)",
            "PIP (species)",
            "LM (site)", "Elastic Net (site)", "Random Forest (site)",
            "PIP (site-aggregated)"),
  RMSE  = c(
    get_best_rmse(all_results$log_map$LM),
    get_best_rmse(all_results$log_map$ENet),
    get_best_rmse(all_results$log_map$RF),
    rmse(pip$dat_imputed_map$log_map, map_loo),
    get_best_rmse(site_results$log_map$LM),
    get_best_rmse(site_results$log_map$ENet),
    get_best_rmse(site_results$log_map$RF),
    rmse_pip_site_map
  ),
  Level = c(rep("Species", 3), "Phylogenetic", rep("Site", 3), "Site + Phylogenetic")
)
comparison_map <- comparison_map[order(comparison_map$RMSE), ]

cat("\nMAT model comparison (RMSE):\n");     print(comparison_mat, row.names = FALSE)
cat("\nlog(MAP) model comparison (RMSE):\n"); print(comparison_map, row.names = FALSE)

write.csv(comparison_mat, "tables/comparison_mat_rmse.csv", row.names = FALSE)
write.csv(comparison_map, "tables/comparison_map_rmse.csv", row.names = FALSE)

# ==============================================================================
# 5. SCATTER PLOTS: OBSERVED VS. PREDICTED (RF and PIP side by side)
# ==============================================================================

# RF: use stored out-of-fold CV predictions (avoids re-predicting on training data
# with mismatched column set — dat_imputed_mat only has the active PGLS traits,
# not the full fossil-trait set the RF was trained on)
rf_cv_mat <- all_results$mat$RF$pred[order(all_results$mat$RF$pred$rowIndex), ]
rf_cv_map <- all_results$log_map$RF$pred[order(all_results$log_map$RF$pred$rowIndex), ]

rf_preds_mat <- rf_cv_mat$pred
obs_mat_rf   <- rf_cv_mat$obs
rf_preds_map <- rf_cv_map$pred
obs_map_rf   <- rf_cv_map$obs

obs_mat  <- pip$dat_imputed_mat$mat
obs_map  <- pip$dat_imputed_map$log_map
lim_mat  <- range(c(obs_mat_rf, rf_preds_mat, obs_mat, mat_loo))
lim_map  <- range(c(obs_map_rf, rf_preds_map, obs_map, map_loo))

rmse_rf_mat  <- rmse(obs_mat_rf, rf_preds_mat)
rmse_pip_mat <- rmse(obs_mat, mat_loo)
rmse_rf_map  <- rmse(obs_map_rf, rf_preds_map)
rmse_pip_map <- rmse(obs_map, map_loo)

png("plots/model_comparison.png", width = 2000, height = 2000, res = 200)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(obs_mat_rf, rf_preds_mat,
     main = "Random Forest — MAT",
     xlab = "Observed MAT (°C)", ylab = "Predicted MAT (°C)",
     pch = 16, col = rgb(0.8, 0.2, 0.2, 0.3), xlim = lim_mat, ylim = lim_mat)
abline(0, 1, lwd = 2)
legend("topleft", legend = paste("RMSE:", round(rmse_rf_mat, 2)), bty = "n")

plot(obs_mat, mat_loo,
     main = "PIP — MAT",
     xlab = "Observed MAT (°C)", ylab = "Predicted MAT (°C)",
     pch = 16, col = rgb(0.2, 0.2, 0.8, 0.3), xlim = lim_mat, ylim = lim_mat)
abline(0, 1, lwd = 2)
legend("topleft", legend = paste("RMSE:", round(rmse_pip_mat, 2)), bty = "n")

plot(obs_map_rf, rf_preds_map,
     main = "Random Forest — log(MAP)",
     xlab = "Observed log(MAP)", ylab = "Predicted log(MAP)",
     pch = 16, col = rgb(0.8, 0.2, 0.2, 0.3), xlim = lim_map, ylim = lim_map)
abline(0, 1, lwd = 2)
legend("topleft", legend = paste("RMSE:", round(rmse_rf_map, 2)), bty = "n")

plot(obs_map, map_loo,
     main = "PIP — log(MAP)",
     xlab = "Observed log(MAP)", ylab = "Predicted log(MAP)",
     pch = 16, col = rgb(0.2, 0.2, 0.8, 0.3), xlim = lim_map, ylim = lim_map)
abline(0, 1, lwd = 2)
legend("topleft", legend = paste("RMSE:", round(rmse_pip_map, 2)), bty = "n")

dev.off()
cat("\nSaved plots/model_comparison.png\n")
