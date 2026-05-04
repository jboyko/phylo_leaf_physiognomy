setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(caret)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

dat <- read.csv("data/lma_species.csv")

# ==============================================================================
# 2. FIT LM: log10(LMA) ~ log10(petiole_metric)
# ==============================================================================

# Dana Royer (pers. comm.): LMA scales linearly with PW^2/A in log-log space.
# Single predictor — no imputation or variable selection needed.

complete_rows <- complete.cases(dat[, c("log10_lma", "log10_petiole_metric")])
cat("Complete cases:", sum(complete_rows), "of", nrow(dat), "species\n")

dat_fit <- dat[complete_rows, ]

ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")

lm_lma <- train(
  x          = dat_fit[, "log10_petiole_metric", drop = FALSE],
  y          = dat_fit$log10_lma,
  method     = "lm",
  trControl  = ctrl,
  preProcess = c("center", "scale")
)

cat("10-fold CV RMSE:", lm_lma$results$RMSE, "\n")
cat("10-fold CV R²:  ", lm_lma$results$Rsquared, "\n")
cat("Coefficient (log10_petiole_metric):",
    coef(lm_lma$finalModel)[["log10_petiole_metric"]], "\n")

# ==============================================================================
# 3. SAVE
# ==============================================================================

lma_nophy <- list(
  model      = lm_lma,
  pred_names = "log10_petiole_metric",
  dat_fit    = dat_fit
)

saveRDS(lma_nophy, file = "models/lma_nophy_models.rds")
cat("Saved models/lma_nophy_models.rds\n")
