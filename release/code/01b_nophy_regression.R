source("code/_setup.R")

library(caret)

# 1. LOAD DATA

dat <- read.csv("data/lma_species.csv")

# 2. UNIVARIATE MODEL: log10(LMA) ~ log10(pw2.a.ratio)

# Dana Royer (pers. comm.): LMA scales linearly with PW^2/A in log-log space.
complete_rows_uni <- complete.cases(dat[, c("log10_lma", "log10_pw2a_ratio")])
cat("Univariate complete cases:", sum(complete_rows_uni), "of", nrow(dat), "species\n")

dat_fit_uni <- dat[complete_rows_uni, ]

ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")

lm_uni <- train(
  x          = dat_fit_uni[, "log10_pw2a_ratio", drop = FALSE],
  y          = dat_fit_uni$log10_lma,
  method     = "lm",
  trControl  = ctrl,
  preProcess = c("center", "scale")
)

cat("Univariate 10-fold CV RMSE:", lm_uni$results$RMSE, "\n")
cat("Univariate 10-fold CV R²:  ", lm_uni$results$Rsquared, "\n")

# 3. MULTIVARIATE MODEL: log10(LMA) ~ all 12 fossil-measurable traits

lma_fossil_traits <- c(
  "pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
  "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
  "avt.tooth.area", "tooth.area.blade.area.ratio",
  "tooth.area.perimeter", "tooth.area.interior", "teeth.blade.area.ratio"
)

# Drop predictors that are still too sparse after tooth-trait filling
NA_THRESHOLD  <- 0.40
na_pct        <- colSums(is.na(dat)) / nrow(dat)
pred_names    <- lma_fossil_traits[lma_fossil_traits %in% names(dat) &
                                     na_pct[lma_fossil_traits] < NA_THRESHOLD]
dropped       <- setdiff(lma_fossil_traits, pred_names)
if (length(dropped)) cat("Dropped (>", NA_THRESHOLD * 100, "% NA):", paste(dropped, collapse = ", "), "\n")
cat("Multivariate predictors (", length(pred_names), "):", paste(pred_names, collapse = ", "), "\n")

predictors <- dat[, pred_names, drop = FALSE]

# Pre-impute once — imputation uses predictor correlations only (not response),
# so leakage into CV estimates is minimal and far cheaper than per-fold imputation.
cat("Pre-imputing predictors...\n")
impute_preproc <- preProcess(predictors, method = "bagImpute")
predictors_imp <- predict(impute_preproc, predictors)

multi_configs <- list(
  impute = list(impute = TRUE,  desc = "multivariate, bagImpute"),
  cc     = list(impute = FALSE, desc = "multivariate, complete-case")
)

multi_results <- list()

for (cfg_name in names(multi_configs)) {
  cfg <- multi_configs[[cfg_name]]
  cat("\n===", cfg$desc, "===\n")

  if (cfg$impute) {
    complete_rows <- !is.na(dat$log10_lma)
    preds_fit     <- predictors_imp
  } else {
    complete_rows <- complete.cases(predictors) & !is.na(dat$log10_lma)
    preds_fit     <- predictors
  }

  cat("Complete cases:", sum(complete_rows), "of", nrow(dat), "species\n")

  lm_fit <- train(
    x          = preds_fit[complete_rows, , drop = FALSE],
    y          = dat$log10_lma[complete_rows],
    method     = "lm",
    trControl  = ctrl,
    preProcess = c("center", "scale")
  )

  cat("10-fold CV RMSE:", lm_fit$results$RMSE, "\n")
  cat("10-fold CV R²:  ", lm_fit$results$Rsquared, "\n")

  multi_results[[cfg_name]] <- lm_fit
}

# 4. SAVE

lma_nophy <- list(
  univariate   = list(
    model      = lm_uni,
    pred_names = "log10_pw2a_ratio",
    dat_fit    = dat_fit_uni
  ),
  multivariate = list(
    impute       = multi_results$impute,
    cc           = multi_results$cc,
    impute_model = impute_preproc,
    pred_names   = pred_names
  )
)

saveRDS(lma_nophy, file = "models/lma_nophy_models.rds")
cat("\nSaved models/lma_nophy_models.rds\n")
