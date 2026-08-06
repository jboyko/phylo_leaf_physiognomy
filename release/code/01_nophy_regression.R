source("code/_setup.R")

library(caret)

# 1. LOAD DATA

dat <- read.csv("data/data_species.csv")
dat$log_map <- log(dat$map)

# 2. SELECT PREDICTORS

# Restricted to traits measurable on fossil leaves (D. Royer, pers. comm.), so
# the trained model applies to fossils without imputing unmeasurable characters.

fossil_traits <- c(
  # Always measurable from fossil specimens
  "pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
  # Toothed leaves only (set to 0 / 1 for untoothed in 00_)
  "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
  "avt.tooth.area", "tooth.area.blade.area.ratio",
  "tooth.area.perimeter", "tooth.area.interior", "teeth.blade.area.ratio"
)

# Drop any that are still too sparse after the upstream tooth-trait filling
NA_THRESHOLD    <- 0.40
na_pct          <- colSums(is.na(dat)) / nrow(dat)
predictor_names <- fossil_traits[fossil_traits %in% names(dat) &
                                   na_pct[fossil_traits] < NA_THRESHOLD]
dropped <- setdiff(fossil_traits, predictor_names)
if (length(dropped)) cat("Dropped (>", NA_THRESHOLD * 100, "% NA):", paste(dropped, collapse = ", "), "\n")
cat("Predictors (", length(predictor_names), "):", paste(predictor_names, collapse = ", "), "\n")

predictors <- dat[, predictor_names]

# 3. PRE-IMPUTE ONCE BEFORE CARET

# Pre-impute once rather than per CV fold: much faster, and imputation uses
# predictor correlations only (not the response), so leakage is negligible.
cat("Pre-imputing predictors...\n")
impute_preproc <- preProcess(predictors, method = "bagImpute")
predictors_imp <- predict(impute_preproc, predictors)

# 4. FIT SPECIES-LEVEL LM — IMPUTE AND COMPLETE-CASE VARIANTS

sp_configs <- list(
  impute = list(impute = TRUE,  desc = "species, bagImpute"),
  cc     = list(impute = FALSE, desc = "species, complete-case")
)

ctrl        <- trainControl(method = "cv", number = 10, savePredictions = "final")
target_vars <- c("mat", "log_map")
sp_results  <- list()

for (cfg_name in names(sp_configs)) {
  cfg <- sp_configs[[cfg_name]]
  cat("\n===", cfg$desc, "===\n")

  cfg_res <- list(desc = cfg$desc)

  for (target in target_vars) {
    if (cfg$impute) {
      complete_rows <- !is.na(dat[[target]])
      preds_fit     <- predictors_imp
    } else {
      complete_rows <- complete.cases(predictors) & !is.na(dat[[target]])
      preds_fit     <- predictors
    }
    cat("  Training LM for", target, "| N:", sum(complete_rows), "\n")
    cfg_res[[target]] <- list(LM = train(
      x          = preds_fit[complete_rows, , drop = FALSE],
      y          = dat[[target]][complete_rows],
      method     = "lm",
      trControl  = ctrl,
      preProcess = c("center", "scale")
    ))
  }

  sp_results[[cfg_name]] <- cfg_res
}

# Backward-compatible top-level keys — impute variant, used by 02_ and 05_
imp_sp      <- sp_results$impute
all_results <- list(
  configs      = sp_results,
  pred_names   = predictor_names,
  impute_model = impute_preproc,
  mat          = imp_sp$mat,
  log_map      = imp_sp$log_map
)

saveRDS(all_results, file = "models/nophy_models.rds")
cat("\nSaved models/nophy_models.rds\n")

# 5. SITE-LEVEL LM — IMPUTE AND COMPLETE-CASE VARIANTS ACROSS THREE DATASETS

# Six site-level LM configs: three aggregations x bagImpute/complete-case.
# All are stored under $configs; top-level keys mirror specimen_impute.

site_configs <- list(
  specimen_impute = list(
    file   = "data/dat_site.csv",
    impute = TRUE,
    desc   = "specimen -> site, bagImpute"
  ),
  specimen_cc = list(
    file   = "data/dat_site.csv",
    impute = FALSE,
    desc   = "specimen -> site, complete-case"
  ),
  sp_zero_impute = list(
    file   = "data/dat_site_sp_zero.csv",
    impute = TRUE,
    desc   = "species -> site (zero-fill), bagImpute"
  ),
  sp_zero_cc = list(
    file   = "data/dat_site_sp_zero.csv",
    impute = FALSE,
    desc   = "species -> site (zero-fill), complete-case"
  ),
  peppe_impute = list(
    file   = "data/dat_site_peppe.csv",
    impute = TRUE,
    desc   = "Peppe: species -> site (excl. untoothed), bagImpute"
  ),
  peppe_cc = list(
    file   = "data/dat_site_peppe.csv",
    impute = FALSE,
    desc   = "Peppe: species -> site (excl. untoothed), complete-case"
  )
)

site_results <- list(configs = list())

for (cfg_name in names(site_configs)) {
  cfg <- site_configs[[cfg_name]]
  cat("\n===", cfg$desc, "===\n")

  d           <- read.csv(cfg$file)
  d$log_map   <- log(d$map)

  na_pct_s <- colSums(is.na(d)) / nrow(d)
  pnames   <- fossil_traits[fossil_traits %in% names(d) &
                               na_pct_s[fossil_traits] < NA_THRESHOLD]
  cat("Predictors (", length(pnames), "):", paste(pnames, collapse = ", "), "\n")

  preds <- d[, pnames, drop = FALSE]

  if (cfg$impute) {
    imp_obj   <- preProcess(preds, method = "bagImpute")
    preds_fit <- predict(imp_obj, preds)
  } else {
    imp_obj   <- NULL
    preds_fit <- preds
  }

  cfg_res <- list(impute_model = imp_obj, pred_names = pnames, desc = cfg$desc)

  for (target in target_vars) {
    cat("  Training LM for", target, "...\n")
    complete_rows    <- complete.cases(preds_fit) & !is.na(d[[target]])
    cfg_res[[target]] <- list(LM = train(
      x          = preds_fit[complete_rows, , drop = FALSE],
      y          = d[[target]][complete_rows],
      method     = "lm",
      trControl  = ctrl,
      preProcess = c("center", "scale")
    ))
  }

  site_results$configs[[cfg_name]] <- cfg_res
}

# Backward-compatible top-level keys (used by 03_ and 05_) — point to the
# original specimen -> site bagImpute config.
orig                       <- site_results$configs$specimen_impute
site_results$mat           <- orig$mat
site_results$log_map       <- orig$log_map
site_results$impute_model  <- orig$impute_model
site_results$pred_names    <- orig$pred_names

saveRDS(site_results, file = "models/site_models.rds")
cat("\nSaved models/site_models.rds\n")
