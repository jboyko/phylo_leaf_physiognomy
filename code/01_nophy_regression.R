source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(caret)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

dat <- read.csv("data/data_species.csv")
dat$log_map <- log(dat$map)

# ==============================================================================
# 2. SELECT PREDICTORS
# ==============================================================================

# Restricted to traits measurable in fossil leaves (Dana Royer, pers. comm.).
# Training on only these traits ensures the model can be applied to fossils
# without imputing unmeasurable characters (e.g., evergreen/deciduous status).
#
# ln.leaf.area.mm2 represents Dana's leaf.area.cm2 (log-scale, same measurement).
# Tooth traits are filled with 0 (or 1 for perim.ratio) for untoothed leaves
# in 00_data_cleaning.R, so they have valid values for all species.

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

# ==============================================================================
# 2b. REPRODUCIBILITY / SAFETY HELPERS (issues #8, #14)
# ==============================================================================

SEED <- 42  # bagImpute (bagged trees) and caret's CV fold assignment are
            # stochastic; seed before every such call for bit-reproducible
            # imputed values and CV numbers (issue #8).

# aggregate(..., na.rm = TRUE) upstream (00_) can leave NaN (not NA) when a
# species/trait combination has zero observations. bagImpute's predict() is
# not guaranteed to treat NaN as missing, so coerce NaN -> NA before any
# imputation step (issue #14).
nan_to_na <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) col[is.nan(col)] <- NA
    col
  })
  df
}

# Fail loudly rather than silently propagate a corrupted (non-finite) value
# into a design matrix (issue #14).
assert_finite <- function(x, label) {
  m <- as.matrix(x)
  if (any(!is.finite(m))) {
    bad <- colnames(m)[colSums(!is.finite(m)) > 0]
    stop(sprintf("%s: non-finite values remain after imputation in columns: %s",
                 label, paste(bad, collapse = ", ")))
  }
}

predictors <- nan_to_na(predictors)

# ==============================================================================
# 3. PRE-IMPUTE ONCE ON FULL DATA — SAVED FOR FOSSIL APPLICATION ONLY
# ==============================================================================

# This full-data imputer is saved to all_results$impute_model for use on new
# fossil specimens in 04_ (no observed response there, so there is nothing to
# leak). It is NOT used to fit or score the CV loop below — that would leak
# held-out folds into the imputation model (issue #9). CV training in section
# 4 refits bagImpute inside each fold via train(preProcess = ...) instead,
# matching the leakage-free approach already used in 03_loso_cv.R.
cat("Pre-imputing predictors (full-data fit for fossil application)...\n")
set.seed(SEED)
impute_preproc <- preProcess(predictors, method = "bagImpute")
predictors_imp <- predict(impute_preproc, predictors)
assert_finite(predictors_imp, "predictors_imp (01_, full-data impute)")

# ==============================================================================
# 4. FIT SPECIES-LEVEL LM — IMPUTE AND COMPLETE-CASE VARIANTS
# ==============================================================================

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
      # Imputation is fit inside each caret CV fold via train(preProcess = ...)
      # so held-out folds cannot leak into the imputation model (issue #9).
      # predictors_imp (fit on all data, section 3) is deliberately NOT used
      # here — only for the final fossil-application model saved below.
      complete_rows  <- !is.na(dat[[target]])
      preds_fit      <- predictors
      preproc_method <- c("bagImpute", "center", "scale")
    } else {
      complete_rows  <- complete.cases(predictors) & !is.na(dat[[target]])
      preds_fit      <- predictors
      preproc_method <- c("center", "scale")
    }
    cat("  Training LM for", target, "| N:", sum(complete_rows), "\n")
    set.seed(SEED)  # seeds both per-fold bagImpute and CV fold assignment (issue #8)
    cfg_res[[target]] <- list(LM = train(
      x          = preds_fit[complete_rows, , drop = FALSE],
      y          = dat[[target]][complete_rows],
      method     = "lm",
      trControl  = ctrl,
      preProcess = preproc_method
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

# ==============================================================================
# 5. SITE-LEVEL LM — IMPUTE AND COMPLETE-CASE VARIANTS ACROSS THREE DATASETS
# ==============================================================================

# Each config specifies a site-level dataset and whether to use bagImpute or
# complete-case analysis. All configs use LM only.
#
# Backward-compatible top-level keys ($mat$LM, $log_map$LM, $impute_model,
# $pred_names) are set from specimen_impute to match the original pipeline.
# All configs are also stored under $configs for comparison.

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
  untoothed_excl_impute = list(
    file   = "data/dat_site_untoothed_excl.csv",
    impute = TRUE,
    desc   = "Peppe: species -> site (excl. untoothed), bagImpute"
  ),
  untoothed_excl_cc = list(
    file   = "data/dat_site_untoothed_excl.csv",
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

  preds <- nan_to_na(d[, pnames, drop = FALSE])

  if (cfg$impute) {
    # Full-data imputer, saved for fossil application only (e.g. 04_'s
    # site_mods$impute_model) — NOT used to fit or score the CV below, which
    # refits bagImpute inside each fold via train(preProcess = ...) instead
    # (issue #9).
    set.seed(SEED)
    imp_obj <- preProcess(preds, method = "bagImpute")
    assert_finite(predict(imp_obj, preds), paste0("site impute_model (", cfg_name, ")"))
    preproc_method <- c("bagImpute", "center", "scale")
  } else {
    imp_obj        <- NULL
    preproc_method <- c("center", "scale")
  }

  cfg_res <- list(impute_model = imp_obj, pred_names = pnames, desc = cfg$desc)

  for (target in target_vars) {
    cat("  Training LM for", target, "...\n")
    if (cfg$impute) {
      complete_rows <- !is.na(d[[target]])
    } else {
      complete_rows <- complete.cases(preds) & !is.na(d[[target]])
    }
    set.seed(SEED)  # seeds both per-fold bagImpute and CV fold assignment (issue #8)
    cfg_res[[target]] <- list(LM = train(
      x          = preds[complete_rows, , drop = FALSE],
      y          = d[[target]][complete_rows],
      method     = "lm",
      trControl  = ctrl,
      preProcess = preproc_method
    ))
  }

  site_results$configs[[cfg_name]] <- cfg_res
}

# Backward-compatible top-level keys (used by 03_ and 05_) — point to the
# original specimen -> site bagImpute config.
orig                       <- site_results$configs$specimen_impute
site_results$mat           <- orig$mat
site_results$log_map       <- orig$log_map
site_results$impute_model  <- orig$impute_model  # needed by 05_sample_size_analysis.R
site_results$pred_names    <- orig$pred_names

saveRDS(site_results, file = "models/site_models.rds")
cat("\nSaved models/site_models.rds\n")
