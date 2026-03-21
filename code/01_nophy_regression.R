setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(caret)
library(glmnet)
library(ranger)

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
# 3. PRE-IMPUTE ONCE BEFORE CARET
# ==============================================================================

# bagImpute inside each caret CV fold would repeat 60 times (10 folds x 3
# methods x 2 targets). Pre-imputing once is ~60x faster with negligible
# effect on CV estimates (imputation uses predictor correlations only, not
# the response, so leakage is minimal).
cat("Pre-imputing predictors...\n")
impute_preproc <- preProcess(predictors, method = "bagImpute")
predictors_imp <- predict(impute_preproc, predictors)

# ==============================================================================
# 4. FIT MODELS VIA 10-FOLD CV
# ==============================================================================

ctrl <- trainControl(method = "cv", number = 10,
                     savePredictions = "final")

target_vars <- c("mat", "log_map")
methods     <- c(LM = "lm", ENet = "glmnet", RF = "ranger")
all_results <- list()

for (target in target_vars) {
  cat("\n--- Target:", target, "---\n")
  all_results[[target]] <- list()
  for (m_label in names(methods)) {
    cat("  Training", m_label, "...\n")
    # Use impurity importance for RF: computed during tree construction at no
    # extra cost. Permutation importance is more rigorous but much slower and
    # unnecessary for predictor ranking / PDP plots.
    tune_len <- if (methods[m_label] == "glmnet") 10 else 1
    args     <- list(x = predictors_imp, y = dat[[target]],
                     method = methods[m_label], trControl = ctrl,
                     tuneLength = tune_len, preProcess = c("center", "scale"))
    if (methods[m_label] == "ranger") args$importance <- "impurity"
    all_results[[target]][[m_label]] <- do.call(train, args)
  }
}

saveRDS(all_results, file = "models/nophy_models.rds")
cat("\nSaved models/nophy_models.rds\n")
