# Add the refitted DiLP benchmark to an existing site-grouped CV run.
# PIP and 12-trait LM predictions are reused, not refitted or resplit.
source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")
source("code/dilp_cv.R")
source("code/site_grouping.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  # Explicit local package checkout for environments without installed dilp.
  devtools::load_all(args[1], quiet = TRUE)
} else {
  pip_require_dilp()
}
raw <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv",
                fileEncoding = "latin1")
raw$site <- normalise_calibration_site(raw$site)
processed <- dilp(raw)$processed_site_data
cv <- read.csv("tables/loso_cv_site_predictions.csv", stringsAsFactors = FALSE)
result <- dilp_site_cv(processed, cv)
write.csv(result$predictions, "tables/dilp_cv_site_predictions.csv", row.names = FALSE)
write.csv(result$coefficients, "tables/dilp_cv_coefficients.csv", row.names = FALSE)

# Use the exact same eligible sites for every model within each target.
models <- c("PIP" = "pip_sp_site_impute",
            "12-trait site regression" = "lm_site_site_sp_zero_impute",
            "DiLP regression" = "dilp_cv_site")
scores <- list()
for (target in c("mat", "log_map")) {
  obs <- cv[[paste0("obs_", target)]]
  preds <- cbind(cv[, paste0(models[1:2], "_", target)],
                 result$predictions[, paste0(models[3], "_", target)])
  names(preds) <- names(models)
  common <- is.finite(obs) & apply(preds, 1, function(x) all(is.finite(x)))
  for (model in names(models)) {
    scores[[paste(model, target)]] <- data.frame(
      model = model, target = target, n_sites = sum(common),
      rmse = sqrt(mean((preds[[model]][common] - obs[common])^2)),
      evaluation = "10-fold site-grouped CV, common sites")
  }
  excluded <- cv$site[!common]
  cat(target, "excluded from common-site comparison:", paste(excluded, collapse = ", "), "\n")
}
scores <- do.call(rbind, scores)
write.csv(scores, "tables/dana_cv_comparison.csv", row.names = FALSE)
print(scores, row.names = FALSE)
