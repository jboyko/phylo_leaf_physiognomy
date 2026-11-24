# Compare held-out site errors with the predictions and fold assignments fixed.
source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

paired_rmse <- function(data, observed, pip, baseline, target, baseline_name) {
  scores <- data[, c(observed, pip, baseline)]
  paired <- scores[Reduce(`&`, lapply(scores, is.finite)), ]
  n_sites <- nrow(paired)
  stopifnot(!anyNA(data$site), !anyDuplicated(data$site), n_sites > 1L)

  pip_error <- (paired[[pip]] - paired[[observed]])^2
  baseline_error <- (paired[[baseline]] - paired[[observed]])^2
  pip_rmse <- sqrt(mean(pip_error))
  baseline_rmse <- sqrt(mean(baseline_error))

  # Resample sites and use the same draw for both models in each replicate.
  # This does not account for dependence due to overlapping training folds.
  set.seed(42L)
  differences <- replicate(10000L, {
    sites <- sample.int(n_sites, replace = TRUE)
    sqrt(mean(pip_error[sites])) - sqrt(mean(baseline_error[sites]))
  })
  interval <- quantile(differences, c(0.025, 0.975), names = FALSE)

  data.frame(
    target = target, baseline = baseline_name, n_sites = n_sites,
    pip_rmse = pip_rmse, baseline_rmse = baseline_rmse,
    difference = pip_rmse - baseline_rmse,
    lower_95 = interval[1], upper_95 = interval[2]
  )
}

climate <- read.csv("tables/loso_cv_site_predictions.csv")
lma <- read.csv("tables/lma_cv_site_predictions.csv")

climate_baselines <- c(
  "Refitted DiLP" = "dilp_cv_site",
  "12-trait site LM" = "lm_site_site_sp_zero_impute",
  "Untoothed-excluded site LM" = "lm_site_site_untoothed_excl_impute",
  "PGLS" = "pgls_sp_site_impute"
)

climate_results <- lapply(c("mat", "log_map"), function(target) {
  do.call(rbind, lapply(names(climate_baselines), function(baseline_name) {
    paired_rmse(
      climate, paste0("obs_", target), paste0("pip_sp_site_impute_", target),
      paste0(climate_baselines[[baseline_name]], "_", target), target, baseline_name
    )
  }))
})

results <- rbind(
  do.call(rbind, climate_results),
  paired_rmse(lma, "obs", "pred_uni_pip", "pred_uni_lm", "log10_lma", "Univariate LM"),
  paired_rmse(lma, "obs", "pred_uni_pip", "pred_uni_pgls", "log10_lma", "Univariate PGLS")
)
write.csv(results, "tables/cv_paired_rmse.csv", row.names = FALSE)
print(results, row.names = FALSE)
