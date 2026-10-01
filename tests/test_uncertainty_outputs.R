# Validate generated intervals against the point-prediction pipeline and saved
# joint matrices, including geometric MAP aggregation and species covariance.
source("code/pip_uncertainty.R")
cv <- read.csv("tables/loso_cv_site_predictions.csv")
ci <- read.csv("tables/pip_cv_site_uncertainty.csv")
stopifnot(nrow(ci) == 184L, !anyDuplicated(ci[c("site", "target")]))
for (target in c("mat", "log_map")) {
  d <- ci[ci$target == target, ]
  point <- cv[match(d$site, cv$site), paste0("pip_sp_site_impute_", target)]
  stopifnot(max(abs(d$estimate - point)) < 1e-8,
            all(d$covered == (d$observed >= d$lower & d$observed <= d$upper)))
}

# tau^2 recomputed directly from held-out residuals, not via the helper.
tau <- do.call(rbind, lapply(c("mat", "log_map"), function(target) {
  d <- ci[ci$target == target, ]
  data.frame(target = target,
    tau2 = max(0, mean((d$observed - d$estimate)^2) - mean(d$se^2)))
}))
stopifnot(all(tau$tau2 > 0))

for (scenario in c("formal_only", "include_informal")) {
  sp <- read.csv(paste0("tables/fossil_predictions_", scenario, ".csv"))
  site <- read.csv(paste0("tables/fossil_site_uncertainty_", scenario, ".csv"))
  sp_ci <- read.csv(paste0("tables/fossil_species_uncertainty_", scenario, ".csv"))
  matrices <- readRDS(paste0("models/fossil_prediction_covariance_", scenario, ".rds"))
  stopifnot(nrow(site) == 20L, nrow(sp_ci) == 2L * nrow(sp))
  for (target in c("mat", "log_map")) {
    S <- matrices[[target]]
    stopifnot(identical(rownames(S), sp$fossil_name),
              min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
    prediction <- if (target == "mat") sp$mat_pip_site else log(sp$map_pip_site)
    d <- sp_ci[sp_ci$target == target, ]
    stopifnot(max(abs(d$se^2 - diag(S))) < 1e-8)
    for (s in unique(sp$site)) {
      ids <- which(sp$site == s)
      row <- site[site$site == s & site$target == target, ]
      stopifnot(abs(row$estimate - mean(prediction[ids])) < 1e-8,
                abs(row$se^2 - sum(S[ids, ids, drop = FALSE]) / length(ids)^2) < 1e-8)
      if (target == "log_map") stopifnot(
        abs(row$lower_response - exp(row$lower)) < 1e-8,
        abs(row$upper_response - exp(row$upper)) < 1e-8)
      # Calibrated interval: shared CV site discrepancy plus conditional SE.
      tau2 <- tau[tau$target == target, "tau2"]
      half <- qnorm(0.975) * sqrt(tau2 + row$se^2)
      back <- if (target == "mat") identity else exp
      stopifnot(abs(row$site_discrepancy_sd^2 - tau2) < 1e-8,
                abs(row$total_se^2 - (tau2 + row$se^2)) < 1e-8,
                abs(row$calibrated_lower - (row$estimate - half)) < 1e-8,
                abs(row$calibrated_upper - (row$estimate + half)) < 1e-8,
                abs(row$calibrated_lower_response - back(row$calibrated_lower)) < 1e-8,
                abs(row$calibrated_upper_response - back(row$calibrated_upper)) < 1e-8)
    }
  }
}
cat("Generated CV and fossil uncertainty checks passed.\n")
