# Summarize conditional intervals; never rescale them using these same errors.
source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")
source("code/pip_uncertainty.R")
sites <- read.csv("tables/pip_cv_site_uncertainty.csv")
species <- read.csv("tables/pip_cv_species_uncertainty.csv")
stopifnot(!anyDuplicated(sites[c("site", "target")]),
          all(is.finite(sites$se)), all(sites$lower <= sites$upper))

summarize_sites <- function(d, group) {
  d <- d[is.finite(d$observed), ]
  data.frame(target = unique(d$target), group = group, n_sites = nrow(d),
    coverage_95 = mean(d$covered),
    mean_width = mean(d$upper - d$lower),
    rmse = sqrt(mean((d$estimate - d$observed)^2)),
    median_se = median(d$se),
    median_se_over_independence = median(d$se / d$independence_se))
}
summary_rows <- list()
for (target in unique(sites$target)) {
  d <- sites[sites$target == target, ]
  summary_rows[[paste(target, "all")]] <- summarize_sites(d, "all")
  for (field in c("n_species", "fraction_represented", "fraction_missing_traits", "observed")) {
    threshold <- median(d[[field]], na.rm = TRUE)
    for (side in c("below_or_equal_median", "above_median")) {
      keep <- if (side == "above_median") d[[field]] > threshold else d[[field]] <= threshold
      subset <- d[!is.na(keep) & keep, ]
      if (nrow(subset)) summary_rows[[paste(target, field, side)]] <-
        summarize_sites(subset, paste0(field, ":", side, " (", signif(threshold, 4), ")"))
    }
  }
}
summary <- do.call(rbind, summary_rows)
write.csv(summary, "tables/pip_cv_interval_coverage.csv", row.names = FALSE)

# Species intervals evaluated against local site climate are a diagnostic of
# the species-association/site-climate mismatch, not independent replicates.
species_summary <- do.call(rbind, lapply(split(species,
  list(species$target, species$represented), drop = TRUE), function(d) {
    data.frame(target = d$target[1], represented_in_training = d$represented[1],
      n_occurrences = nrow(d), n_sites = length(unique(d$site)),
      fraction_covering_local_site_climate = mean(d$covers_site_climate, na.rm = TRUE),
      mean_width = mean(d$upper - d$lower))
  }))
write.csv(species_summary, "tables/pip_cv_species_coverage.csv", row.names = FALSE)

# Calibrated intervals: shared site discrepancy plus conditional SE. tau^2 is
# estimated from the other nine folds only, so coverage is out-of-sample.
calibrated <- do.call(rbind, lapply(split(sites, sites$target), function(d) {
  d <- d[is.finite(d$observed), ]
  do.call(rbind, lapply(split(d, d$fold), function(h) {
    o <- d[d$fold != h$fold[1], ]
    tau2 <- pip_site_discrepancy_var(o$observed - o$estimate, o$se)
    cbind(h, tau2_other_folds = tau2,
          pip_calibrated_interval(h$estimate, h$se, tau2))
  }))
}))
calibrated$calibrated_covered <- calibrated$observed >= calibrated$calibrated_lower &
  calibrated$observed <= calibrated$calibrated_upper
calibrated_summary <- do.call(rbind, lapply(split(calibrated, calibrated$target), function(d) {
  full <- pip_load_site_discrepancy()
  full <- full[full$target == d$target[1], ]
  data.frame(target = d$target[1], n_sites = nrow(d),
    site_discrepancy_sd_all_folds = sqrt(full$tau2),
    rmse = full$rmse, mean_bias = mean(d$estimate - d$observed),
    conditional_coverage_95 = mean(d$covered),
    calibrated_coverage_95 = mean(d$calibrated_covered),
    calibrated_mean_width = mean(d$calibrated_upper - d$calibrated_lower))
}))
write.csv(calibrated_summary, "tables/pip_cv_calibrated_coverage.csv", row.names = FALSE)
write.csv(calibrated, "tables/pip_cv_calibrated_site_intervals.csv", row.names = FALSE)

png("plots/pip_cv_conditional_intervals.png", width = 1800, height = 1000, res = 160)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
for (target in c("mat", "log_map")) {
  d <- sites[sites$target == target, ]
  d <- d[order(d$observed), ]
  unit <- if (target == "mat") "MAT (degrees C)" else "ln(MAP in cm)"
  color <- ifelse(d$covered, "#48788b", "#bd4037")
  plot(seq_len(nrow(d)), d$estimate, ylim = range(d$lower, d$upper, d$observed),
    pch = 16, cex = 0.45, col = color, xlab = "Held-out sites, ordered by observed climate",
    ylab = unit, main = paste0(round(100 * mean(d$covered), 1), "% coverage of nominal 95% intervals"))
  segments(seq_len(nrow(d)), d$lower, seq_len(nrow(d)), d$upper, col = color)
  points(seq_len(nrow(d)), d$observed, pch = 16, cex = 0.55)
  legend("topleft", c("Observed climate", "Interval covers observation", "Interval misses observation"),
    col = c("black", "#48788b", "#bd4037"), pch = c(16, NA, NA),
    lty = c(NA, 1, 1), bty = "n", cex = 0.65)
}
dev.off()
print(summary[summary$group == "all", ], row.names = FALSE)
print(species_summary, row.names = FALSE)
print(calibrated_summary, row.names = FALSE)
