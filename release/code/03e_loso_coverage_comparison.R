# Compare genuinely leave-one-site-out coverage with the existing ten-fold run.
source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")
output_dir <- "tables/leave_one_site_out"
# Require every fold checkpoint before publishing combined results.
paths <- file.path("models/leave_one_site_out", sprintf("coverage_fold_%02d.rds", 1:92))
stopifnot(all(file.exists(paths)))
checkpoints <- lapply(paths, readRDS)
loo <- do.call(rbind, lapply(checkpoints, function(x) do.call(rbind, x$sites)))
sp <- do.call(rbind, lapply(checkpoints, function(x) do.call(rbind, x$species)))
old <- read.csv("tables/pip_cv_site_uncertainty.csv")

truth <- read.csv("data/dat_site.csv")
truth$log_map <- log(truth$map)
rownames(truth) <- truth$Site
stopifnot(nrow(loo) == 184L, length(unique(loo$site)) == 92L,
          !anyDuplicated(loo[c("site", "target")]),
          setequal(loo$site, old$site), length(unique(loo$fold)) == 92L)
manifest <- list()
for (fold in sort(unique(loo$fold))) {
  checkpoint <- readRDS(file.path("models/leave_one_site_out",
                                 sprintf("coverage_fold_%02d.rds", fold)))
  held <- checkpoint$held_sites
  stopifnot(checkpoint$fold == fold, checkpoint$seed == 42L + fold,
            length(held) == 1L, length(checkpoint$train_sites) == 91L,
            !held %in% checkpoint$train_sites,
            setequal(c(held, checkpoint$train_sites), unique(loo$site)),
            all(loo$site[loo$fold == fold] == held))
  for (target in c("mat", "log_map")) {
    d <- loo[loo$fold == fold & loo$target == target, ]
    species <- sp[sp$fold == fold & sp$target == target, ]
    fit <- checkpoint$model_parameters[[target]]
    degrees <- length(fit$training_species) - length(fit$beta)
    half <- qt(.975, degrees) * d$se
    stopifnot(nrow(d) == 1L, nrow(species) == d$n_species,
              abs(mean(species$estimate) - d$estimate) < 1e-8,
              abs(truth[held, target] - d$observed) < 1e-8,
              abs(d$lower - (d$estimate - half)) < 1e-8,
              abs(d$upper - (d$estimate + half)) < 1e-8)
    manifest[[paste(fold, target)]] <- data.frame(fold, held_site = held,
      n_training_sites = length(checkpoint$train_sites), target,
      n_training_species = length(fit$training_species), lambda = fit$lambda)
  }
}
write.csv(loo, file.path(output_dir, "pip_cv_site_uncertainty.csv"), row.names = FALSE)
write.csv(sp, file.path(output_dir, "pip_cv_species_uncertainty.csv"), row.names = FALSE)
write.csv(do.call(rbind, manifest), file.path(output_dir, "validation_manifest.csv"), row.names = FALSE)
summary_rows <- list()
for (scheme in c("ten_fold", "leave_one_site_out")) {
  frame <- if (scheme == "ten_fold") old else loo
  for (target in c("mat", "log_map")) {
    d <- frame[frame$target == target, ]
    covered <- d$observed >= d$lower & d$observed <= d$upper
    stopifnot(all(is.finite(d$observed)), all(is.finite(d$se)), all(d$se >= 0),
              all(is.finite(d$lower)), all(is.finite(d$upper)),
              all(covered == d$covered))
    summary_rows[[paste(scheme, target)]] <- data.frame(scheme, target,
      n_sites = nrow(d), covered = sum(covered), coverage = mean(covered),
      rmse = sqrt(mean((d$estimate - d$observed)^2)),
      mean_interval_width = mean(d$upper-d$lower),
      mean_error = mean(d$estimate-d$observed))
  }
}
comparison <- do.call(rbind, summary_rows)
write.csv(comparison, file.path(output_dir, "coverage_comparison.csv"), row.names = FALSE)
paired <- merge(old, loo, by = c("site", "target"), suffixes = c("_ten_fold", "_leave_one_site_out"))
paired$prediction_change <- paired$estimate_leave_one_site_out - paired$estimate_ten_fold
paired$se_change <- paired$se_leave_one_site_out - paired$se_ten_fold
write.csv(paired, file.path(output_dir, "paired_site_comparison.csv"), row.names = FALSE)

png("plots/pip_leave_one_site_out_coverage.png", width = 1800, height = 1000, res = 160)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
for (target in c("mat", "log_map")) {
  d <- loo[loo$target == target, ]; d <- d[order(d$observed), ]
  color <- ifelse(d$covered, "#48788b", "#bd4037")
  plot(seq_len(nrow(d)), d$estimate, ylim = range(d$lower, d$upper, d$observed),
       xlab = "Sites, ordered by observed climate", pch = 16, cex = .45, col = color,
       ylab = if (target == "mat") "MAT (degrees C)" else "ln(MAP in cm)",
       main = paste0("Leave-one-site-out: ", round(100*mean(d$covered),1), "% coverage"))
  segments(seq_len(nrow(d)), d$lower, seq_len(nrow(d)), d$upper, col = color)
  points(seq_len(nrow(d)), d$observed, pch = 16, cex = .55)
  legend("topleft", c("Observed climate", "Nominal 95% interval covers", "Nominal 95% interval misses"),
         col = c("black", "#48788b", "#bd4037"), pch = c(16, NA, NA),
         lty = c(NA,1,1), bty = "n", cex = .65)
}
dev.off()
print(comparison, row.names = FALSE)
cat("Verified: exactly one excluded site and 91 training sites per fit.\n")
