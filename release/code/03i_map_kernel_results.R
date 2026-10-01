# Audit and report the complete nested site-grouped MAP kernel experiment.
source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")
paths <- file.path("models/map_kernel_experiment",
                   sprintf("outer_%02d.rds", 1:10))
if (!all(file.exists(paths))) stop("All ten outer fold checkpoints are required")
folds <- lapply(paths, readRDS)
out <- "tables/map_kernel_experiment"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
truth <- read.csv("data/dat_site.csv", stringsAsFactors = FALSE)
rownames(truth) <- truth$Site
occurrence <- do.call(rbind, lapply(folds, `[[`, "occurrence"))
site <- do.call(rbind, lapply(folds, `[[`, "site"))
comparator <- do.call(rbind, lapply(folds, `[[`, "comparator"))
stopifnot(length(unique(site$site)) == 92L,
          setequal(site$site, truth$Site),
          !anyDuplicated(site[, c("site", "candidate")]),
          all(is.finite(site$prediction)),
          max(abs(site$observed - log(truth[site$site, "map"]))) < 1e-8,
          max(abs(site$map_cm - exp(site$prediction))) < 1e-10)
for (i in seq_along(folds)) {
  f <- folds[[i]]
  stopifnot(f$outer == i,
            !length(intersect(f$train_sites, f$held_sites)),
            setequal(c(f$train_sites, f$held_sites), truth$Site),
            setequal(f$site$site, f$held_sites),
            length(unique(f$site$candidate)) == 4L,
            length(unique(f$occurrence$candidate)) == 4L,
            all(f$selected$config_id %in% f$grid$config_id))
  members <- unlist(lapply(f$inner_membership, `[[`, "held_sites"))
  stopifnot(setequal(members, f$train_sites), !anyDuplicated(members))
  for (m in f$inner_membership) {
    stopifnot(!length(intersect(c(m$train_sites, m$held_sites),
                                f$held_sites)),
              setequal(c(m$train_sites, m$held_sites), f$train_sites))
  }
  stopifnot(max(abs(f$inner_scores$mse -
                    rowSums(f$inner_losses) / length(f$train_sites))) < 1e-12)
  for (j in seq_len(nrow(f$selected))) {
    s <- f$selected[j, ]
    group <- f$inner_scores[f$inner_scores$candidate == s$candidate, ]
    stopifnot(s$mse <= min(group$mse) + 1e-12)
  }
  for (candidate in unique(f$site$candidate)) {
    o <- f$occurrence[f$occurrence$candidate == candidate, ]
    s <- f$site[f$site$candidate == candidate, ]
    means <- tapply(o$prediction, o$site, mean)
    stopifnot(setequal(names(means), f$held_sites),
              max(abs(s$prediction - means[s$site])) < 1e-12,
              max(abs(o$prediction - o$linear_mean -
                        o$trait_adjustment - o$phylo_adjustment)) < 1e-12)
  }
  old <- readRDS(sprintf("models/map_model_experiment/outer_%02d.rds", i))
  p <- f$comparator[f$comparator$candidate == "pip", ]
  r <- f$comparator[f$comparator$candidate == "pip_recalibrated", ]
  stopifnot(max(abs(p$prediction - old$predictions$pip[
    match(p$site, old$predictions$site)])) < 1e-12,
    max(abs(r$prediction - old$predictions$pip_recalibrated[
      match(r$site, old$predictions$site)])) < 1e-12)
}
site_all <- rbind(site[, c("fold", "candidate", "site", "prediction",
                            "observed", "map_cm")], comparator)
metrics <- do.call(rbind, lapply(split(site_all, site_all$candidate), function(d) {
  e <- d$prediction - d$observed
  data.frame(candidate = d$candidate[1], n_sites = nrow(d),
             rmse = sqrt(mean(e^2)), mae = mean(abs(e)),
             mean_signed_error = mean(e),
             correlation = cor(d$prediction, d$observed),
             prediction_sd = sd(d$prediction), observed_sd = sd(d$observed),
             spread_ratio = sd(d$prediction) / sd(d$observed))
}))
stopifnot(all(metrics$n_sites == 92L))
metrics <- metrics[order(metrics$rmse), ]; rownames(metrics) <- NULL
fold_metrics <- do.call(rbind, lapply(split(site_all,
  interaction(site_all$fold, site_all$candidate, drop = TRUE)), function(d)
    data.frame(fold = d$fold[1], candidate = d$candidate[1],
               n_sites = nrow(d),
               rmse = sqrt(mean((d$prediction - d$observed)^2)))))
selection <- do.call(rbind, lapply(folds, function(f)
  cbind(fold = f$outer, f$selected)))
manifest <- do.call(rbind, lapply(folds, function(f)
  data.frame(fold = f$outer, site = c(f$train_sites, f$held_sites),
             role = c(rep("outer_train", length(f$train_sites)),
                      rep("outer_test", length(f$held_sites))),
             root_height = f$root_height,
             preprocessing_seed = f$outer_preprocessing$seed)))
stopifnot(length(unique(manifest$root_height)) == 1L)
write.csv(occurrence, file.path(out, "occurrence_predictions.csv"), row.names = FALSE)
write.csv(site_all, file.path(out, "site_predictions.csv"), row.names = FALSE)
write.csv(metrics, file.path(out, "model_comparison.csv"), row.names = FALSE)
write.csv(fold_metrics, file.path(out, "fold_metrics.csv"), row.names = FALSE)
write.csv(selection, file.path(out, "selected_parameters.csv"), row.names = FALSE)
write.csv(manifest, file.path(out, "validation_manifest.csv"), row.names = FALSE)

labels <- c(linear_kernel = "Linear kernel + phylogeny",
            rbf_kernel = "RBF kernel + phylogeny",
            linear_mean_rbf = "Linear mean + RBF + phylogeny",
            linear_mean_phylogeny = "Linear mean + phylogeny",
            pip = "Existing PIP", pip_recalibrated = "Nested recalibrated PIP")
png("plots/map_kernel_experiment_comparison.png", width = 1800,
    height = 1200, res = 160)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
limits <- range(site_all[, c("observed", "prediction")])
for (candidate in names(labels)) {
  d <- site_all[site_all$candidate == candidate, ]
  score <- metrics$rmse[metrics$candidate == candidate]
  plot(d$observed, d$prediction, xlim = limits, ylim = limits,
       pch = 16, col = "#406a84aa", cex = .85,
       xlab = "Observed ln(MAP in cm)", ylab = "Predicted ln(MAP in cm)",
       main = sprintf("%s | RMSE %.3f", labels[[candidate]], score))
  abline(0, 1, col = "#999999", lty = 2)
}
dev.off()
print(metrics, row.names = FALSE)
cat("All nested membership, comparator, aggregation, and selection checks passed.\n")
