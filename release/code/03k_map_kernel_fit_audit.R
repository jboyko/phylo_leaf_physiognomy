# Refit selected outer models to audit jitter, rank, and saved predictions.
Sys.setenv(PIP_KERNEL_SETUP_ONLY = "1")
source("code/03h_map_kernel_cv.R")
rows <- list()
for (outer in 1:10) {
  checkpoint <- readRDS(file.path(out_dir,
    sprintf("outer_%02d.rds", outer)))
  d <- prepare_kernel_split(checkpoint$train_sites, checkpoint$held_sites,
                            42 + outer)
  for (k in seq_len(nrow(checkpoint$selected))) {
    cfg <- checkpoint$selected[k, ]
    ell <- if (is.na(cfg$ell_factor)) NA_real_ else
      cfg$ell_factor * d$reference
    fit <- map_kernel_fit(d$y, d$z, d$b, cfg$lambda_phy, cfg$eta,
                          cfg$kind, cfg$mean_kind, ell)
    pred <- map_kernel_predict(fit, d$new_z, d$cross)$prediction
    saved <- checkpoint$occurrence[
      checkpoint$occurrence$candidate == cfg$candidate, ]
    stopifnot(nrow(saved) == length(pred),
              identical(saved$site, d$occ$site),
              identical(saved$species, d$occ$species),
              max(abs(pred - saved$prediction)) < 1e-10)
    rows[[length(rows) + 1L]] <- data.frame(
      fold = outer, candidate = cfg$candidate, numerical_jitter = fit$jitter,
      retained_mean_columns = paste(fit$retained, collapse = ";"),
      removed_mean_columns = paste(fit$removed, collapse = ";"),
      removed_constant_traits = paste(d$scaling$removed, collapse = ";"),
      n_species = length(d$ids), n_traits = ncol(d$z))
  }
}
out <- do.call(rbind, rows)
dir.create("tables/map_kernel_experiment", recursive = TRUE,
           showWarnings = FALSE)
write.csv(out, "tables/map_kernel_experiment/fit_diagnostics.csv",
          row.names = FALSE)
print(table(out$numerical_jitter))
cat("All selected outer fits reproduce saved predictions.\n")
