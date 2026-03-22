setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# USER SETTINGS
# ==============================================================================

SAMPLE_SIZES <- c(1, 2, 3, 5, 8, 12, 20, 30, 50)
N_BOOT       <- 200   # bootstrap replicates per sample size
set.seed(42)

# ==============================================================================
# 1. LOAD MODELS AND DATA
# ==============================================================================

pip       <- readRDS("models/pip_components.rds")
site_mods <- readRDS("models/site_models.rds")

raw_dat <- read.csv("data/RoyerLeafShapeClimateDataFixedNames_June2012.csv",
                    stringsAsFactors = FALSE)
raw_dat$genusSpecies[raw_dat$genusSpecies == " Dialyanthera sp."] <- "Dialyanthera sp."
raw_dat$genusSpecies <- gsub(" ", "_", raw_dat$genusSpecies)

dat_site         <- read.csv("data/dat_site.csv")
dat_site$log_map <- log(dat_site$map)
rownames(dat_site) <- dat_site$Site

# ==============================================================================
# 2. RECOMPUTE SPECIES-LEVEL PIP LOOCV PREDICTIONS
# ==============================================================================

pip_loocv <- function(y, X, beta, V_lam) {
  K     <- solve(V_lam)
  resid <- y - as.numeric(X %*% beta)
  y - as.numeric(K %*% resid) / diag(K)
}

mat_loo <- pip_loocv(pip$dat_imputed_mat$mat,
                     pip$X_mat, pip$beta_mat, pip$V_lam_mat)
map_loo <- pip_loocv(pip$dat_imputed_map$log_map,
                     pip$X_map, pip$beta_map, pip$V_lam_map)
names(mat_loo) <- rownames(pip$dat_imputed_mat)
names(map_loo) <- rownames(pip$dat_imputed_map)

# ==============================================================================
# 3. PREPARE SITE-SPECIMEN INDEX
# ==============================================================================

pred_names   <- site_mods$pred_names
impute_site  <- site_mods$impute_model

# Only keep raw specimens whose site has an observed climate value
raw_dat      <- raw_dat[raw_dat$Site %in% rownames(dat_site), ]
site_sp_idx  <- split(seq_len(nrow(raw_dat)), raw_dat$Site)

# ==============================================================================
# 4. BOOTSTRAP RMSE ACROSS SAMPLE SIZES
# ==============================================================================

results <- vector("list", length(SAMPLE_SIZES))

for (s in seq_along(SAMPLE_SIZES)) {
  n <- SAMPLE_SIZES[s]
  cat("Sample size N =", n, "... ")

  # Sites with at least n specimens and observed climate
  valid_sites <- names(site_sp_idx)[sapply(site_sp_idx, length) >= n]
  cat(length(valid_sites), "sites\n")

  lm_mat_boot  <- numeric(N_BOOT)
  pip_mat_boot <- numeric(N_BOOT)
  lm_map_boot  <- numeric(N_BOOT)
  pip_map_boot <- numeric(N_BOOT)

  for (b in seq_len(N_BOOT)) {
    lm_mat_pred  <- numeric(length(valid_sites))
    pip_mat_pred <- numeric(length(valid_sites))
    lm_map_pred  <- numeric(length(valid_sites))
    pip_map_pred <- numeric(length(valid_sites))
    obs_mat      <- numeric(length(valid_sites))
    obs_map      <- numeric(length(valid_sites))

    for (s_idx in seq_along(valid_sites)) {
      site_name <- valid_sites[s_idx]
      rows      <- site_sp_idx[[site_name]]
      sub_rows  <- sample(rows, n, replace = FALSE)
      sub_dat   <- raw_dat[sub_rows, ]

      # ── LM prediction ────────────────────────────────────────────────────────
      # Compute subsample-mean traits, impute NAs, predict with trained LM
      avail_traits <- intersect(pred_names, names(sub_dat))
      trait_means  <- colMeans(sub_dat[, avail_traits, drop = FALSE], na.rm = TRUE)
      trait_df     <- as.data.frame(as.list(trait_means))
      # Pad any missing columns with NA (will be imputed)
      for (tr in pred_names[!pred_names %in% avail_traits]) trait_df[[tr]] <- NA
      trait_df  <- trait_df[, pred_names, drop = FALSE]
      trait_imp <- predict(impute_site, newdata = trait_df)

      lm_mat_pred[s_idx] <- predict(site_mods$mat$LM,     newdata = trait_imp)
      lm_map_pred[s_idx] <- predict(site_mods$log_map$LM, newdata = trait_imp)

      # ── PIP prediction ────────────────────────────────────────────────────────
      # Average species-level LOOCV predictions for species in the subsample
      sub_spp <- intersect(unique(sub_dat$genusSpecies), names(mat_loo))
      pip_mat_pred[s_idx] <- if (length(sub_spp) > 0) mean(mat_loo[sub_spp]) else NA
      pip_map_pred[s_idx] <- if (length(sub_spp) > 0) mean(map_loo[sub_spp]) else NA

      obs_mat[s_idx] <- dat_site[site_name, "mat"]
      obs_map[s_idx] <- dat_site[site_name, "log_map"]
    }

    # RMSE: for PIP drop sites where no training species appeared in the subsample
    pip_ok <- !is.na(pip_mat_pred)
    lm_mat_boot[b]  <- sqrt(mean((obs_mat        - lm_mat_pred)^2))
    pip_mat_boot[b] <- sqrt(mean((obs_mat[pip_ok] - pip_mat_pred[pip_ok])^2))
    lm_map_boot[b]  <- sqrt(mean((obs_map        - lm_map_pred)^2))
    pip_map_boot[b] <- sqrt(mean((obs_map[pip_ok] - pip_map_pred[pip_ok])^2))
  }

  results[[s]] <- data.frame(
    n            = n,
    n_sites      = length(valid_sites),
    lm_mat_mean  = mean(lm_mat_boot),
    lm_mat_lo    = quantile(lm_mat_boot,  0.025),
    lm_mat_hi    = quantile(lm_mat_boot,  0.975),
    pip_mat_mean = mean(pip_mat_boot),
    pip_mat_lo   = quantile(pip_mat_boot, 0.025),
    pip_mat_hi   = quantile(pip_mat_boot, 0.975),
    lm_map_mean  = mean(lm_map_boot),
    lm_map_lo    = quantile(lm_map_boot,  0.025),
    lm_map_hi    = quantile(lm_map_boot,  0.975),
    pip_map_mean = mean(pip_map_boot),
    pip_map_lo   = quantile(pip_map_boot, 0.025),
    pip_map_hi   = quantile(pip_map_boot, 0.975),
    row.names    = NULL
  )
}

ss_results <- do.call(rbind, results)
write.csv(ss_results, "tables/sample_size_results.csv", row.names = FALSE)
cat("\nSaved tables/sample_size_results.csv\n")
print(ss_results[, c("n", "n_sites", "lm_mat_mean", "pip_mat_mean",
                      "lm_map_mean", "pip_map_mean")])

# ==============================================================================
# 5. SAVE PER-SITE SPECIMEN COUNTS (used to position full-site CV reference)
# ==============================================================================

site_counts <- sapply(site_sp_idx, length)
write.csv(
  data.frame(median_n = median(site_counts), mean_n = mean(site_counts)),
  "tables/site_specimen_counts.csv", row.names = FALSE
)
cat("Median specimens per site:", median(site_counts), "\n")
cat("Mean specimens per site:  ", round(mean(site_counts), 1), "\n")
