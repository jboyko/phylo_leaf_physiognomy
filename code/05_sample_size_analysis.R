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

  # All sites with observed climate; small sites use all their specimens
  valid_sites <- names(site_sp_idx)
  cat(length(valid_sites), "sites\n")

  lm_mat_boot    <- numeric(N_BOOT)
  pip_mat_boot   <- numeric(N_BOOT)
  lm_map_boot    <- numeric(N_BOOT)
  pip_map_boot   <- numeric(N_BOOT)
  pip_ok_n_boot  <- numeric(N_BOOT)

  for (b in seq_len(N_BOOT)) {
    pip_mat_pred <- numeric(length(valid_sites))
    pip_map_pred <- numeric(length(valid_sites))
    obs_mat      <- numeric(length(valid_sites))
    obs_map      <- numeric(length(valid_sites))
    # Trait matrix for batched LM imputation + prediction
    trait_mat    <- matrix(NA_real_, nrow = length(valid_sites),
                           ncol = length(pred_names),
                           dimnames = list(valid_sites, pred_names))

    for (s_idx in seq_along(valid_sites)) {
      site_name <- valid_sites[s_idx]
      rows      <- site_sp_idx[[site_name]]
      sub_rows  <- sample(rows, min(n, length(rows)), replace = FALSE)
      sub_dat   <- raw_dat[sub_rows, ]

      # ── Collect trait means for LM (batch imputed below) ─────────────────────
      avail_traits <- intersect(pred_names, names(sub_dat))
      trait_mat[s_idx, avail_traits] <-
        colMeans(sub_dat[, avail_traits, drop = FALSE], na.rm = TRUE)

      # ── PIP prediction ────────────────────────────────────────────────────────
      # Weighted average of species-level LOOCV predictions, weighted by
      # specimen count in the subsample (matches LM's implicit abundance-weighting)
      sub_spp <- intersect(unique(sub_dat$genusSpecies), names(mat_loo))
      if (length(sub_spp) > 0) {
        spp_wts <- as.numeric(table(sub_dat$genusSpecies)[sub_spp])
        pip_mat_pred[s_idx] <- weighted.mean(mat_loo[sub_spp], w = spp_wts)
        pip_map_pred[s_idx] <- weighted.mean(map_loo[sub_spp], w = spp_wts)
      } else {
        pip_mat_pred[s_idx] <- NA
        pip_map_pred[s_idx] <- NA
      }

      obs_mat[s_idx] <- dat_site[site_name, "mat"]
      obs_map[s_idx] <- dat_site[site_name, "log_map"]
    }

    # ── Batch impute and predict for LM (one call per bootstrap) ─────────────
    trait_imp    <- predict(impute_site, newdata = as.data.frame(trait_mat))
    lm_mat_pred  <- as.numeric(predict(site_mods$mat$LM,     newdata = trait_imp))
    lm_map_pred  <- as.numeric(predict(site_mods$log_map$LM, newdata = trait_imp))

    # RMSE: restrict both models to sites where PIP has a prediction
    # (sites where at least one training species appeared in the subsample)
    pip_ok           <- !is.na(pip_mat_pred)
    pip_ok_n_boot[b] <- sum(pip_ok)
    lm_mat_boot[b]   <- sqrt(mean((obs_mat[pip_ok] - lm_mat_pred[pip_ok])^2))
    pip_mat_boot[b]  <- sqrt(mean((obs_mat[pip_ok] - pip_mat_pred[pip_ok])^2))
    lm_map_boot[b]   <- sqrt(mean((obs_map[pip_ok] - lm_map_pred[pip_ok])^2))
    pip_map_boot[b]  <- sqrt(mean((obs_map[pip_ok] - pip_map_pred[pip_ok])^2))
  }

  results[[s]] <- data.frame(
    n            = n,
    n_sites      = length(valid_sites),
    pip_n_sites  = round(mean(pip_ok_n_boot)),
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
print(ss_results[, c("n", "n_sites", "pip_n_sites", "lm_mat_mean", "pip_mat_mean",
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
