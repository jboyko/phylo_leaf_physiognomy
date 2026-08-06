setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(caret)
library(dilp)
library(dplyr)
library(tibble)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# CONSTANTS
# ==============================================================================

K_FOLDS      <- 10
NA_THRESHOLD <- 0.40
SEED         <- 42

agg_cols <- c("pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
              "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
              "avt.tooth.area", "tooth.area.blade.area.ratio", "tooth.area.perimeter",
              "tooth.area.interior", "teeth.blade.area.ratio",
              "mat", "map", "latitude_n", "longitude_e")

fossil_traits <- c(
  "pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
  "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
  "avt.tooth.area", "tooth.area.blade.area.ratio",
  "tooth.area.perimeter", "tooth.area.interior", "teeth.blade.area.ratio"
)

vars_tooth_0 <- c(
  "primary.teeth.number", "secondary.teeth.number", "teeth.number",
  "teeth.perimeter.percm", "teeth.interior.percm", "tooth.area.cm2",
  "avt.tooth.area", "tooth.area.perimeter", "tooth.area.interior",
  "tooth.area.blade.area.ratio", "teeth.blade.area.ratio"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

fill_tooth_traits <- function(d) {
  untoothed    <- !is.na(d$margin.score) & d$margin.score == 1
  tooth_in_d   <- intersect(vars_tooth_0, names(d))
  for (v in tooth_in_d) d[[v]][untoothed] <- 0
  if ("perim.ratio" %in% names(d)) d[["perim.ratio"]][untoothed] <- 1
  d
}

# Extract coefficient table and fit summary from a plain lm object.
# method/train_level/pred_level/config/na_handling identify the model configuration.
summarise_lm <- function(fit, method, train_level, pred_level, config, na_handling, target, fold) {
  if (is.null(fit)) return(list(coefs = NULL, fit_stats = NULL))
  s   <- summary(fit)
  ct  <- as.data.frame(s$coefficients)
  names(ct) <- c("estimate", "std_error", "t_value", "p_value")
  ct$predictor <- rownames(ct)
  rownames(ct) <- NULL
  meta <- data.frame(method = method, train_level = train_level, pred_level = pred_level,
                     config = config, na_handling = na_handling, target = target, fold = fold,
                     stringsAsFactors = FALSE)
  coef_df <- cbind(meta[rep(1L, nrow(ct)), , drop = FALSE], ct, row.names = NULL)
  fit_df  <- cbind(meta, data.frame(
    r_squared     = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    residual_se   = s$sigma,
    n             = nrow(fit$model),
    df_residual   = s$df[2],
    f_stat        = unname(s$fstatistic[1]),
    f_pvalue      = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
    lambda        = NA_real_,
    stringsAsFactors = FALSE
  ), row.names = NULL)
  list(coefs = coef_df, fit_stats = fit_df)
}

# Extract coefficient table and fit summary from a pgls_res entry.
# PGLS and PIP share the same fit; this captures it once under the given method label.
# R² is computed from stored y_fit and epsilon. SEs are attempted from pgls_fit$coefficients.
summarise_pgls <- function(pc, method, config, na_handling, target, fold) {
  if (is.null(pc)) return(list(coefs = NULL, fit_stats = NULL))
  meta <- data.frame(method = method, train_level = "sp", pred_level = "site",
                     config = config, na_handling = na_handling, target = target, fold = fold,
                     stringsAsFactors = FALSE)
  ct <- tryCatch({
    raw <- coef(summary(pc$pgls_fit))
    df  <- as.data.frame(raw)
    names(df) <- c("estimate", "std_error", "t_value", "p_value")[seq_len(ncol(df))]
    if (ncol(df) < 4) for (col in c("std_error","t_value","p_value")[ncol(df):3])
      df[[col]] <- NA_real_
    df$predictor <- rownames(df)
    rownames(df) <- NULL
    df
  }, error = function(e) {
    data.frame(predictor = rownames(pc$beta), estimate = as.numeric(pc$beta),
               std_error = NA_real_, t_value = NA_real_, p_value = NA_real_,
               stringsAsFactors = FALSE)
  })
  coef_df <- cbind(meta[rep(1L, nrow(ct)), , drop = FALSE], ct, row.names = NULL)
  r2 <- tryCatch({
    y      <- pc$y_fit
    ss_res <- sum(pc$epsilon^2)
    ss_tot <- sum((y - mean(y))^2)
    1 - ss_res / ss_tot
  }, error = function(e) NA_real_)
  n      <- length(pc$epsilon)
  p      <- length(pc$beta) - 1L
  adj_r2 <- if (!is.na(r2)) 1 - (1 - r2) * (n - 1) / (n - p - 1) else NA_real_
  fit_df <- cbind(meta, data.frame(
    r_squared     = r2,
    adj_r_squared = adj_r2,
    residual_se   = sqrt(sum(pc$epsilon^2) / (n - p - 1)),
    n             = n,
    df_residual   = n - p - 1L,
    f_stat        = NA_real_,
    f_pvalue      = NA_real_,
    lambda        = pc$lambda,
    stringsAsFactors = FALSE
  ), row.names = NULL)
  list(coefs = coef_df, fit_stats = fit_df)
}

active_pred_names <- function(dat_sp, pnames) {
  present <- pnames[pnames %in% names(dat_sp)]
  na_pct  <- colSums(is.na(dat_sp[, present, drop = FALSE])) / nrow(dat_sp)
  present[na_pct < NA_THRESHOLD]
}

# Aggregate morphotype means to species-level means. Mirrors 00_data_cleaning.R section 3.
agg_species <- function(d_filled) {
  dat_sp <- aggregate(d_filled[, agg_cols],
    by  = list(genusSpecies = d_filled$genusSpecies,
               genus        = d_filled$genus,
               Family       = d_filled$Family,
               Order        = d_filled$Order),
    FUN = mean, na.rm = TRUE)
  dat_sp <- dat_sp[!is.na(dat_sp$Order) & dat_sp$Order != "unknown", ]
  # Collapse duplicate species names (same species with inconsistent Family/Order
  # across sites — a data quality issue in subsets).
  if (anyDuplicated(dat_sp$genusSpecies)) {
    dat_sp <- aggregate(dat_sp[, agg_cols],
                        by  = list(genusSpecies = dat_sp$genusSpecies),
                        FUN = mean, na.rm = TRUE)
  }
  dat_sp$log_map <- log(dat_sp$map)
  rownames(dat_sp) <- dat_sp$genusSpecies
  dat_sp
}

# Aggregate morphotype means directly to site means (tooth-filled).
agg_site_direct <- function(d_filled) {
  d_agg              <- aggregate(d_filled[, agg_cols],
                                  by  = list(Site = d_filled$Site),
                                  FUN = mean, na.rm = TRUE)
  d_agg$log_map      <- log(d_agg$map)
  rownames(d_agg)    <- d_agg$Site
  d_agg
}

# With dilp data already at morphotype (species×site) level, this collapses
# to the same result as agg_site_direct.
agg_site_sp_zero <- function(d_filled) {
  agg_site_direct(d_filled)
}

# Aggregate morphotype means → site without tooth-fill (Peppe 2011): untoothed
# species have NA tooth traits and are excluded via na.rm = TRUE.
# d_prefill is the morphotype data before tooth-trait filling.
agg_site_peppe <- function(d_prefill) {
  d_agg              <- aggregate(d_prefill[, agg_cols],
                                  by  = list(Site = d_prefill$Site),
                                  FUN = mean, na.rm = TRUE)
  d_agg$log_map      <- log(d_agg$map)
  rownames(d_agg)    <- d_agg$Site
  d_agg
}

# ==============================================================================
# 1. LOAD RAW DATA
# ==============================================================================

cat("Loading raw data via dilp...\n")
raw_leaf    <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv",
                        fileEncoding = "latin1", stringsAsFactors = FALSE)
dilp_out    <- dilp(raw_leaf)
dat_leaf    <- dilp_out$processed_leaf_data
dat_morpho  <- dilp_out$processed_morphotype_data

# Re-attach taxonomy — mirrors 00_data_cleaning.R section 1
dat_morpho <- subset(dat_morpho, select = -measurer_comments)
taxonomy_lookup <- dat_leaf %>%
  dplyr::select(site, morphotype, order, family, genus, species) %>%
  distinct()
dat_morpho <- dat_morpho %>%
  left_join(taxonomy_lookup, by = c("site", "morphotype"))

dat_morpho$genusSpecies <- gsub("[() ]", "_", paste(dat_morpho$genus, dat_morpho$species, sep = "_"))
dat_morpho$map          <- dat_morpho$map_mm / 10   # mm -> cm (pipeline convention; CLAUDE.md §8)

# Rename dilp columns to pipeline convention — mirrors 00_data_cleaning.R section 2
rename_pipeline <- function(df) {
  df %>% rename(
    pw2.a.ratio                 = petiole_metric,
    ln.leaf.area.mm2            = ln_leaf_area,
    feret.diam.ratio            = fdr,
    margin.score                = margin,
    perim.ratio                 = perimeter_ratio,
    teeth.perimeter.percm       = tc_p,
    teeth.interior.percm        = tc_ip,
    avt.tooth.area              = avg_ta,
    tooth.area.blade.area.ratio = ta_ba,
    tooth.area.perimeter        = ta_p,
    tooth.area.interior         = ta_ip,
    teeth.blade.area.ratio      = tc_ba,
    mat                         = mat_c,
    Site                        = site,
    Family                      = family,
    Order                       = order
  )
}
raw_dat <- rename_pipeline(dat_morpho)  # pre-fill; fill_tooth_traits applied per fold

dat_site_obs         <- read.csv("data/dat_site.csv", stringsAsFactors = FALSE)
dat_site_obs$log_map <- log(dat_site_obs$map)
rownames(dat_site_obs) <- dat_site_obs$Site

all_sites <- unique(raw_dat$Site)
cat("Total sites:", length(all_sites), "\n")

# ==============================================================================
# 2. PRE-COMPUTE FULL VCV (once, before the CV loop)
# ==============================================================================

cat("Pre-computing full VCV from tre_pruned.tre...\n")
full_phy              <- read.tree("data/tre_pruned.tre")
full_vcv              <- vcv(full_phy)
diag(full_vcv)        <- diag(full_vcv) + 1e-6
cat("Full VCV:", nrow(full_vcv), "x", ncol(full_vcv), "\n")

# ==============================================================================
# 3. ASSIGN CV FOLDS (stratified by site MAT, snake-pattern)
# ==============================================================================

set.seed(SEED)
site_mat_order        <- dat_site_obs[all_sites, "mat"]
names(site_mat_order) <- all_sites
mat_rank              <- rank(site_mat_order, ties.method = "first")
fold_assignment       <- ((mat_rank - 1) %% K_FOLDS) + 1
names(fold_assignment) <- all_sites

cat("Fold sizes:", table(fold_assignment), "\n")

# ==============================================================================
# 4. LOSO CV LOOP
# ==============================================================================

cv_results      <- list()
model_coefs     <- list()
model_fit_stats <- list()

for (fold in seq_len(K_FOLDS)) {

  fold_start   <- proc.time()["elapsed"]
  held_sites   <- names(fold_assignment)[fold_assignment == fold]
  train_sites  <- names(fold_assignment)[fold_assignment != fold]

  cat("\n--- Fold", fold, "/", K_FOLDS, "| held-out sites:", length(held_sites), "---\n")

  # --------------------------------------------------------------------------
  # 4a. Build training data
  # --------------------------------------------------------------------------

  train_raw    <- raw_dat[raw_dat$Site %in% train_sites, ]
  train_filled <- fill_tooth_traits(train_raw)
  train_prefill <- train_raw   # no zero-fill, for Peppe variant

  dat_sp_train <- agg_species(train_filled)

  # Restrict to species present in the full VCV
  train_sp_all <- intersect(rownames(dat_sp_train), rownames(full_vcv))
  dat_sp_train <- dat_sp_train[train_sp_all, ]

  V_train      <- full_vcv[train_sp_all, train_sp_all]

  pred_names   <- active_pred_names(dat_sp_train, fossil_traits)
  cat("  Training species:", length(train_sp_all), "| Active predictors:", length(pred_names), "\n")

  # Site-level training aggregations
  cat("  Aggregating site-level training data...\n")
  dat_site_direct_train <- agg_site_direct(train_filled)
  dat_site_sp0_train    <- agg_site_sp_zero(train_filled)
  dat_site_peppe_train  <- agg_site_peppe(train_prefill)

  # --------------------------------------------------------------------------
  # 4b. Species-level imputation (traits only)
  # --------------------------------------------------------------------------

  cat("  Fitting species-level imputation model...\n")
  imp_sp_traits <- preProcess(dat_sp_train[, pred_names, drop = FALSE],
                               method = "bagImpute")
  dat_sp_imp_X  <- predict(imp_sp_traits,
                            dat_sp_train[, pred_names, drop = FALSE])
  rownames(dat_sp_imp_X) <- rownames(dat_sp_train)

  # --------------------------------------------------------------------------
  # 4c. Fit species-level LM
  # --------------------------------------------------------------------------

  cat("  Fitting species-level LM...\n")
  lm_sp <- list()
  for (target in c("mat", "log_map")) {
    target_vec <- dat_sp_train[[target]]

    # impute variant
    complete_imp     <- !is.na(target_vec)
    df_imp           <- cbind(setNames(data.frame(target_vec[complete_imp]), target),
                               dat_sp_imp_X[complete_imp, , drop = FALSE])
    lm_sp$impute[[target]] <- tryCatch(
      lm(as.formula(paste(target, "~ .")), data = df_imp),
      error = function(e) { cat("  LM sp impute", target, "failed:", conditionMessage(e), "\n"); NULL }
    )

    # cc variant
    complete_cc      <- complete.cases(dat_sp_train[, c(target, pred_names)]) &
                          !is.na(target_vec)
    df_cc            <- dat_sp_train[complete_cc, c(target, pred_names), drop = FALSE]
    lm_sp$cc[[target]] <- tryCatch(
      lm(as.formula(paste(target, "~ .")), data = df_cc),
      error = function(e) { cat("  LM sp cc", target, "failed:", conditionMessage(e), "\n"); NULL }
    )
  }

  # --------------------------------------------------------------------------
  # 4d. Fit site-level LM (6 configs: 3 datasets x 2 impute/cc)
  # --------------------------------------------------------------------------

  site_lm_configs <- list(
    specimen_impute = list(dat = dat_site_direct_train, impute = TRUE),
    specimen_cc     = list(dat = dat_site_direct_train, impute = FALSE),
    sp_zero_impute  = list(dat = dat_site_sp0_train,    impute = TRUE),
    sp_zero_cc      = list(dat = dat_site_sp0_train,    impute = FALSE),
    peppe_impute    = list(dat = dat_site_peppe_train,  impute = TRUE),
    peppe_cc        = list(dat = dat_site_peppe_train,  impute = FALSE)
  )

  cat("  Fitting site-level LM (6 configs)...\n")
  site_lm <- list()
  for (cfg_name in names(site_lm_configs)) {
    cfg    <- site_lm_configs[[cfg_name]]
    d_site <- cfg$dat
    pnames <- active_pred_names(d_site, fossil_traits)
    preds  <- d_site[, pnames, drop = FALSE]

    if (cfg$impute) {
      imp_site  <- preProcess(preds, method = "bagImpute")
      preds_fit <- predict(imp_site, preds)
    } else {
      imp_site  <- NULL
      preds_fit <- preds
    }

    cfg_res <- list(imp_site = imp_site, pnames = pnames)
    for (target in c("mat", "log_map")) {
      complete_rows    <- complete.cases(preds_fit) & !is.na(d_site[[target]])
      df_fit           <- cbind(setNames(data.frame(d_site[[target]][complete_rows]), target),
                                 preds_fit[complete_rows, , drop = FALSE])
      cfg_res[[target]] <- tryCatch(
        lm(as.formula(paste(target, "~ .")), data = df_fit),
        error = function(e) { cat("  LM site", cfg_name, target, "failed:", conditionMessage(e), "\n"); NULL }
      )
    }
    site_lm[[cfg_name]] <- cfg_res
  }

  # --------------------------------------------------------------------------
  # 4e. Fit PGLS / PIP (impute and cc)
  # --------------------------------------------------------------------------

  cat("  Fitting PGLS/PIP (impute + cc, lambda optimisation per config)...\n")
  pgls_res <- list()
  for (pgls_cfg in c("impute", "cc")) {
    pgls_res[[pgls_cfg]] <- list()
    for (target in c("mat", "log_map")) {

      d_raw <- dat_sp_train[, c(target, pred_names), drop = FALSE]

      if (pgls_cfg == "impute") {
        imp_full      <- preProcess(d_raw, method = "bagImpute")
        d_fit         <- predict(imp_full, d_raw)
        rownames(d_fit) <- rownames(d_raw)
        sp_fit        <- rownames(d_fit)
        V_fit         <- V_train[sp_fit, sp_fit]
        # Ensure row order matches VCV
        d_fit         <- d_fit[rownames(V_fit), , drop = FALSE]
        sp_fit        <- rownames(d_fit)
      } else {
        complete_rows  <- complete.cases(d_raw)
        d_fit          <- d_raw[complete_rows, , drop = FALSE]
        sp_fit         <- rownames(d_fit)
        V_fit_raw      <- V_train[sp_fit, sp_fit]
        # Remove jitter for consistent ordering, re-apply after pruning
        diag(V_fit_raw) <- diag(V_fit_raw) - 1e-6
        V_fit          <- V_fit_raw
        diag(V_fit)    <- diag(V_fit) + 1e-6
        d_fit          <- d_fit[rownames(V_fit), , drop = FALSE]
        sp_fit         <- rownames(d_fit)
      }

      formula <- as.formula(paste(target, "~", paste(pred_names, collapse = " + ")))

      pc <- tryCatch({
        cat("    PGLS", pgls_cfg, target, "| N =", nrow(d_fit), "...\n")
        pgls_fit <- pglmEstLambda(formula = formula, data = d_fit, phylomat = V_fit)
        lambda   <- pgls_fit$lambda
        cat("    -> lambda =", round(lambda, 4), "\n")

        m    <- model.frame(formula, d_fit, na.action = na.pass)
        X    <- model.matrix(formula, m)
        beta_raw    <- coef(pgls_fit)
        beta        <- t(as.matrix(beta_raw))
        common_vars <- intersect(colnames(X), rownames(beta))
        X    <- X[, common_vars, drop = FALSE]
        beta <- beta[common_vars, , drop = FALSE]

        V_lam       <- V_fit * lambda
        diag(V_lam) <- diag(V_fit)

        epsilon      <- d_fit[[target]] - as.numeric(X %*% beta)
        names(epsilon) <- rownames(d_fit)

        K_train <- solve(V_lam)

        list(
          pgls_fit    = pgls_fit,
          beta        = beta,
          lambda      = lambda,
          K_train     = K_train,
          epsilon     = epsilon,
          y_fit       = d_fit[[target]],
          sp_fit      = sp_fit,
          common_vars = common_vars,
          formula     = formula,
          pred_names  = pred_names
        )
      }, error = function(e) {
        cat("  PGLS", pgls_cfg, target, "failed:", conditionMessage(e), "\n")
        NULL
      })

      pgls_res[[pgls_cfg]][[target]] <- pc
    }
  }

  # --------------------------------------------------------------------------
  # 4f. Extract model summaries for this fold
  # --------------------------------------------------------------------------

  cat("  Extracting model summaries...\n")

  # Species-level LM
  for (lm_variant in c("impute", "cc")) {
    for (target in c("mat", "log_map")) {
      res <- summarise_lm(lm_sp[[lm_variant]][[target]],
                          method = "LM", train_level = "sp", pred_level = "site",
                          config = "grand_mean", na_handling = lm_variant,
                          target = target, fold = fold)
      model_coefs     <- c(model_coefs,     list(res$coefs))
      model_fit_stats <- c(model_fit_stats, list(res$fit_stats))
    }
  }

  # Site-level LM (6 configs — parse aggregation method and NA handling from cfg_name)
  for (cfg_name in names(site_lm)) {
    na_handling <- if (grepl("_impute$", cfg_name)) "impute" else "cc"
    config      <- sub("_(impute|cc)$", "", cfg_name)
    for (target in c("mat", "log_map")) {
      res <- summarise_lm(site_lm[[cfg_name]][[target]],
                          method = "LM", train_level = "site", pred_level = "site",
                          config = config, na_handling = na_handling,
                          target = target, fold = fold)
      model_coefs     <- c(model_coefs,     list(res$coefs))
      model_fit_stats <- c(model_fit_stats, list(res$fit_stats))
    }
  }

  # PGLS / PIP — same fit, saved once as "PGLS" (PIP adds a prediction-time correction only)
  for (pgls_cfg in c("impute", "cc")) {
    for (target in c("mat", "log_map")) {
      res <- summarise_pgls(pgls_res[[pgls_cfg]][[target]],
                            method = "PGLS", config = "grand_mean",
                            na_handling = pgls_cfg, target = target, fold = fold)
      model_coefs     <- c(model_coefs,     list(res$coefs))
      model_fit_stats <- c(model_fit_stats, list(res$fit_stats))
    }
  }

  # --------------------------------------------------------------------------
  # 4g. Predict held-out sites
  # --------------------------------------------------------------------------

  cat("  Predicting", length(held_sites), "held-out sites...\n")
  for (s in held_sites) {
    obs_mat     <- dat_site_obs[s, "mat"]
    obs_log_map <- dat_site_obs[s, "log_map"]
    if (is.na(obs_mat) && is.na(obs_log_map)) next

    held_raw    <- raw_dat[raw_dat$Site == s, ]
    if (nrow(held_raw) == 0) next

    held_filled  <- fill_tooth_traits(held_raw)
    held_prefill <- held_raw

    rec <- list(
      site        = s,
      fold        = fold,
      obs_mat     = obs_mat,
      obs_log_map = obs_log_map
    )

    # Within-site species means
    held_sp_agg <- aggregate(held_filled[, pred_names, drop = FALSE],
                              by  = list(species = held_filled$genusSpecies),
                              FUN = mean, na.rm = TRUE)
    rownames(held_sp_agg) <- held_sp_agg$species
    held_sp_names         <- held_sp_agg$species

    # Apply species-level training imputation to held-out species
    held_sp_imp <- tryCatch(
      predict(imp_sp_traits,
              held_sp_agg[, pred_names, drop = FALSE]),
      error = function(e) NULL
    )
    if (!is.null(held_sp_imp)) rownames(held_sp_imp) <- held_sp_names

    # ---- Species-level LM predictions ----
    for (lm_variant in c("impute", "cc")) {
      for (target in c("mat", "log_map")) {
        col_name <- paste0("lm_sp_site_", lm_variant, "_", target)
        fit      <- lm_sp[[lm_variant]][[target]]
        if (is.null(fit)) { rec[[col_name]] <- NA; next }

        pred_vals <- tryCatch({
          if (lm_variant == "impute" && !is.null(held_sp_imp)) {
            as.numeric(predict(fit, newdata = held_sp_imp))
          } else {
            as.numeric(predict(fit, newdata = held_sp_agg[, pred_names, drop = FALSE]))
          }
        }, error = function(e) NULL)

        rec[[col_name]] <- if (!is.null(pred_vals)) mean(pred_vals, na.rm = TRUE) else NA
      }
    }

    # ---- Site-level LM predictions ----
    held_site_direct <- tryCatch(agg_site_direct(held_filled)[s, , drop = FALSE],
                                  error = function(e) NULL)
    held_site_sp0    <- tryCatch(agg_site_sp_zero(held_filled)[s, , drop = FALSE],
                                  error = function(e) NULL)
    held_site_peppe  <- tryCatch(agg_site_peppe(held_prefill)[s, , drop = FALSE],
                                  error = function(e) NULL)

    held_data_map <- list(
      specimen_impute = held_site_direct,
      specimen_cc     = held_site_direct,
      sp_zero_impute  = held_site_sp0,
      sp_zero_cc      = held_site_sp0,
      peppe_impute    = held_site_peppe,
      peppe_cc        = held_site_peppe
    )

    for (cfg_name in names(site_lm)) {
      cfg        <- site_lm[[cfg_name]]
      held_site  <- held_data_map[[cfg_name]]
      pnames_cfg <- cfg$pnames

      for (target in c("mat", "log_map")) {
        col_name <- paste0("lm_site_site_", cfg_name, "_", target)
        fit      <- cfg[[target]]
        if (is.null(fit) || is.null(held_site)) { rec[[col_name]] <- NA; next }

        pred_val <- tryCatch({
          avail <- intersect(pnames_cfg, names(held_site))
          row   <- held_site[, avail, drop = FALSE]
          # Apply site-level imputation if available
          if (!is.null(cfg$imp_site)) {
            row <- predict(cfg$imp_site, row)
          }
          as.numeric(predict(fit, newdata = row))
        }, error = function(e) NA)

        rec[[col_name]] <- if (length(pred_val) == 1 && !is.na(pred_val)) pred_val else NA
      }
    }

    # ---- PGLS / PIP predictions ----
    for (pgls_cfg in c("impute", "cc")) {
      for (target in c("mat", "log_map")) {
        col_pgls <- paste0("pgls_sp_site_", pgls_cfg, "_", target)
        col_pip  <- paste0("pip_sp_site_",  pgls_cfg, "_", target)

        pc <- pgls_res[[pgls_cfg]][[target]]
        if (is.null(pc)) { rec[[col_pgls]] <- NA; rec[[col_pip]] <- NA; next }

        # Build design matrix for held-out species
        X_new <- tryCatch({
          # Use imputed traits (both impute and cc variants predict using imputed inputs)
          if (is.null(held_sp_imp)) stop("no imputed held-out traits")
          newdf      <- held_sp_imp[, intersect(pc$pred_names, names(held_sp_imp)), drop = FALSE]
          # Add dummy response so model.frame works
          newdf[[target]] <- 0
          mm         <- model.matrix(pc$formula, newdf)
          avail_vars <- intersect(pc$common_vars, colnames(mm))
          mm[, avail_vars, drop = FALSE]
        }, error = function(e) NULL)

        if (is.null(X_new) || nrow(X_new) == 0) {
          rec[[col_pgls]] <- NA; rec[[col_pip]] <- NA; next
        }
        held_sp_in_X <- held_sp_names[seq_len(nrow(X_new))]

        # pc$beta is p x 1 (variable names in rownames); align X columns to beta rows
        beta_vars <- intersect(colnames(X_new), rownames(pc$beta))
        if (length(beta_vars) == 0) {
          rec[[col_pgls]] <- NA; rec[[col_pip]] <- NA; next
        }
        X_aligned    <- X_new[, beta_vars, drop = FALSE]
        beta_aligned <- pc$beta[beta_vars, , drop = FALSE]

        y_pgls <- tryCatch(
          as.numeric(X_aligned %*% beta_aligned),
          error = function(e) NULL
        )
        if (is.null(y_pgls)) { rec[[col_pgls]] <- NA; rec[[col_pip]] <- NA; next }
        names(y_pgls) <- held_sp_in_X

        rec[[col_pgls]] <- mean(y_pgls, na.rm = TRUE)

        # PIP correction: covariance between held-out species and training species
        y_pip <- y_pgls
        held_in_vcv <- intersect(held_sp_in_X, rownames(full_vcv))
        if (length(held_in_vcv) > 0) {
          tryCatch({
            C_cross  <- full_vcv[held_in_vcv, pc$sp_fit, drop = FALSE] * pc$lambda
            eps_ord  <- pc$epsilon[pc$sp_fit]
            pip_adj  <- as.numeric(C_cross %*% pc$K_train %*% eps_ord)
            names(pip_adj)           <- held_in_vcv
            y_pip[held_in_vcv]       <- y_pgls[held_in_vcv] + pip_adj
          }, error = function(e) {
            cat("  PIP correction failed for site", s, ":", conditionMessage(e), "\n")
          })
        }
        rec[[col_pip]] <- mean(y_pip, na.rm = TRUE)
      }
    }

    cv_results[[s]] <- rec
    cat("    Site:", s, "| obs MAT:", round(obs_mat, 1),
        "| pred LM-sp:", round(rec[["lm_sp_site_impute_mat"]], 1),
        "| pred PIP:", round(rec[["pip_sp_site_impute_mat"]], 1), "\n")
  }

  fold_time <- proc.time()["elapsed"] - fold_start
  cat("  Fold", fold, "complete in", round(fold_time, 1), "s\n")

  # Save per-fold model components and accumulated predictions
  fold_out <- list(
    fold          = fold,
    held_sites    = held_sites,
    train_sites   = train_sites,
    pred_names    = pred_names,
    lm_sp         = lm_sp,
    site_lm       = site_lm,
    pgls_res      = pgls_res,
    cv_results_so_far = cv_results
  )
  saveRDS(fold_out, sprintf("models/loso_cv_fold_%02d.rds", fold))
  cat(sprintf("  Saved models/loso_cv_fold_%02d.rds\n", fold))
}

# ==============================================================================
# 5. COMPILE RESULTS INTO DATA FRAME
# ==============================================================================

cat("\nCompiling results...\n")

all_keys   <- unique(unlist(lapply(cv_results, names)))
results_df <- as.data.frame(
  t(sapply(cv_results, function(rec) {
    out <- setNames(rep(NA_real_, length(all_keys)), all_keys)
    for (k in names(rec)) if (length(rec[[k]]) == 1) out[[k]] <- as.numeric(rec[[k]])
    out
  })),
  stringsAsFactors = FALSE
)
results_df$site <- sapply(cv_results, `[[`, "site")

cat("Results:", nrow(results_df), "sites\n")

# Column naming convention: {method}_{train_level}_{pred_level}_{config}_{na}_{target}
#   train_level: sp  = trained on species-level means
#                site = trained on site-level means
#   pred_level:  site = prediction aggregated to site (always, in this LOSO script)
#   config:      impute/cc for sp-trained models; specimen/sp_zero/peppe + impute/cc for site LMs
# Examples: lm_sp_site_impute_mat, lm_site_site_peppe_cc_log_map, pip_sp_site_cc_mat
# Names are set directly in the loop so fold RDS files also carry the correct names.

# ==============================================================================
# 6. RMSE TABLE
# ==============================================================================

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

pred_cols <- grep("_(mat|log_map)$", names(results_df), value = TRUE)
pred_cols <- setdiff(pred_cols, c("obs_mat", "obs_log_map"))

rmse_rows <- lapply(pred_cols, function(col) {
  target   <- if (grepl("_mat$", col)) "mat" else "log_map"
  obs      <- results_df[[paste0("obs_", target)]]
  pred     <- results_df[[col]]
  complete <- !is.na(obs) & !is.na(pred)
  if (sum(complete) == 0) return(NULL)
  data.frame(
    column   = col,
    target   = target,
    rmse     = rmse(obs[complete], pred[complete]),
    n_sites  = sum(complete),
    stringsAsFactors = FALSE
  )
})
rmse_rows      <- do.call(rbind, Filter(Negate(is.null), rmse_rows))
rmse_rows      <- rmse_rows[order(rmse_rows$target, rmse_rows$rmse), ]
rownames(rmse_rows) <- NULL

cat("\nLOSO CV RMSE summary:\n")
print(rmse_rows, row.names = FALSE)

# ==============================================================================
# 7. WRITE OUTPUTS
# ==============================================================================

model_coefs_df     <- do.call(rbind, Filter(Negate(is.null), model_coefs))
model_fit_stats_df <- do.call(rbind, Filter(Negate(is.null), model_fit_stats))
rownames(model_coefs_df)     <- NULL
rownames(model_fit_stats_df) <- NULL

write.csv(rmse_rows,          "tables/loso_cv_rmse.csv",             row.names = FALSE)
write.csv(results_df,         "tables/loso_cv_site_predictions.csv", row.names = FALSE)
write.csv(model_coefs_df,     "tables/loso_cv_model_coefs.csv",      row.names = FALSE)
write.csv(model_fit_stats_df, "tables/loso_cv_model_fit.csv",        row.names = FALSE)

cat("\nSaved tables/loso_cv_rmse.csv\n")
cat("Saved tables/loso_cv_site_predictions.csv\n")
cat("Saved tables/loso_cv_model_coefs.csv\n")
cat("Saved tables/loso_cv_model_fit.csv\n")
