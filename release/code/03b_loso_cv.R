source("code/_setup.R")

library(ape)
library(dplyr)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# CONSTANTS

K_FOLDS <- 10
SEED    <- 42
target  <- "log10_lma"

uni_pred      <- "log10_pw2a_ratio"
uni_formula   <- as.formula(paste(target, "~", uni_pred))

lma_fossil_traits <- c(
  "pw2.a.ratio", "ln.leaf.area.mm2", "feret.diam.ratio", "margin.score",
  "perim.ratio", "teeth.perimeter.percm", "teeth.interior.percm",
  "avt.tooth.area", "tooth.area.blade.area.ratio",
  "tooth.area.perimeter", "tooth.area.interior", "teeth.blade.area.ratio"
)
multi_formula <- as.formula(paste(target, "~", paste(lma_fossil_traits, collapse = " + ")))

# 1. LOAD AND PROCESS RAW DATA  (mirrors 00b_data_cleaning.R)

lma_cols <- c("specimen_number", "site", "morphotype", "order", "family", "genus", "species",
              "margin", "lma_g_m2",
              "petiole_width", "petiole_area", "blade_area", "blade_perimeter",
              "feret", "minimum_feret", "raw_blade_area", "raw_blade_perimeter",
              "internal_raw_blade_area", "internal_raw_blade_perimeter",
              "length_of_cut_perimeter", "no_of_primary_teeth", "no_of_subsidiary_teeth",
              "MAT_C", "MAP_mm", "latitude_N", "longitude_E")

input_1 <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv",
                    fileEncoding = "latin1", stringsAsFactors = FALSE)
input_2 <- read.csv("data/extra_calibration_data_for_LMA.csv",
                    fileEncoding = "latin1", stringsAsFactors = FALSE)

input <- bind_rows(
  input_1 %>% dplyr::select(any_of(lma_cols)),
  input_2 %>% dplyr::select(any_of(lma_cols))
)

dilp_out    <- dilp(input)
leaf_data   <- dilp_out$processed_leaf_data
morpho_data <- dilp_out$processed_morphotype_data

taxonomy_lookup <- leaf_data %>%
  dplyr::select(site, morphotype, order, family, genus, species) %>%
  distinct()

morpho_data <- morpho_data %>%
  dplyr::select(-any_of("measurer_comments")) %>%
  left_join(taxonomy_lookup, by = c("site", "morphotype")) %>%
  filter(!is.na(genus), genus != "", !is.na(species), species != "")

# Tooth-trait filling for confirmed untoothed leaves (mirrors 00b_data_cleaning.R)
morpho_data <- morpho_data %>%
  mutate(
    tc_p            = ifelse(margin == 1, 0, tc_p),
    tc_ip           = ifelse(margin == 1, 0, tc_ip),
    tc_ba           = ifelse(margin == 1, 0, tc_ba),
    ta_p            = ifelse(margin == 1, 0, ta_p),
    ta_ip           = ifelse(margin == 1, 0, ta_ip),
    ta_ba           = ifelse(margin == 1, 0, ta_ba),
    avg_ta          = ifelse(margin == 1, 0, avg_ta),
    perimeter_ratio = ifelse(margin == 1, 1, perimeter_ratio)
  )

# Rename dilp columns to pipeline convention (mirrors 00b_data_cleaning.R)
morpho_data <- morpho_data %>%
  rename(
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
    teeth.blade.area.ratio      = tc_ba
  ) %>%
  mutate(
    genusSpecies     = paste(genus, species, sep = "_"),
    log10_lma        = log10(lma_g_m2),
    log10_pw2a_ratio = log10(pw2.a.ratio)
  ) %>%
  filter(!is.na(log10_lma))

cat("Complete morphotype rows:", nrow(morpho_data), "\n")
cat("Sites:", length(unique(morpho_data$site)), "\n")

# 2. SITE-LEVEL OBSERVED LMA

site_obs <- morpho_data %>%
  group_by(site, genusSpecies) %>%
  summarise(log10_lma = mean(log10_lma, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  summarise(obs_log10_lma = mean(log10_lma, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

rownames(site_obs) <- site_obs$site

# 3. PRE-COMPUTE FULL VCV

phy      <- read.tree("data/tre_lma_pruned.tre")
full_vcv <- vcv(phy)
diag(full_vcv) <- diag(full_vcv) + 1e-6
cat("Full VCV:", nrow(full_vcv), "x", ncol(full_vcv), "\n")

# 4. ASSIGN LOSO FOLDS (stratified by site-mean log10_lma, snake pattern)

set.seed(SEED)
foldable_sites <- site_obs$site[!is.na(site_obs$obs_log10_lma)]
lma_rank       <- rank(site_obs[foldable_sites, "obs_log10_lma"], ties.method = "first")
fold_assignment <- ((lma_rank - 1) %% K_FOLDS) + 1
names(fold_assignment) <- foldable_sites

cat("Fold sizes:", paste(table(fold_assignment), collapse = " "), "\n")

# 5. LOSO CV LOOP

cv_results      <- list()
model_coefs     <- list()
model_fit_stats <- list()

for (fold in seq_len(K_FOLDS)) {

  fold_start  <- proc.time()["elapsed"]
  held_sites  <- names(fold_assignment)[fold_assignment == fold]
  train_sites <- names(fold_assignment)[fold_assignment != fold]

  cat("\n--- Fold", fold, "/", K_FOLDS,
      "| held-out sites:", length(held_sites), "---\n")

  # --------------------------------------------------------------------------
  # 5a. Aggregate training morphotypes to species means
  # --------------------------------------------------------------------------

  train_morpho <- morpho_data[morpho_data$site %in% train_sites, ]

  all_pred_cols <- c(target, uni_pred, lma_fossil_traits)
  agg_cols      <- intersect(all_pred_cols, names(train_morpho))

  dat_sp_train <- train_morpho %>%
    group_by(genusSpecies) %>%
    summarise(across(all_of(agg_cols), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    filter(!is.na(.data[[target]])) %>%
    as.data.frame()

  if (anyDuplicated(dat_sp_train$genusSpecies)) {
    dat_sp_train <- aggregate(
      dat_sp_train[, agg_cols],
      by = list(genusSpecies = dat_sp_train$genusSpecies),
      FUN = mean, na.rm = TRUE
    )
  }

  rownames(dat_sp_train) <- dat_sp_train$genusSpecies

  # Restrict to species present in the full VCV
  train_sp     <- intersect(rownames(dat_sp_train), rownames(full_vcv))
  dat_sp_train <- dat_sp_train[train_sp, ]
  V_train      <- full_vcv[train_sp, train_sp]

  cat("  Training species:", length(train_sp), "\n")

  # --------------------------------------------------------------------------
  # 5b. Multivariate imputation (fit on training only to avoid leakage)
  # --------------------------------------------------------------------------

  multi_preds_train <- dat_sp_train[, lma_fossil_traits, drop = FALSE]
  multi_imp_obj     <- tryCatch(
    preProcess(multi_preds_train, method = "bagImpute"),
    error = function(e) { cat("  Imputation failed:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(multi_imp_obj)) {
    multi_preds_imp <- predict(multi_imp_obj, multi_preds_train)
    rownames(multi_preds_imp) <- rownames(multi_preds_train)
  }

  # --------------------------------------------------------------------------
  # 5c. Fit univariate LM (complete-case)
  # --------------------------------------------------------------------------

  dat_uni_cc <- dat_sp_train[complete.cases(dat_sp_train[, c(target, uni_pred)]), ]

  lm_uni <- tryCatch(
    lm(uni_formula, data = dat_uni_cc),
    error = function(e) { cat("  Uni LM failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(lm_uni)) {
    s <- summary(lm_uni)
    ct <- as.data.frame(s$coefficients)
    names(ct) <- c("estimate", "std_error", "t_value", "p_value")
    ct$predictor <- rownames(ct); rownames(ct) <- NULL
    meta <- data.frame(method = "LM_uni", fold = fold, target = target, stringsAsFactors = FALSE)
    model_coefs     <- c(model_coefs,     list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
    model_fit_stats <- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
      residual_se = s$sigma, n = nrow(lm_uni$model), lambda = NA_real_
    ), row.names = NULL)))
  }

  # --------------------------------------------------------------------------
  # 5d. Fit multivariate LM — impute variant
  # --------------------------------------------------------------------------

  lm_multi_imp <- NULL
  if (!is.null(multi_imp_obj)) {
    dat_multi_imp        <- dat_sp_train
    dat_multi_imp[, lma_fossil_traits] <- multi_preds_imp
    complete_multi_imp   <- complete.cases(dat_multi_imp[, c(target, lma_fossil_traits)])

    lm_multi_imp <- tryCatch(
      lm(multi_formula, data = dat_multi_imp[complete_multi_imp, ]),
      error = function(e) { cat("  Multi LM impute failed:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(lm_multi_imp)) {
      s <- summary(lm_multi_imp)
      ct <- as.data.frame(s$coefficients)
      names(ct) <- c("estimate", "std_error", "t_value", "p_value")
      ct$predictor <- rownames(ct); rownames(ct) <- NULL
      meta <- data.frame(method = "LM_multi_imp", fold = fold, target = target, stringsAsFactors = FALSE)
      model_coefs     <- c(model_coefs,     list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
      model_fit_stats <- c(model_fit_stats, list(cbind(meta, data.frame(
        r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
        residual_se = s$sigma, n = nrow(lm_multi_imp$model), lambda = NA_real_
      ), row.names = NULL)))
    }
  }

  # --------------------------------------------------------------------------
  # 5e. Fit multivariate LM — complete-case variant
  # --------------------------------------------------------------------------

  complete_multi_cc <- complete.cases(dat_sp_train[, c(target, lma_fossil_traits)])
  lm_multi_cc <- tryCatch(
    lm(multi_formula, data = dat_sp_train[complete_multi_cc, ]),
    error = function(e) { cat("  Multi LM cc failed:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(lm_multi_cc)) {
    s <- summary(lm_multi_cc)
    ct <- as.data.frame(s$coefficients)
    names(ct) <- c("estimate", "std_error", "t_value", "p_value")
    ct$predictor <- rownames(ct); rownames(ct) <- NULL
    meta <- data.frame(method = "LM_multi_cc", fold = fold, target = target, stringsAsFactors = FALSE)
    model_coefs     <- c(model_coefs,     list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
    model_fit_stats <- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
      residual_se = s$sigma, n = nrow(lm_multi_cc$model), lambda = NA_real_
    ), row.names = NULL)))
  }

  # --------------------------------------------------------------------------
  # 5f. Fit univariate PGLS / PIP
  # --------------------------------------------------------------------------

  pc_uni <- tryCatch({
    sp_uni     <- intersect(rownames(dat_uni_cc), train_sp)
    d_uni_phy  <- dat_uni_cc[sp_uni, ]
    V_uni      <- V_train[sp_uni, sp_uni]

    pgls_uni   <- pglmEstLambda(formula = uni_formula, data = d_uni_phy, phylomat = V_uni)
    lambda_uni <- pgls_uni$lambda
    cat("  Uni lambda:", round(lambda_uni, 4), "\n")

    m_uni    <- model.frame(uni_formula, d_uni_phy, na.action = na.pass)
    X_uni    <- model.matrix(uni_formula, m_uni)
    beta_uni <- t(as.matrix(coef(pgls_uni)))
    cv_uni   <- intersect(colnames(X_uni), rownames(beta_uni))
    X_uni    <- X_uni[, cv_uni, drop = FALSE]
    beta_uni <- beta_uni[cv_uni, , drop = FALSE]

    V_lam_uni    <- V_uni * lambda_uni
    diag(V_lam_uni) <- diag(V_uni)
    eps_uni      <- d_uni_phy[[target]] - as.numeric(X_uni %*% beta_uni)
    names(eps_uni) <- sp_uni
    K_uni        <- solve(V_lam_uni)

    n <- length(eps_uni); p <- length(cv_uni) - 1L
    r2 <- 1 - sum(eps_uni^2) / sum((d_uni_phy[[target]] - mean(d_uni_phy[[target]]))^2)
    ct <- tryCatch({
      df <- as.data.frame(coef(summary(pgls_uni)))
      names(df) <- c("estimate", "std_error", "t_value", "p_value")[seq_len(ncol(df))]
      df$predictor <- rownames(df); rownames(df) <- NULL; df
    }, error = function(e) {
      data.frame(predictor = rownames(beta_uni), estimate = as.numeric(beta_uni),
                 std_error = NA_real_, t_value = NA_real_, p_value = NA_real_)
    })
    meta <- data.frame(method = "PGLS_uni", fold = fold, target = target, stringsAsFactors = FALSE)
    model_coefs     <<- c(model_coefs,     list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
    model_fit_stats <<- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = r2, adj_r_squared = 1 - (1 - r2) * (n - 1) / (n - p - 1),
      residual_se = sqrt(sum(eps_uni^2) / (n - p - 1)), n = n, lambda = lambda_uni
    ), row.names = NULL)))

    list(beta = beta_uni, lambda = lambda_uni, K_train = K_uni,
         epsilon = eps_uni, common_vars = cv_uni, sp = sp_uni)
  }, error = function(e) {
    cat("  Uni PGLS failed:", conditionMessage(e), "\n"); NULL
  })

  # --------------------------------------------------------------------------
  # 5g. Fit multivariate PGLS / PIP — impute variant
  # --------------------------------------------------------------------------

  pc_multi_imp <- NULL
  if (!is.null(multi_imp_obj)) {
    pc_multi_imp <- tryCatch({
      dat_multi_imp_phy        <- dat_sp_train[train_sp, ]
      dat_multi_imp_phy[, lma_fossil_traits] <- multi_preds_imp[train_sp, ]
      complete_phy <- complete.cases(dat_multi_imp_phy[, c(target, lma_fossil_traits)])
      d_multi_imp  <- dat_multi_imp_phy[complete_phy, ]
      sp_multi_imp <- rownames(d_multi_imp)
      V_multi_imp  <- V_train[sp_multi_imp, sp_multi_imp]

      pgls_multi_imp   <- pglmEstLambda(formula = multi_formula, data = d_multi_imp, phylomat = V_multi_imp)
      lambda_multi_imp <- pgls_multi_imp$lambda
      cat("  Multi impute lambda:", round(lambda_multi_imp, 4), "\n")

      m_mi    <- model.frame(multi_formula, d_multi_imp, na.action = na.pass)
      X_mi    <- model.matrix(multi_formula, m_mi)
      b_mi    <- t(as.matrix(coef(pgls_multi_imp)))
      cv_mi   <- intersect(colnames(X_mi), rownames(b_mi))
      X_mi    <- X_mi[, cv_mi, drop = FALSE]
      b_mi    <- b_mi[cv_mi, , drop = FALSE]

      V_lam_mi    <- V_multi_imp * lambda_multi_imp
      diag(V_lam_mi) <- diag(V_multi_imp)
      eps_mi      <- d_multi_imp[[target]] - as.numeric(X_mi %*% b_mi)
      names(eps_mi) <- sp_multi_imp
      K_mi        <- solve(V_lam_mi)

      n <- length(eps_mi); p <- length(cv_mi) - 1L
      r2 <- 1 - sum(eps_mi^2) / sum((d_multi_imp[[target]] - mean(d_multi_imp[[target]]))^2)
      meta <- data.frame(method = "PGLS_multi_imp", fold = fold, target = target, stringsAsFactors = FALSE)
      model_fit_stats <<- c(model_fit_stats, list(cbind(meta, data.frame(
        r_squared = r2, adj_r_squared = 1 - (1 - r2) * (n - 1) / (n - p - 1),
        residual_se = sqrt(sum(eps_mi^2) / (n - p - 1)), n = n, lambda = lambda_multi_imp
      ), row.names = NULL)))

      list(beta = b_mi, lambda = lambda_multi_imp, K_train = K_mi,
           epsilon = eps_mi, common_vars = cv_mi, sp = sp_multi_imp,
           imp_obj = multi_imp_obj)
    }, error = function(e) {
      cat("  Multi PGLS impute failed:", conditionMessage(e), "\n"); NULL
    })
  }

  # --------------------------------------------------------------------------
  # 5h. Fit multivariate PGLS / PIP — complete-case variant
  # --------------------------------------------------------------------------

  pc_multi_cc <- tryCatch({
    complete_cc  <- complete.cases(dat_sp_train[train_sp, c(target, lma_fossil_traits)])
    d_multi_cc   <- dat_sp_train[train_sp, ][complete_cc, ]
    sp_multi_cc  <- rownames(d_multi_cc)
    V_multi_cc   <- V_train[sp_multi_cc, sp_multi_cc]

    pgls_multi_cc   <- pglmEstLambda(formula = multi_formula, data = d_multi_cc, phylomat = V_multi_cc)
    lambda_multi_cc <- pgls_multi_cc$lambda
    cat("  Multi cc lambda:", round(lambda_multi_cc, 4), "\n")

    m_mc    <- model.frame(multi_formula, d_multi_cc, na.action = na.pass)
    X_mc    <- model.matrix(multi_formula, m_mc)
    b_mc    <- t(as.matrix(coef(pgls_multi_cc)))
    cv_mc   <- intersect(colnames(X_mc), rownames(b_mc))
    X_mc    <- X_mc[, cv_mc, drop = FALSE]
    b_mc    <- b_mc[cv_mc, , drop = FALSE]

    V_lam_mc    <- V_multi_cc * lambda_multi_cc
    diag(V_lam_mc) <- diag(V_multi_cc)
    eps_mc      <- d_multi_cc[[target]] - as.numeric(X_mc %*% b_mc)
    names(eps_mc) <- sp_multi_cc
    K_mc        <- solve(V_lam_mc)

    n <- length(eps_mc); p <- length(cv_mc) - 1L
    r2 <- 1 - sum(eps_mc^2) / sum((d_multi_cc[[target]] - mean(d_multi_cc[[target]]))^2)
    meta <- data.frame(method = "PGLS_multi_cc", fold = fold, target = target, stringsAsFactors = FALSE)
    model_fit_stats <<- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = r2, adj_r_squared = 1 - (1 - r2) * (n - 1) / (n - p - 1),
      residual_se = sqrt(sum(eps_mc^2) / (n - p - 1)), n = n, lambda = lambda_multi_cc
    ), row.names = NULL)))

    list(beta = b_mc, lambda = lambda_multi_cc, K_train = K_mc,
         epsilon = eps_mc, common_vars = cv_mc, sp = sp_multi_cc)
  }, error = function(e) {
    cat("  Multi PGLS cc failed:", conditionMessage(e), "\n"); NULL
  })

  # --------------------------------------------------------------------------
  # 5i. Predict held-out sites
  # --------------------------------------------------------------------------

  for (s in held_sites) {
    obs <- as.numeric(site_obs[s, "obs_log10_lma"])
    if (is.na(obs)) next

    held_morpho <- morpho_data[morpho_data$site == s, ]
    if (nrow(held_morpho) == 0) next

    # Within-site species means for held-out site
    held_sp_agg <- held_morpho %>%
      group_by(genusSpecies) %>%
      summarise(across(all_of(intersect(c(uni_pred, lma_fossil_traits), names(held_morpho))),
                       \(x) mean(x, na.rm = TRUE)),
                .groups = "drop") %>%
      as.data.frame()
    rownames(held_sp_agg) <- held_sp_agg$genusSpecies
    held_sp_names         <- held_sp_agg$genusSpecies

    # -- Univariate LM prediction --
    pred_uni_lm <- if (!is.null(lm_uni)) {
      uni_complete <- !is.na(held_sp_agg[[uni_pred]])
      if (any(uni_complete))
        tryCatch(mean(as.numeric(predict(lm_uni, newdata = held_sp_agg[uni_complete, ])),
                      na.rm = TRUE),
                 error = function(e) NA_real_)
      else NA_real_
    } else NA_real_

    # -- Multivariate LM predictions (impute held-out traits using training model) --
    pred_multi_lm_imp <- NA_real_
    if (!is.null(lm_multi_imp) && !is.null(multi_imp_obj)) {
      tryCatch({
        held_imp        <- predict(multi_imp_obj, held_sp_agg[, lma_fossil_traits, drop = FALSE])
        rownames(held_imp) <- held_sp_names
        pred_multi_lm_imp <- mean(as.numeric(predict(lm_multi_imp, newdata = held_imp)),
                                   na.rm = TRUE)
      }, error = function(e) NULL)
    }

    pred_multi_lm_cc <- NA_real_
    if (!is.null(lm_multi_cc)) {
      tryCatch({
        cc_held <- complete.cases(held_sp_agg[, lma_fossil_traits])
        if (any(cc_held))
          pred_multi_lm_cc <- mean(as.numeric(
            predict(lm_multi_cc, newdata = held_sp_agg[cc_held, ])), na.rm = TRUE)
      }, error = function(e) NULL)
    }

    # -- Univariate PIP prediction --
    pred_uni_pgls <- NA_real_
    pred_uni_pip  <- NA_real_
    if (!is.null(pc_uni)) {
      tryCatch({
        uni_complete_held <- !is.na(held_sp_agg[[uni_pred]])
        held_uni          <- held_sp_agg[uni_complete_held, , drop = FALSE]
        held_uni[[target]] <- 0
        mm_uni    <- model.matrix(uni_formula, model.frame(uni_formula, held_uni))
        avail_uni <- intersect(pc_uni$common_vars, colnames(mm_uni))
        X_h_uni   <- mm_uni[, avail_uni, drop = FALSE]
        b_h_uni   <- pc_uni$beta[avail_uni, , drop = FALSE]

        y_pgls_uni        <- as.numeric(X_h_uni %*% b_h_uni)
        names(y_pgls_uni) <- rownames(held_uni)
        pred_uni_pgls     <- mean(y_pgls_uni, na.rm = TRUE)

        y_pip_uni <- y_pgls_uni
        in_vcv    <- intersect(names(y_pgls_uni), rownames(full_vcv))
        if (length(in_vcv) > 0) {
          C_cross       <- full_vcv[in_vcv, pc_uni$sp, drop = FALSE] * pc_uni$lambda
          pip_adj       <- as.numeric(C_cross %*% pc_uni$K_train %*% pc_uni$epsilon[pc_uni$sp])
          names(pip_adj) <- in_vcv
          y_pip_uni[in_vcv] <- y_pgls_uni[in_vcv] + pip_adj
        }
        pred_uni_pip <- mean(y_pip_uni, na.rm = TRUE)
      }, error = function(e) NULL)
    }

    # -- Multivariate PIP — impute --
    pred_multi_pgls_imp <- NA_real_
    pred_multi_pip_imp  <- NA_real_
    if (!is.null(pc_multi_imp)) {
      tryCatch({
        held_mi_raw        <- held_sp_agg[, lma_fossil_traits, drop = FALSE]
        held_mi            <- predict(pc_multi_imp$imp_obj, held_mi_raw)
        rownames(held_mi)  <- held_sp_names
        held_mi[[target]]  <- 0
        mm_mi   <- model.matrix(multi_formula, model.frame(multi_formula, held_mi))
        av_mi   <- intersect(pc_multi_imp$common_vars, colnames(mm_mi))
        X_h_mi  <- mm_mi[, av_mi, drop = FALSE]
        b_h_mi  <- pc_multi_imp$beta[av_mi, , drop = FALSE]

        y_pgls_mi        <- as.numeric(X_h_mi %*% b_h_mi)
        names(y_pgls_mi) <- held_sp_names
        pred_multi_pgls_imp <- mean(y_pgls_mi, na.rm = TRUE)

        y_pip_mi <- y_pgls_mi
        in_vcv   <- intersect(held_sp_names, rownames(full_vcv))
        if (length(in_vcv) > 0) {
          C_cross      <- full_vcv[in_vcv, pc_multi_imp$sp, drop = FALSE] * pc_multi_imp$lambda
          pip_adj      <- as.numeric(C_cross %*% pc_multi_imp$K_train %*% pc_multi_imp$epsilon[pc_multi_imp$sp])
          names(pip_adj) <- in_vcv
          y_pip_mi[in_vcv] <- y_pgls_mi[in_vcv] + pip_adj
        }
        pred_multi_pip_imp <- mean(y_pip_mi, na.rm = TRUE)
      }, error = function(e) NULL)
    }

    # -- Multivariate PIP — complete-case --
    pred_multi_pgls_cc <- NA_real_
    pred_multi_pip_cc  <- NA_real_
    if (!is.null(pc_multi_cc)) {
      tryCatch({
        cc_held_mc  <- complete.cases(held_sp_agg[, lma_fossil_traits])
        held_mc     <- held_sp_agg[cc_held_mc, lma_fossil_traits, drop = FALSE]
        held_mc_sp  <- rownames(held_mc)
        held_mc[[target]] <- 0
        mm_mc   <- model.matrix(multi_formula, model.frame(multi_formula, held_mc))
        av_mc   <- intersect(pc_multi_cc$common_vars, colnames(mm_mc))
        X_h_mc  <- mm_mc[, av_mc, drop = FALSE]
        b_h_mc  <- pc_multi_cc$beta[av_mc, , drop = FALSE]

        y_pgls_mc        <- as.numeric(X_h_mc %*% b_h_mc)
        names(y_pgls_mc) <- held_mc_sp
        pred_multi_pgls_cc <- mean(y_pgls_mc, na.rm = TRUE)

        y_pip_mc <- y_pgls_mc
        in_vcv   <- intersect(held_mc_sp, rownames(full_vcv))
        if (length(in_vcv) > 0) {
          C_cross      <- full_vcv[in_vcv, pc_multi_cc$sp, drop = FALSE] * pc_multi_cc$lambda
          pip_adj      <- as.numeric(C_cross %*% pc_multi_cc$K_train %*% pc_multi_cc$epsilon[pc_multi_cc$sp])
          names(pip_adj) <- in_vcv
          y_pip_mc[in_vcv] <- y_pgls_mc[in_vcv] + pip_adj
        }
        pred_multi_pip_cc <- mean(y_pip_mc, na.rm = TRUE)
      }, error = function(e) NULL)
    }

    cv_results[[s]] <- list(
      site                = s,
      fold                = fold,
      obs                 = obs,
      pred_uni_lm         = pred_uni_lm,
      pred_uni_pgls       = pred_uni_pgls,
      pred_uni_pip        = pred_uni_pip,
      pred_multi_lm_imp   = pred_multi_lm_imp,
      pred_multi_lm_cc    = pred_multi_lm_cc,
      pred_multi_pgls_imp = pred_multi_pgls_imp,
      pred_multi_pip_imp  = pred_multi_pip_imp,
      pred_multi_pgls_cc  = pred_multi_pgls_cc,
      pred_multi_pip_cc   = pred_multi_pip_cc
    )

    cat("    Site:", s,
        "| obs:", round(obs, 3),
        "| uni_lm:", round(pred_uni_lm, 3),
        "| uni_pip:", round(pred_uni_pip, 3),
        "| multi_pip_imp:", round(pred_multi_pip_imp, 3), "\n")
  }

  fold_time <- proc.time()["elapsed"] - fold_start
  cat("  Fold", fold, "complete in", round(fold_time, 1), "s\n")
}

# 6. COMPILE RESULTS

cat("\nCompiling results...\n")

results_df <- do.call(rbind, lapply(cv_results, function(r) {
  data.frame(
    site                = r$site,
    fold                = r$fold,
    obs                 = r$obs,
    pred_uni_lm         = r$pred_uni_lm,
    pred_uni_pgls       = r$pred_uni_pgls,
    pred_uni_pip        = r$pred_uni_pip,
    pred_multi_lm_imp   = r$pred_multi_lm_imp,
    pred_multi_lm_cc    = r$pred_multi_lm_cc,
    pred_multi_pgls_imp = r$pred_multi_pgls_imp,
    pred_multi_pip_imp  = r$pred_multi_pip_imp,
    pred_multi_pgls_cc  = r$pred_multi_pgls_cc,
    pred_multi_pip_cc   = r$pred_multi_pip_cc,
    stringsAsFactors    = FALSE
  )
}))
rownames(results_df) <- NULL

# 7. RMSE TABLE

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

pred_cols <- c(
  "pred_uni_lm", "pred_uni_pgls", "pred_uni_pip",
  "pred_multi_lm_imp", "pred_multi_lm_cc",
  "pred_multi_pgls_imp", "pred_multi_pip_imp",
  "pred_multi_pgls_cc", "pred_multi_pip_cc"
)

rmse_rows <- data.frame(
  model   = pred_cols,
  target  = target,
  rmse    = sapply(pred_cols, function(col) rmse(results_df$obs, results_df[[col]])),
  n_sites = sapply(pred_cols, function(col) sum(!is.na(results_df[[col]]))),
  row.names = NULL
)
rmse_rows <- rmse_rows[order(rmse_rows$rmse), ]

cat("\n10-fold LOSO CV RMSE (log10 LMA, site level):\n")
print(rmse_rows, row.names = FALSE)

# 8. WRITE OUTPUTS

model_coefs_df     <- do.call(rbind, Filter(Negate(is.null), model_coefs))
model_fit_stats_df <- do.call(rbind, Filter(Negate(is.null), model_fit_stats))
rownames(model_coefs_df)     <- NULL
rownames(model_fit_stats_df) <- NULL

write.csv(results_df,         "tables/lma_cv_site_predictions.csv", row.names = FALSE)
write.csv(rmse_rows,          "tables/lma_cv_rmse.csv",             row.names = FALSE)
write.csv(model_coefs_df,     "tables/lma_cv_model_coefs.csv",      row.names = FALSE)
write.csv(model_fit_stats_df, "tables/lma_cv_model_fit.csv",        row.names = FALSE)

cat("\nSaved tables/lma_cv_site_predictions.csv\n")
cat("Saved tables/lma_cv_rmse.csv\n")
cat("Saved tables/lma_cv_model_coefs.csv\n")
cat("Saved tables/lma_cv_model_fit.csv\n")

# 9. SAME-SITE COMPARISON

# The CC multivariate models only predict for the 63 sites with complete trait
# data. To fairly compare all models, compute RMSE restricted to those sites.

cc_sites <- !is.na(results_df$pred_multi_pip_cc)
cat("\nSame-site comparison: RMSE restricted to the", sum(cc_sites),
    "sites with complete multivariate CC predictions:\n")

rmse_cc63 <- data.frame(
  model        = pred_cols,
  target       = target,
  rmse_all     = sapply(pred_cols, function(col) rmse(results_df$obs, results_df[[col]])),
  n_sites_all  = sapply(pred_cols, function(col) sum(!is.na(results_df[[col]]))),
  rmse_cc63    = sapply(pred_cols, function(col)
                   rmse(results_df$obs[cc_sites], results_df[[col]][cc_sites])),
  n_sites_cc63 = sum(cc_sites),
  row.names    = NULL
)
rmse_cc63 <- rmse_cc63[order(rmse_cc63$rmse_cc63), ]

print(rmse_cc63, row.names = FALSE)

write.csv(rmse_cc63, "tables/lma_cv_rmse_same_site.csv", row.names = FALSE)
cat("Saved tables/lma_cv_rmse_same_site.csv\n")
