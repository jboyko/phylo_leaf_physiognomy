setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# CONSTANTS
# ==============================================================================

K_FOLDS <- 10
SEED    <- 42

pred_name <- "log10_petiole_metric"
target    <- "log10_lma"
formula   <- as.formula(paste(target, "~", pred_name))

# ==============================================================================
# 1. LOAD DATA AND PHYLOGENY
# ==============================================================================

phy <- read.tree("data/tre_lma_pruned.tre")
dat <- read.csv("data/lma_species.csv")

rownames(dat) <- dat$genusSpecies

# Restrict to complete cases that are in the tree
dat <- dat[complete.cases(dat[, c(target, pred_name)]), ]
dat <- dat[rownames(dat) %in% phy$tip.label, ]
cat("Species for CV:", nrow(dat), "\n")

# Pre-compute full VCV once
full_vcv        <- vcv(phy)
diag(full_vcv)  <- diag(full_vcv) + 1e-6
cat("Full VCV:", nrow(full_vcv), "x", ncol(full_vcv), "\n")

# ==============================================================================
# 2. ASSIGN CV FOLDS (random, stratified by log10_lma)
# ==============================================================================

set.seed(SEED)
lma_rank        <- rank(dat[[target]], ties.method = "first")
fold_assignment <- ((lma_rank - 1) %% K_FOLDS) + 1
names(fold_assignment) <- rownames(dat)

cat("Fold sizes:", paste(table(fold_assignment), collapse = " "), "\n")

# ==============================================================================
# 3. CV LOOP
# ==============================================================================

cv_results      <- list()
model_coefs     <- list()
model_fit_stats <- list()

for (fold in seq_len(K_FOLDS)) {

  fold_start  <- proc.time()["elapsed"]
  held_sp     <- names(fold_assignment)[fold_assignment == fold]
  train_sp    <- names(fold_assignment)[fold_assignment != fold]

  cat("\n--- Fold", fold, "/", K_FOLDS,
      "| held-out:", length(held_sp), "species ---\n")

  # Pass only formula columns to pglmEstLambda — prune() calls complete.cases()
  # on the whole data frame, so extra columns with NAs would silently drop rows.
  dat_train   <- dat[train_sp, c(target, pred_name)]
  dat_held    <- dat[held_sp,  c(target, pred_name)]
  V_train     <- full_vcv[train_sp, train_sp]

  # --------------------------------------------------------------------------
  # 3a. Fit LM
  # --------------------------------------------------------------------------

  lm_fit <- tryCatch(
    lm(formula, data = dat_train),
    error = function(e) { cat("  LM failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(lm_fit)) {
    s  <- summary(lm_fit)
    ct <- as.data.frame(s$coefficients)
    names(ct) <- c("estimate", "std_error", "t_value", "p_value")
    ct$predictor <- rownames(ct); rownames(ct) <- NULL
    meta <- data.frame(method = "LM", fold = fold, target = target,
                       stringsAsFactors = FALSE)
    model_coefs     <- c(model_coefs, list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
    model_fit_stats <- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
      residual_se = s$sigma, n = nrow(lm_fit$model), lambda = NA_real_
    ), row.names = NULL)))
  }

  # --------------------------------------------------------------------------
  # 3b. Fit PGLS / PIP
  # --------------------------------------------------------------------------

  pc <- tryCatch({
    pgls_fit <- pglmEstLambda(formula = formula, data = dat_train, phylomat = V_train)
    lambda   <- pgls_fit$lambda
    cat("  lambda:", round(lambda, 4), "\n")

    m    <- model.frame(formula, dat_train, na.action = na.pass)
    X    <- model.matrix(formula, m)
    # coef() returns a 1×p data frame; transpose to p×1 for X %*% beta
    beta        <- t(as.matrix(coef(pgls_fit)))
    common_vars <- intersect(colnames(X), rownames(beta))
    X_fit       <- X[, common_vars, drop = FALSE]
    beta_fit    <- beta[common_vars, , drop = FALSE]

    V_lam       <- V_train * lambda
    diag(V_lam) <- diag(V_train)
    epsilon      <- dat_train[[target]] - as.numeric(X_fit %*% beta_fit)
    names(epsilon) <- train_sp
    K_train      <- solve(V_lam)

    # Model summary
    n <- length(epsilon); p <- length(common_vars) - 1L
    r2 <- 1 - sum(epsilon^2) / sum((dat_train[[target]] - mean(dat_train[[target]]))^2)
    ct <- tryCatch({
      raw <- coef(summary(pgls_fit))
      df  <- as.data.frame(raw)
      names(df) <- c("estimate", "std_error", "t_value", "p_value")[seq_len(ncol(df))]
      df$predictor <- rownames(df); rownames(df) <- NULL
      df
    }, error = function(e) {
      data.frame(predictor = rownames(beta_fit), estimate = as.numeric(beta_fit),
                 std_error = NA_real_, t_value = NA_real_, p_value = NA_real_)
    })
    meta <- data.frame(method = "PGLS", fold = fold, target = target, stringsAsFactors = FALSE)
    model_coefs     <<- c(model_coefs, list(cbind(meta[rep(1L, nrow(ct)), ], ct, row.names = NULL)))
    model_fit_stats <<- c(model_fit_stats, list(cbind(meta, data.frame(
      r_squared = r2, adj_r_squared = 1 - (1 - r2) * (n - 1) / (n - p - 1),
      residual_se = sqrt(sum(epsilon^2) / (n - p - 1)), n = n, lambda = lambda
    ), row.names = NULL)))

    list(beta = beta_fit, lambda = lambda, K_train = K_train,
         epsilon = epsilon, common_vars = common_vars)
  }, error = function(e) {
    cat("  PGLS failed:", conditionMessage(e), "\n"); NULL
  })

  # --------------------------------------------------------------------------
  # 3c. Predict held-out species
  # --------------------------------------------------------------------------

  for (sp in held_sp) {
    obs  <- dat_held[sp, target]
    newrow <- dat_held[sp, , drop = FALSE]

    pred_lm <- if (!is.null(lm_fit))
      tryCatch(as.numeric(predict(lm_fit, newdata = newrow)), error = function(e) NA)
    else NA

    pred_pgls <- NA
    pred_pip  <- NA

    if (!is.null(pc)) {
      tryCatch({
        # Build design matrix for held-out species
        newrow_dm          <- newrow[, pred_name, drop = FALSE]
        newrow_dm[[target]] <- 0  # dummy response for model.frame
        mm        <- model.matrix(formula, model.frame(formula, newrow_dm))
        avail     <- intersect(pc$common_vars, colnames(mm))
        X_new     <- mm[, avail, drop = FALSE]
        beta_sub  <- pc$beta[avail, , drop = FALSE]

        pred_pgls <- as.numeric(X_new %*% beta_sub)

        # PIP correction
        pred_pip <- pred_pgls
        if (sp %in% rownames(full_vcv)) {
          C_cross  <- full_vcv[sp, train_sp, drop = FALSE] * pc$lambda
          pip_adj  <- as.numeric(C_cross %*% pc$K_train %*% pc$epsilon[train_sp])
          pred_pip <- pred_pgls + pip_adj
        }
      }, error = function(e) NULL)
    }

    cv_results[[sp]] <- list(
      species   = sp,
      fold      = fold,
      obs       = obs,
      pred_lm   = pred_lm,
      pred_pgls = pred_pgls,
      pred_pip  = pred_pip
    )
  }

  fold_time <- proc.time()["elapsed"] - fold_start
  cat("  Fold", fold, "complete in", round(fold_time, 1), "s\n")
}

# ==============================================================================
# 4. COMPILE RESULTS
# ==============================================================================

cat("\nCompiling results...\n")

results_df <- do.call(rbind, lapply(cv_results, function(r) {
  data.frame(species = r$species, fold = r$fold, obs = r$obs,
             pred_lm = r$pred_lm, pred_pgls = r$pred_pgls, pred_pip = r$pred_pip,
             stringsAsFactors = FALSE)
}))
rownames(results_df) <- NULL

# ==============================================================================
# 5. RMSE TABLE
# ==============================================================================

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

rmse_rows <- data.frame(
  model  = c("LM", "PGLS", "PIP"),
  target = target,
  rmse   = c(
    rmse(results_df$obs, results_df$pred_lm),
    rmse(results_df$obs, results_df$pred_pgls),
    rmse(results_df$obs, results_df$pred_pip)
  ),
  n_species = c(
    sum(!is.na(results_df$pred_lm)),
    sum(!is.na(results_df$pred_pgls)),
    sum(!is.na(results_df$pred_pip))
  )
)
rmse_rows <- rmse_rows[order(rmse_rows$rmse), ]
rownames(rmse_rows) <- NULL

cat("\n10-fold CV RMSE (log10 LMA):\n")
print(rmse_rows, row.names = FALSE)

# ==============================================================================
# 6. WRITE OUTPUTS
# ==============================================================================

model_coefs_df     <- do.call(rbind, Filter(Negate(is.null), model_coefs))
model_fit_stats_df <- do.call(rbind, Filter(Negate(is.null), model_fit_stats))
rownames(model_coefs_df)     <- NULL
rownames(model_fit_stats_df) <- NULL

write.csv(results_df,         "tables/lma_cv_species_predictions.csv", row.names = FALSE)
write.csv(rmse_rows,          "tables/lma_cv_rmse.csv",                row.names = FALSE)
write.csv(model_coefs_df,     "tables/lma_cv_model_coefs.csv",         row.names = FALSE)
write.csv(model_fit_stats_df, "tables/lma_cv_model_fit.csv",           row.names = FALSE)

cat("\nSaved tables/lma_cv_species_predictions.csv\n")
cat("Saved tables/lma_cv_rmse.csv\n")
cat("Saved tables/lma_cv_model_coefs.csv\n")
cat("Saved tables/lma_cv_model_fit.csv\n")
