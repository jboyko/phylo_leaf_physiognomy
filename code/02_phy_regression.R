source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(ape)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# 1. LOAD DATA AND PHYLOGENY
# ==============================================================================

phy         <- read.tree("data/tre_pruned.tre")
dat         <- read.csv("data/data_species.csv")
all_results <- readRDS("models/nophy_models.rds")

if (!"log_map" %in% names(dat)) dat$log_map <- log(dat$map)
rownames(dat) <- dat$genusSpecies  # species names as rownames throughout

phylomat <- vcv(phy)
diag(phylomat) <- diag(phylomat) + 1e-6

# ==============================================================================
# 2. BUILD FULL-TREE TAXONOMY TABLE
# ==============================================================================

# Covers all families/orders in the WCVP backbone for fossil placement in 04_.
# Written by 00_data_cleaning.R — no need to reload the full 123k-tip tree.
name_table_full <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)
cat("name_table_full:", nrow(name_table_full), "rows\n")

# ==============================================================================
# 3. SELECT TRAITS
# ==============================================================================

# Use the same predictor set as the LM in 01_ for a fair comparison.
# No ENet-based variable selection — all fossil-measurable traits that pass
# the NA threshold are included in both phylogenetic and non-phylogenetic models.
pred_names  <- all_results$pred_names
target_vars <- c("mat", "log_map")

cat("Predictors (", length(pred_names), "):", paste(pred_names, collapse = ", "), "\n")

# ==============================================================================
# 3b. REPRODUCIBILITY / SAFETY HELPERS (issues #8, #14)
# ==============================================================================

SEED <- 42  # bagImpute (bagged trees) is stochastic; seed before every call
            # for bit-reproducible imputed values and PGLS fits (issue #8).

# aggregate(..., na.rm = TRUE) upstream (00_) can leave NaN (not NA) when a
# species/trait combination has zero observations. bagImpute's predict() is
# not guaranteed to treat NaN as missing, so coerce NaN -> NA before any
# imputation step (issue #14).
nan_to_na <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) col[is.nan(col)] <- NA
    col
  })
  df
}

# Fail loudly rather than silently propagate a corrupted (non-finite) value
# into a design matrix (issue #14).
assert_finite <- function(x, label) {
  m <- as.matrix(x)
  if (any(!is.finite(m))) {
    bad <- colnames(m)[colSums(!is.finite(m)) > 0]
    stop(sprintf("%s: non-finite values remain after imputation in columns: %s",
                 label, paste(bad, collapse = ", ")))
  }
}

# ==============================================================================
# 4. FIT PGLS — IMPUTE AND COMPLETE-CASE VARIANTS
# ==============================================================================

pgls_configs <- list(
  impute = list(impute = TRUE,  desc = "bagImpute"),
  cc     = list(impute = FALSE, desc = "complete-case")
)

pip_configs <- list()

for (cfg_name in names(pgls_configs)) {
  cfg <- pgls_configs[[cfg_name]]
  cat("\n===", cfg$desc, "===\n")

  pip_configs[[cfg_name]] <- list()

  for (target in target_vars) {
    cat("  Fitting PGLS for", target, "...\n")

    cols <- c(target, pred_names)
    d    <- nan_to_na(dat[, cols, drop = FALSE])

    if (cfg$impute) {
      # Impute response + predictors for PGLS fitting
      set.seed(SEED)
      imp_obj    <- preProcess(d, method = "bagImpute")
      d_fit      <- predict(imp_obj, d)
      assert_finite(d_fit, paste0("PGLS d_fit (", cfg_name, ", ", target, ")"))
      # Traits-only imputation for fossil application (no observed response)
      traits_only <- nan_to_na(dat[, pred_names, drop = FALSE])
      set.seed(SEED)
      imp_traits <- preProcess(traits_only, method = "bagImpute")
      assert_finite(predict(imp_traits, traits_only),
                    paste0("impute_model_traits (", target, ")"))
      phy_fit    <- phy
      pm_fit     <- phylomat
      # Ensure row order matches VCV matrix
      d_fit      <- d_fit[rownames(pm_fit), , drop = FALSE]
    } else {
      imp_obj       <- NULL
      imp_traits    <- NULL
      complete_rows <- complete.cases(d)
      d_fit         <- d[complete_rows, , drop = FALSE]
      phy_fit       <- keep.tip(phy, rownames(d_fit))
      pm_fit        <- vcv(phy_fit)
      diag(pm_fit)  <- diag(pm_fit) + 1e-6
      # Ensure row order matches VCV matrix
      d_fit         <- d_fit[rownames(pm_fit), , drop = FALSE]
    }

    cat("  N species:", nrow(d_fit), "\n")

    formula  <- as.formula(paste(target, "~", paste(pred_names, collapse = " + ")))
    pgls_fit <- pglmEstLambda(formula = formula, data = d_fit, phylomat = pm_fit)
    cat("  lambda:", pgls_fit$lambda, "\n")

    # Design matrix and coefficients
    m    <- model.frame(formula, d_fit, na.action = na.pass)
    X    <- model.matrix(formula, m)
    beta_raw    <- coef(pgls_fit)
    beta        <- t(as.matrix(beta_raw))
    common_vars <- intersect(colnames(X), rownames(beta))
    X    <- X[, common_vars, drop = FALSE]
    beta <- beta[common_vars, , drop = FALSE]

    lambda <- pgls_fit$lambda
    V_lam  <- pm_fit * lambda
    diag(V_lam) <- diag(pm_fit)

    resid        <- d_fit[[target]] - as.numeric(X %*% beta)
    names(resid) <- rownames(d_fit)

    pip_configs[[cfg_name]][[target]] <- list(
      beta                = beta,
      lambda              = lambda,
      V_lam               = V_lam,
      resid               = resid,
      X                   = X,
      formula             = formula,
      vcv                 = pgls_fit$vcv,
      dat_fit             = d_fit,
      impute_model        = imp_obj,
      impute_model_traits = imp_traits,
      tree                = phy_fit,
      phylomat            = pm_fit
    )
  }
}

# ==============================================================================
# 5. SAVE PIP COMPONENTS
# ==============================================================================

taxonomy <- dat[, c("genusSpecies", "genus", "Family", "Order")]
colnames(taxonomy) <- c("species", "genus", "family", "order")

# Backward-compatible top-level keys point to the impute variant, matching
# the original pipeline. All configs available under $configs.
imp <- pip_configs$impute

pip_components <- list(
  configs         = pip_configs,
  pred_names      = pred_names,
  taxonomy        = taxonomy,
  name_table_full = name_table_full,
  tree_pruned     = phy,

  # Backward-compatible keys — impute variant
  beta_mat                = imp$mat$beta,
  lambda_mat              = imp$mat$lambda,
  V_lam_mat               = imp$mat$V_lam,
  resid_mat               = imp$mat$resid,
  X_mat                   = imp$mat$X,
  formula_mat             = imp$mat$formula,
  vcv_mat                 = imp$mat$vcv,

  beta_map                = imp$log_map$beta,
  lambda_map              = imp$log_map$lambda,
  V_lam_map               = imp$log_map$V_lam,
  resid_map               = imp$log_map$resid,
  X_map                   = imp$log_map$X,
  formula_map             = imp$log_map$formula,
  vcv_map                 = imp$log_map$vcv,

  dat_imputed_mat         = imp$mat$dat_fit,
  dat_imputed_map         = imp$log_map$dat_fit,
  impute_model_mat_traits = imp$mat$impute_model_traits,
  impute_model_map_traits = imp$log_map$impute_model_traits
)

saveRDS(pip_components, file = "models/pip_components.rds")
cat("Saved models/pip_components.rds\n")
