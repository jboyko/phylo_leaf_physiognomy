source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

# Generates the fitted model objects for the dilp package's dilp_pgls()
# function, including the calibration quantities for site intervals.
# Outputs land in dilp_update/ — copy them into the package as described
# in the printed instructions at the end.

dir.create("dilp_update", showWarnings = FALSE)

# ==============================================================================
# 1. LOAD PIPELINE OUTPUTS
# ==============================================================================

pip  <- readRDS("models/pip_components.rds")
rmse <- read.csv("tables/loso_cv_rmse.csv", stringsAsFactors = FALSE)
source("code/pip_uncertainty.R")
site_discrepancy <- pip_load_site_discrepancy()

# ==============================================================================
# 2. REGRESSION COEFFICIENTS (pgls_beta_mat, pgls_beta_map)
#    Stored as named numeric vectors so X %*% pgls_beta_* works in dilp_pgls().
#    Source: impute variant of PGLS fit (same as pip_components backward-compat
#    top-level keys).
# ==============================================================================

# Pipeline uses renamed columns; dilp_pgls() expects raw dilp names.
# Map pipeline names back to dilp column names (inverse of rename_pipeline()).
pipeline_to_dilp <- c(
  pw2.a.ratio                 = "petiole_metric",
  ln.leaf.area.mm2            = "ln_leaf_area",
  feret.diam.ratio            = "fdr",
  margin.score                = "margin",
  perim.ratio                 = "perimeter_ratio",
  teeth.perimeter.percm       = "tc_p",
  teeth.interior.percm        = "tc_ip",
  avt.tooth.area              = "avg_ta",
  tooth.area.blade.area.ratio = "ta_ba",
  tooth.area.perimeter        = "ta_p",
  tooth.area.interior         = "ta_ip",
  teeth.blade.area.ratio      = "tc_ba"
)

dilp_column_names <- function(nms) {
  ifelse(nms %in% names(pipeline_to_dilp), pipeline_to_dilp[nms], nms)
}

pgls_beta_mat <- setNames(as.numeric(pip$beta_mat),
                         dilp_column_names(rownames(pip$beta_mat)))
pgls_beta_map <- setNames(as.numeric(pip$beta_map),
                         dilp_column_names(rownames(pip$beta_map)))

# Fossils have no response variable, so retain the traits-only imputers fitted
# during calibration. The MAT and MAP fitters are intentionally separate.
pgls_impute_model_mat_traits <- pip$impute_model_mat_traits
pgls_impute_model_map_traits <- pip$impute_model_map_traits

cat("pgls_beta_mat predictors:", paste(names(pgls_beta_mat), collapse = ", "), "\n")
cat("pgls_beta_map predictors:", paste(names(pgls_beta_map), collapse = ", "), "\n\n")

saveRDS(pgls_beta_mat, "dilp_update/pgls_beta_mat.rds")
saveRDS(pgls_beta_map, "dilp_update/pgls_beta_map.rds")
cat("Saved pgls_beta_mat.rds and pgls_beta_map.rds\n")

# ==============================================================================
# 3. PHYLOGENETIC WEIGHT VECTORS (pgls_Ke_mat, pgls_Ke_map)
#    Ke = solve(V_lam) %*% resid
#    Pre-multiplied so dilp_pgls() only needs t(V_cross) %*% Ke at prediction
#    time. Named by extant training species (rownames of V_lam).
# ==============================================================================

cat("\nComputing Ke vectors (solve(V_lam) %*% resid) ...\n")

K_mat <- solve(pip$V_lam_mat)
K_map <- solve(pip$V_lam_map)
Ke_mat <- K_mat %*% pip$resid_mat
Ke_map <- K_map %*% pip$resid_map

pgls_Ke_mat <- setNames(as.numeric(Ke_mat), rownames(Ke_mat))
pgls_Ke_map <- setNames(as.numeric(Ke_map), rownames(Ke_map))

cat("  pgls_Ke_mat: length", length(pgls_Ke_mat),
    "| range", round(range(pgls_Ke_mat), 4), "\n")
cat("  pgls_Ke_map: length", length(pgls_Ke_map),
    "| range", round(range(pgls_Ke_map), 4), "\n\n")

saveRDS(pgls_Ke_mat, "dilp_update/pgls_Ke_mat.rds")
saveRDS(pgls_Ke_map, "dilp_update/pgls_Ke_map.rds")
cat("Saved pgls_Ke_mat.rds and pgls_Ke_map.rds\n")

# ==============================================================================
# 4. SCAFFOLD TREE (inst/extdata/tre_scaffold.tre)
#    Built by 00_data_cleaning.R. The dilp package loads this via
#    system.file("extdata", "tre_scaffold.tre", package = "dilp").
#    Copy whenever the training species list or WCVP tree changes.
# ==============================================================================

file.copy("data/tre_scaffold.tre", "dilp_update/tre_scaffold.tre", overwrite = TRUE)
cat("Copied tre_scaffold.tre\n")

# Stage only the univariate LMA fit used in the fossil analysis, with the
# original placement tree used by 04b. The grafted covariance is calculated
# at prediction time, as in that analysis.
lma_fit <- readRDS("models/lma_pip_components.rds")$univariate
lma_taxa <- rownames(lma_fit$dat_fit)
lma_tree <- ape::read.tree("data/tre_lma_pruned.tre")
stopifnot(all(lma_taxa %in% lma_tree$tip.label),
          identical(colnames(lma_fit$X), c("(Intercept)", "log10_pw2a_ratio")))
lma_pgls_calibration <- list(
  beta = as.numeric(lma_fit$beta),
  lambda = lma_fit$lambda,
  weights = setNames(as.numeric(solve(lma_fit$V_lam[lma_taxa, lma_taxa],
                                       lma_fit$resid[lma_taxa])), lma_taxa),
  tree = lma_tree
)
saveRDS(lma_pgls_calibration, "dilp_update/lma_pgls_calibration.rds", compress = "xz")
cat("Staged univariate LMA calibration for", length(lma_taxa), "extant species\n")

# ==============================================================================
# 5. UNCERTAINTY CALIBRATION FROM THE FITTED MODELS AND SITE CV
#    The package reconstructs the calibration covariance from the rooted
#    scaffold. Store the design and small fitted summaries rather than K.
# ==============================================================================

stage_uncertainty <- function(X, K, residuals, target, beta) {
  ids <- rownames(X)
  fit <- pip_uncertainty_fit(X, K[ids, ids, drop = FALSE], residuals[ids])
  colnames(X) <- dilp_column_names(colnames(X))
  stopifnot(identical(colnames(X), names(beta)))
  list(X = X, sigma2 = fit$sigma2, beta_vcov = unname(fit$beta_vcov),
       tau2 = site_discrepancy[target, "tau2"])
}

pgls_uncertainty_mat <- stage_uncertainty(
  pip$X_mat, K_mat, pip$resid_mat, "mat", pgls_beta_mat
)
pgls_uncertainty_map <- stage_uncertainty(
  pip$X_map, K_map, pip$resid_map, "log_map", pgls_beta_map
)

# Retain the older RMSE-based error summaries for existing package callers.
rmse_mat_pip <- rmse$rmse[rmse$column == "pip_sp_site_impute_mat"]
rmse_map_pip <- rmse$rmse[rmse$column == "pip_sp_site_impute_log_map"]
if (length(rmse_mat_pip) != 1L || length(rmse_map_pip) != 1L ||
    !is.finite(rmse_mat_pip) || !is.finite(rmse_map_pip)) {
  stop("Could not find one finite impute PIP MAT and log(MAP) RMSE in ",
       "tables/loso_cv_rmse.csv. Regenerate final CV outputs before staging dilp.")
}
pgls_rmse_mat <- unname(rmse_mat_pip)
pgls_rmse_log_map <- unname(rmse_map_pip)

cat("\nNew RMSE constants:\n")
cat("  MAT RMSE (pgls_Ke_mat) :", round(rmse_mat_pip, 4), "\n")
cat("  log(MAP) RMSE          :", round(rmse_map_pip, 4), "\n")

# ==============================================================================
# 6. PATCH sysdata.rda
#    Update only PIP objects in the package's existing sysdata and
#    save to dilp_update/sysdata.rda.
#    Then copy that file to R/sysdata.rda in the dilp package.
# ==============================================================================

DILP_SYSDATA <- Sys.getenv("DILP_SYSDATA", "../dilp/R/sysdata.rda")

if (!file.exists(DILP_SYSDATA)) {
  warning("Could not find ", DILP_SYSDATA,
          " — skipping sysdata patch. Set DILP_SYSDATA to the correct path.")
} else {
  sysdata_env <- new.env(parent = emptyenv())
  load(DILP_SYSDATA, envir = sysdata_env)
  cat("\nLoaded", DILP_SYSDATA, "(", length(ls(sysdata_env)), "objects)\n")

  sysdata_env$pgls_beta_mat <- pgls_beta_mat
  sysdata_env$pgls_beta_map <- pgls_beta_map
  sysdata_env$pgls_Ke_mat   <- pgls_Ke_mat
  sysdata_env$pgls_Ke_map   <- pgls_Ke_map
  sysdata_env$pgls_lambda_mat <- pip$lambda_mat
  sysdata_env$pgls_lambda_map <- pip$lambda_map
  sysdata_env$pgls_name_table <- pip$name_table_full
  sysdata_env$pgls_pred_names <- unname(pipeline_to_dilp[pip$pred_names])
  sysdata_env$pgls_impute_model_mat_traits <- pgls_impute_model_mat_traits
  sysdata_env$pgls_impute_model_map_traits <- pgls_impute_model_map_traits
  sysdata_env$pgls_rmse_mat <- pgls_rmse_mat
  sysdata_env$pgls_rmse_log_map <- pgls_rmse_log_map
  sysdata_env$pgls_uncertainty_mat <- pgls_uncertainty_mat
  sysdata_env$pgls_uncertainty_map <- pgls_uncertainty_map

  save(list = ls(sysdata_env), envir = sysdata_env,
       file = "dilp_update/sysdata.rda", compress = "xz")
  cat("Saved dilp_update/sysdata.rda (", length(ls(sysdata_env)),
      "objects — existing objects preserved)\n")
}

# ==============================================================================
# 7. PRINT REMAINING MANUAL STEPS
# ==============================================================================

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("REMAINING MANUAL STEPS\n")
cat(rep("=", 70), "\n", sep = "")
cat("
1. Copy dilp_update/sysdata.rda  ->  R/sysdata.rda  in the dilp package.

2. Copy dilp_update/tre_scaffold.tre  ->  inst/extdata/tre_scaffold.tre.

3. Copy dilp_update/lma_pgls_calibration.rds  ->  inst/extdata/lma_pgls_calibration.rds.

4. Run devtools::document() and devtools::check() in the package checkout.
")
cat(rep("=", 70), "\n", sep = "")
