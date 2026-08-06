source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

# Generates all files needed to update the pre-baked model objects and
# hard-coded constants in the dilp package's dilp_pgls() function.
# Outputs land in dilp_update/ — copy them into the package as described
# in the printed instructions at the end.

dir.create("dilp_update", showWarnings = FALSE)

# ==============================================================================
# 1. LOAD PIPELINE OUTPUTS
# ==============================================================================

pip  <- readRDS("models/pip_components.rds")
rmse <- read.csv("tables/loso_cv_rmse.csv", stringsAsFactors = FALSE)

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

rename_to_dilp <- function(beta_vec) {
  nms <- rownames(beta_vec)
  nms <- ifelse(nms %in% names(pipeline_to_dilp), pipeline_to_dilp[nms], nms)
  setNames(as.numeric(beta_vec), nms)
}

pgls_beta_mat <- rename_to_dilp(pip$beta_mat)
pgls_beta_map <- rename_to_dilp(pip$beta_map)

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

Ke_mat <- solve(pip$V_lam_mat) %*% pip$resid_mat
Ke_map <- solve(pip$V_lam_map) %*% pip$resid_map

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

# ==============================================================================
# 5. NEW RMSE CONSTANTS FOR HARD-CODED UNCERTAINTY COLUMNS
#    Lines ~193-195 of dilp_pgls.R use literal RMSE values. Extract from the
#    10-fold site-grouped CV results so they stay in sync with the pipeline.
# ==============================================================================

rmse_mat_pip <- rmse$rmse[rmse$column == "pip_sp_site_impute_mat"]
rmse_map_pip <- rmse$rmse[rmse$column == "pip_sp_site_impute_log_map"]

cat("\nNew RMSE constants:\n")
cat("  MAT RMSE (pgls_Ke_mat) :", round(rmse_mat_pip, 4), "\n")
cat("  log(MAP) RMSE          :", round(rmse_map_pip, 4), "\n")

# ==============================================================================
# 6. PATCH sysdata.rda
#    Load the dilp package's existing sysdata, overwrite only the four PIP
#    objects, save to dilp_update/sysdata.rda.
#    Then copy that file to R/sysdata.rda in the dilp package.
# ==============================================================================

DILP_SYSDATA <- "../dilp/R/sysdata.rda"   # adjust if dilp lives elsewhere

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

3. In dilp_pgls() (~line 193) change the hard-coded RMSE constants:
     results$MAT.PIP.error       <- 3.64
     results$MAP.PIP.error.plus  <- exp(log_map_pip + 0.54) - results$MAP.PIP
     results$MAP.PIP.error.minus <- results$MAP.PIP - exp(log_map_pip - 0.54)
   To:
")
cat(sprintf("     results$MAT.PIP.error       <- %.4f\n", rmse_mat_pip))
cat(sprintf("     results$MAP.PIP.error.plus  <- exp(log_map_pip + %.4f) - results$MAP.PIP\n", rmse_map_pip))
cat(sprintf("     results$MAP.PIP.error.minus <- results$MAP.PIP - exp(log_map_pip - %.4f)\n", rmse_map_pip))
cat("
4. devtools::document() && devtools::check()
")
cat(rep("=", 70), "\n", sep = "")
