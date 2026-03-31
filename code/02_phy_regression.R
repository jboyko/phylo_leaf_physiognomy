setwd("/Users/jboyko/phylo_leaf_physiognomy")

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

phylomat <- vcv(phy)
diag(phylomat) <- diag(phylomat) + 1e-6

# ==============================================================================
# 2. BUILD FULL-TREE TAXONOMY TABLE
# ==============================================================================

# Covers all families/orders in the WCVP backbone for fossil placement in 03_.
# Written by 00_data_cleaning.R — no need to reload the full 123k-tip tree.
name_table_full <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)
cat("name_table_full:", nrow(name_table_full), "rows\n")

# ==============================================================================
# 3. SELECT ACTIVE TRAITS (FROM ELASTICNET)
# ==============================================================================

mat_coefs <- coef(all_results$mat$ENet$finalModel,
                  all_results$mat$ENet$bestTune$lambda)
map_coefs <- coef(all_results$log_map$ENet$finalModel,
                  all_results$log_map$ENet$bestTune$lambda)

active_traits_mat <- rownames(mat_coefs)[mat_coefs[, 1] != 0]
active_traits_mat <- active_traits_mat[active_traits_mat != "(Intercept)"]
active_traits_map <- rownames(map_coefs)[map_coefs[, 1] != 0]
active_traits_map <- active_traits_map[active_traits_map != "(Intercept)"]

cat("Active MAT traits (", length(active_traits_mat), "):",
    paste(active_traits_mat, collapse = ", "), "\n")
cat("Active MAP traits (", length(active_traits_map), "):",
    paste(active_traits_map, collapse = ", "), "\n")

mat_formula <- as.formula(paste("mat ~",     paste(active_traits_mat, collapse = " + ")))
map_formula <- as.formula(paste("log_map ~", paste(active_traits_map, collapse = " + ")))

# ==============================================================================
# 4. IMPUTE MISSING TRAITS
# ==============================================================================

# Full model (response + predictors): used for PGLS fitting and LOOCV in 04_
impute_model_mat <- preProcess(as.data.frame(dat[, c("mat",     active_traits_mat)]),
                               method = "bagImpute")
dat_imputed_mat  <- predict(impute_model_mat,
                            newdata = dat[, c("mat", active_traits_mat)])
rownames(dat_imputed_mat) <- dat$Group.1

impute_model_map <- preProcess(as.data.frame(dat[, c("log_map", active_traits_map)]),
                               method = "bagImpute")
dat_imputed_map  <- predict(impute_model_map,
                            newdata = dat[, c("log_map", active_traits_map)])
rownames(dat_imputed_map) <- dat$Group.1

# Traits-only model: used to impute fossil specimens (no observed response)
impute_model_mat_traits <- preProcess(as.data.frame(dat[, active_traits_mat]),
                                      method = "bagImpute")
impute_model_map_traits <- preProcess(as.data.frame(dat[, active_traits_map]),
                                      method = "bagImpute")

# ==============================================================================
# 5. FIT PGLS (MAT AND MAP)
# ==============================================================================

mat_pgls_fit <- pglmEstLambda(formula = mat_formula,
                               data    = dat_imputed_mat,
                               phylomat = phylomat)
cat("MAT lambda:", mat_pgls_fit$lambda, "\n")

map_pgls_fit <- pglmEstLambda(formula = map_formula,
                               data    = dat_imputed_map,
                               phylomat = phylomat)
cat("MAP lambda:", map_pgls_fit$lambda, "\n")

# ==============================================================================
# 6. BUILD DESIGN MATRICES AND VCV
# ==============================================================================

# ── MAT ───────────────────────────────────────────────────────────────────────
m_mat     <- model.frame(mat_formula, dat_imputed_mat, na.action = na.pass)
X_all_mat <- model.matrix(mat_formula, m_mat)

beta_mat_raw    <- coef(mat_pgls_fit)
beta_mat        <- t(as.matrix(beta_mat_raw))
common_vars_mat <- intersect(colnames(X_all_mat), rownames(beta_mat))
X_all_mat       <- X_all_mat[, common_vars_mat, drop = FALSE]
beta_mat        <- beta_mat[common_vars_mat, , drop = FALSE]

lambda_mat <- mat_pgls_fit$lambda
# lamTrans() convention used by pglmEstLambda: scale off-diagonals by lambda,
# leave diagonal unchanged. Our phylomat diagonal = root-to-tip distance (h),
# not 1, so adding (1-lambda) to the diagonal would be incorrect.
V_lam_mat <- phylomat * lambda_mat
diag(V_lam_mat) <- diag(phylomat)

# ── MAP ───────────────────────────────────────────────────────────────────────
m_map     <- model.frame(map_formula, dat_imputed_map, na.action = na.pass)
X_all_map <- model.matrix(map_formula, m_map)

beta_map_raw    <- coef(map_pgls_fit)
beta_map        <- t(as.matrix(beta_map_raw))
common_vars_map <- intersect(colnames(X_all_map), rownames(beta_map))
X_all_map       <- X_all_map[, common_vars_map, drop = FALSE]
beta_map        <- beta_map[common_vars_map, , drop = FALSE]

lambda_map <- map_pgls_fit$lambda
V_lam_map <- phylomat * lambda_map
diag(V_lam_map) <- diag(phylomat)

# ==============================================================================
# 7. SAVE PIP COMPONENTS
# ==============================================================================

resid_mat_full <- dat_imputed_mat$mat - as.numeric(X_all_mat %*% beta_mat)
names(resid_mat_full) <- rownames(dat_imputed_mat)

resid_map_full <- dat_imputed_map$log_map - as.numeric(X_all_map %*% beta_map)
names(resid_map_full) <- rownames(dat_imputed_map)

taxonomy <- dat[, c("Group.1", "Group.2", "Group.3", "Group.4")]
colnames(taxonomy) <- c("species", "genus", "family", "order")

pip_components <- list(
  # MAT model
  beta_mat              = beta_mat,
  lambda_mat            = lambda_mat,
  V_lam_mat             = V_lam_mat,
  resid_mat             = resid_mat_full,
  X_mat                 = X_all_mat,
  formula_mat           = mat_formula,
  vcv_mat               = mat_pgls_fit$vcv,

  # MAP model
  beta_map              = beta_map,
  lambda_map            = lambda_map,
  V_lam_map             = V_lam_map,
  resid_map             = resid_map_full,
  X_map                 = X_all_map,
  formula_map           = map_formula,
  vcv_map               = map_pgls_fit$vcv,

  # Shared
  dat_imputed_mat        = dat_imputed_mat,
  dat_imputed_map        = dat_imputed_map,
  impute_model_mat_traits = impute_model_mat_traits,
  impute_model_map_traits = impute_model_map_traits,
  taxonomy               = taxonomy,
  name_table_full        = name_table_full,
  tree_pruned            = phy
)
saveRDS(pip_components, file = "models/pip_components.rds")
cat("Saved models/pip_components.rds\n")
