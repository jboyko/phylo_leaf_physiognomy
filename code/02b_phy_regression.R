setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# 1. LOAD DATA AND PHYLOGENY
# ==============================================================================

phy        <- read.tree("data/tre_lma_pruned.tre")
dat        <- read.csv("data/lma_species.csv")
lma_nophy  <- readRDS("models/lma_nophy_models.rds")

rownames(dat) <- dat$genusSpecies

phylomat <- vcv(phy)
diag(phylomat) <- diag(phylomat) + 1e-6

target <- "log10_lma"

# ==============================================================================
# 2. UNIVARIATE PGLS: log10_lma ~ log10_pw2a_ratio  (complete-case only)
# ==============================================================================

# bagImpute is meaningless with a single predictor, so complete-case is the
# only sensible option here.

uni_pred  <- "log10_pw2a_ratio"
uni_form  <- as.formula(paste(target, "~", uni_pred))

d_uni           <- dat[, c(target, uni_pred), drop = FALSE]
complete_uni    <- complete.cases(d_uni)
d_uni_fit       <- d_uni[complete_uni, , drop = FALSE]
d_uni_fit       <- d_uni_fit[rownames(d_uni_fit) %in% phy$tip.label, , drop = FALSE]
phy_uni         <- keep.tip(phy, rownames(d_uni_fit))
pm_uni          <- vcv(phy_uni)
diag(pm_uni)    <- diag(pm_uni) + 1e-6
d_uni_fit       <- d_uni_fit[rownames(pm_uni), , drop = FALSE]

cat("Univariate N species:", nrow(d_uni_fit), "\n")

pgls_uni <- pglmEstLambda(formula = uni_form, data = d_uni_fit, phylomat = pm_uni)
cat("Univariate lambda:", pgls_uni$lambda, "\n")

X_uni        <- model.matrix(uni_form, model.frame(uni_form, d_uni_fit, na.action = na.pass))
beta_uni     <- t(as.matrix(coef(pgls_uni)))
lambda_uni   <- pgls_uni$lambda
V_lam_uni    <- pm_uni * lambda_uni
diag(V_lam_uni) <- diag(pm_uni)
resid_uni    <- d_uni_fit[[target]] - as.numeric(X_uni %*% beta_uni)
names(resid_uni) <- rownames(d_uni_fit)

pip_univariate <- list(
  pred_names = uni_pred,
  formula    = uni_form,
  beta       = beta_uni,
  lambda     = lambda_uni,
  V_lam      = V_lam_uni,
  resid      = resid_uni,
  X          = X_uni,
  vcv        = pgls_uni$vcv,
  dat_fit    = d_uni_fit,
  phylomat   = pm_uni,
  tree       = phy_uni
)

# ==============================================================================
# 3. MULTIVARIATE PGLS: log10_lma ~ all 12 fossil-measurable traits
# ==============================================================================

# Use the same predictor set as the multivariate LM in 01b for direct comparison.
pred_names <- lma_nophy$multivariate$pred_names
cat("\nMultivariate predictors (", length(pred_names), "):",
    paste(pred_names, collapse = ", "), "\n")

multi_formula <- as.formula(paste(target, "~", paste(pred_names, collapse = " + ")))

pgls_configs <- list(
  impute = list(impute = TRUE,  desc = "bagImpute"),
  cc     = list(impute = FALSE, desc = "complete-case")
)

pip_multi_configs <- list()

for (cfg_name in names(pgls_configs)) {
  cfg <- pgls_configs[[cfg_name]]
  cat("\n=== Multivariate PGLS,", cfg$desc, "===\n")

  cols <- c(target, pred_names)
  d    <- dat[, cols, drop = FALSE]

  if (cfg$impute) {
    imp_obj         <- preProcess(d, method = "bagImpute")
    d_fit           <- predict(imp_obj, d)
    rownames(d_fit) <- rownames(d)   # predict.preProcess may drop rownames
    imp_traits      <- preProcess(dat[, pred_names, drop = FALSE], method = "bagImpute")
    phy_fit         <- phy
    pm_fit          <- phylomat
    in_both         <- intersect(rownames(d_fit), rownames(pm_fit))
    d_fit           <- d_fit[in_both, , drop = FALSE]
    pm_fit          <- pm_fit[in_both, in_both]
  } else {
    imp_obj       <- NULL
    imp_traits    <- NULL
    complete_rows <- complete.cases(d)
    d_fit         <- d[complete_rows, , drop = FALSE]
    in_tree       <- intersect(rownames(d_fit), rownames(phylomat))
    d_fit         <- d_fit[in_tree, , drop = FALSE]
    phy_fit       <- keep.tip(phy, in_tree)
    pm_fit        <- vcv(phy_fit)
    diag(pm_fit)  <- diag(pm_fit) + 1e-6
    d_fit         <- d_fit[rownames(pm_fit), , drop = FALSE]
  }

  cat("N species:", nrow(d_fit), "\n")

  pgls_fit <- pglmEstLambda(formula = multi_formula, data = d_fit, phylomat = pm_fit)
  cat("lambda:", pgls_fit$lambda, "\n")

  m           <- model.frame(multi_formula, d_fit, na.action = na.pass)
  X           <- model.matrix(multi_formula, m)
  beta_raw    <- coef(pgls_fit)
  beta        <- t(as.matrix(beta_raw))
  common_vars <- intersect(colnames(X), rownames(beta))
  X           <- X[, common_vars, drop = FALSE]
  beta        <- beta[common_vars, , drop = FALSE]

  lambda <- pgls_fit$lambda
  V_lam  <- pm_fit * lambda
  diag(V_lam) <- diag(pm_fit)

  resid        <- d_fit[[target]] - as.numeric(X %*% beta)
  names(resid) <- rownames(d_fit)

  pip_multi_configs[[cfg_name]] <- list(
    beta                = beta,
    lambda              = lambda,
    V_lam               = V_lam,
    resid               = resid,
    X                   = X,
    formula             = multi_formula,
    vcv                 = pgls_fit$vcv,
    dat_fit             = d_fit,
    impute_model        = imp_obj,
    impute_model_traits = imp_traits,
    tree                = phy_fit,
    phylomat            = pm_fit
  )
}

# ==============================================================================
# 4. SAVE PIP COMPONENTS
# ==============================================================================

lma_pip <- list(
  univariate   = pip_univariate,
  multivariate = list(
    configs    = pip_multi_configs,
    pred_names = pred_names
  )
)

saveRDS(lma_pip, file = "models/lma_pip_components.rds")
cat("\nSaved models/lma_pip_components.rds\n")
