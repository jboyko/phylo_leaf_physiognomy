setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ape)
library(caret)
source("code/Phylogenetically-Informed_Predictions_Source.R")

# ==============================================================================
# 1. LOAD DATA AND PHYLOGENY
# ==============================================================================

# tre_lma_pruned.tre is written by 00b_data_cleaning.R with BUILD_PHYLOGENY=TRUE
phy <- read.tree("data/tre_lma_pruned.tre")
dat <- read.csv("data/lma_species.csv")

rownames(dat) <- dat$genusSpecies

phylomat <- vcv(phy)
diag(phylomat) <- diag(phylomat) + 1e-6

# ==============================================================================
# 2. FIT PGLS — COMPLETE-CASE ONLY
# ==============================================================================

# Single predictor: bagImpute is meaningless with one variable, so complete-case
# is the only sensible option here.
pred_name <- "log10_petiole_metric"
target    <- "log10_lma"
formula   <- as.formula(paste(target, "~", pred_name))

cols          <- c(target, pred_name)
d             <- dat[, cols, drop = FALSE]
complete_rows <- complete.cases(d)
d_fit         <- d[complete_rows, , drop = FALSE]
d_fit         <- d_fit[rownames(d_fit) %in% phy$tip.label, , drop = FALSE]
phy_fit       <- keep.tip(phy, rownames(d_fit))
pm_fit        <- vcv(phy_fit)
diag(pm_fit)  <- diag(pm_fit) + 1e-6
d_fit         <- d_fit[rownames(pm_fit), , drop = FALSE]

cat("N species:", nrow(d_fit), "\n")

pgls_fit <- pglmEstLambda(formula = formula, data = d_fit, phylomat = pm_fit)
cat("lambda:", pgls_fit$lambda, "\n")

m    <- model.frame(formula, d_fit, na.action = na.pass)
X    <- model.matrix(formula, m)
# coef() returns a 1×p data frame; transpose to p×1 for X %*% beta
beta   <- t(as.matrix(coef(pgls_fit)))   # p × 1, rownames = predictor names
beta_t <- t(beta)                         # 1 × p, stored for PIP application

lambda <- pgls_fit$lambda
V_lam  <- pm_fit * lambda
diag(V_lam) <- diag(pm_fit)

resid        <- d_fit[[target]] - as.numeric(X %*% beta)
names(resid) <- rownames(d_fit)

# ==============================================================================
# 3. SAVE PIP COMPONENTS
# ==============================================================================

lma_pip <- list(
  pred_names  = pred_name,
  tree_pruned = phy_fit,
  beta        = beta_t,
  lambda      = lambda,
  V_lam       = V_lam,
  resid       = resid,
  X           = X,
  formula     = formula,
  vcv         = pgls_fit$vcv,
  dat_fit     = d_fit,
  phylomat    = pm_fit
)

saveRDS(lma_pip, file = "models/lma_pip_components.rds")
cat("Saved models/lma_pip_components.rds\n")
