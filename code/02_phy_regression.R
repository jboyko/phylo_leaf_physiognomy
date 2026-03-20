# ==============================================================================
# 1. SETUP & DATA PREPARATION
# ==============================================================================

setwd("/Users/jboyko/University\ of\ Michigan\ Dropbox/James\ Boyko/James\ Boyko’s\ files/Home/phylo_leaf_physiognomy")

# imports
library(ape)
library(caret)

source("code/Phylogenetically-Informed_Predictions_Source.R")
phy <- read.tree("data/tre_pruned.tre")
dat <- read.csv("data/data_species.csv")
all_results <- readRDS("models/nophy_models.rds")
phylomat <- vcv(phy)

# extract non-zero coefficients 
mat_coefs <- coef(all_results$mat$ENet$finalModel, all_results$mat$ENet$bestTune$lambda)
map_coefs <- coef(all_results$log_map$ENet$finalModel, all_results$log_map$ENet$bestTune$lambda)

# filter for traits with a non-zero impact
active_traits_mat <- rownames(mat_coefs)[mat_coefs[,1] != 0]
active_traits_mat <- active_traits_mat[active_traits_mat != "(Intercept)"]
active_traits_map <- rownames(map_coefs)[map_coefs[,1] != 0]
active_traits_map <- active_traits_map[active_traits_map != "(Intercept)"]

mat_formula <- as.formula(paste("mat ~", paste(active_traits_mat, collapse = " + ")))
map_formula <- as.formula(paste("log_map ~", paste(active_traits_map, collapse = " + ")))

# impute data
vars_to_impute <- active_traits_mat
impute_model <- preProcess(as.data.frame(dat[, c("mat", vars_to_impute)]), 
  method = "bagImpute")
dat_imputed <- predict(impute_model, newdata = dat[, c("mat", vars_to_impute)])
rownames(dat_imputed) <- dat$Group.1
# get lambda
diag(phylomat) <- diag(phylomat) + 1e-6
mat_pgls_fit <- pglmEstLambda(formula = mat_formula, 
  data = dat_imputed, 
  phylomat = phylomat)
print(mat_pgls_fit$lambda)

# ==============================================================================
# 2. RUN PIP LOOCV 
# ==============================================================================

# 1. Prepare Design Matrix
m <- model.frame(mat_formula, dat_imputed, na.action = na.pass)
X_all <- model.matrix(mat_formula, m)

# 2. Extract and Transpose Coefficients (Fixing the dimension error)
# coef() returns a 1-row data frame; we transpose to get a k-row column vector
beta_raw <- coef(mat_pgls_fit)
beta <- t(as.matrix(beta_raw))

# 3. Align Matrix and Coefficients
# Ensure variables are in the exact same order
common_vars <- intersect(colnames(X_all), rownames(beta))
X_all <- X_all[, common_vars, drop = FALSE]
beta <- beta[common_vars, , drop = FALSE]

# 4. Prepare Phylogenetic Parameters
lambda <- mat_pgls_fit$lambda
V_lam <- phylomat * lambda
diag(V_lam) <- diag(V_lam) + (1 - lambda)

# 5. Initialize Prediction Vector
dat_imputed$mat_pip_pred <- NA

# 6. Run the LOOCV Loop
cat("Starting PIP Loop for", nrow(dat_imputed), "species...\n")
for (i in 1:nrow(dat_imputed)) {
  # Indices
  idx_miss <- i
  idx_incl <- setdiff(1:nrow(dat_imputed), i)
  # Partition VCV
  V_inv <- solve(V_lam[idx_incl, idx_incl])
  V_cov <- V_lam[idx_incl, idx_miss, drop = FALSE]
  # Residuals of Training Data
  y_obs <- dat_imputed$mat[idx_incl]
  y_hat <- X_all[idx_incl, , drop = FALSE] %*% beta
  resid_incl <- y_obs - y_hat
  # Calculate Phylogenetic Adjustment
  phylo_adj <- t(V_cov) %*% V_inv %*% resid_incl
  # PIP Prediction
  # As.numeric ensures we assign a scalar, stripping any stray matrix dimensions
  dat_imputed$mat_pip_pred[i] <- as.numeric((X_all[idx_miss, , drop = FALSE] %*% beta) + phylo_adj)
  if(i %% 100 == 0) cat("Species", i, "completed\n")
}

# ==============================================================================
# 3. GENERATE COMPARISON TABLE
# ==============================================================================

# Function to extract the best RMSE from caret objects
get_best_rmse <- function(model) {
  return(min(model$results$RMSE))
}

# 1. Extract Baseline RMSEs
rmse_lm   <- get_best_rmse(all_results$mat$LM)
rmse_enet <- get_best_rmse(all_results$mat$ENet)
rmse_rf   <- get_best_rmse(all_results$mat$RF)

# 2. Calculate PIP RMSE
rmse_pip <- sqrt(mean((dat_imputed$mat - dat_imputed$mat_pip_pred)^2))

# 3. Create Summary Table
comparison_table <- data.frame(
  Model = c("Linear Model (OLS)", "Elastic Net (ENet)", "Random Forest (RF)", "Phylo-Informed (PIP)"),
  RMSE = c(rmse_lm, rmse_enet, rmse_rf, rmse_pip),
  Type = c("Linear Equation", "Regularized Equation", "Non-Linear/ML", "Evolutionary Model")
)

# Order by Performance (Lower RMSE is better)
comparison_table <- comparison_table[order(comparison_table$RMSE), ]
print(comparison_table)

# ==============================================================================
# 4. VISUALIZATION
# ==============================================================================

# Scatterplot comparing the best non-phylo model (RF) vs PIP
par(mfrow = c(1, 2))

# Random Forest Plot (Extract predictions from caret)
rf_preds <- predict(all_results$mat$RF, dat_imputed)
plot(dat_imputed$mat, rf_preds, 
  main = "Random Forest (No Phylogeny)",
  xlab = "Observed MAT", ylab = "Predicted MAT",
  pch = 16, col = rgb(1, 0, 0, 0.2), xlim=c(0,30), ylim=c(0,30))
abline(0, 1, lwd = 2)
text(5, 25, paste("RMSE:", round(rmse_rf, 2)), col="red", adj=0)

# PIP Plot
plot(dat_imputed$mat, dat_imputed$mat_pip_pred, 
  main = "PIP (With Phylogeny)",
  xlab = "Observed MAT", ylab = "Predicted MAT",
  pch = 16, col = rgb(0, 0, 1, 0.2), xlim=c(0,30), ylim=c(0,30))
abline(0, 1, lwd = 2)
text(5, 25, paste("RMSE:", round(rmse_pip, 2)), col="blue", adj=0)
