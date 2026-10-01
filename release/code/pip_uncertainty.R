# Conditional universal-kriging prediction-error covariance. Lambda, tree,
# ages and completed predictor values are held fixed. Responses are MAT or
# log(MAP); this is not yet a calibrated interval for actual site climate.
pip_uncertainty_fit <- function(X, K, residuals) {
  stopifnot(is.matrix(X), is.matrix(K), nrow(X) == nrow(K),
            ncol(K) == nrow(K), length(residuals) == nrow(X),
            all(is.finite(X)), all(is.finite(K)), all(is.finite(residuals)))
  df <- nrow(X) - ncol(X)
  if (df <= 0 || qr(X)$rank != ncol(X)) stop("Unidentified regression design")
  sigma2 <- as.numeric(crossprod(residuals, K %*% residuals)) / df
  if (!is.finite(sigma2) || sigma2 < 0) stop("Invalid residual variance")
  list(X = X, K = K, df = df, sigma2 = sigma2,
       beta_vcov = sigma2 * solve(crossprod(X, K %*% X)))
}

pip_prediction_covariance <- function(fit, X_new, C, V_new) {
  stopifnot(identical(colnames(X_new), colnames(fit$X)),
            ncol(C) == nrow(fit$X), nrow(C) == nrow(X_new),
            identical(dim(V_new), c(nrow(X_new), nrow(X_new))),
            all(is.finite(X_new)), all(is.finite(C)), all(is.finite(V_new)))
  weights <- C %*% fit$K
  A <- X_new - weights %*% fit$X
  S <- fit$sigma2 * (V_new - weights %*% t(C)) +
    A %*% fit$beta_vcov %*% t(A)
  S <- (S + t(S)) / 2
  # Do not silently repair an invalid joint covariance model.
  tol <- 1e-8 * max(1, max(abs(S)))
  if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < -tol)
    stop("Prediction-error covariance is not positive semidefinite")
  dimnames(S) <- list(rownames(X_new), rownames(X_new))
  S
}

pip_interval <- function(estimate, variance, df, level = 0.95) {
  if (any(!is.finite(variance)) || any(variance < -1e-8))
    stop("Invalid prediction variance")
  se <- sqrt(pmax(variance, 0))
  # Approximate Student intervals: lambda is estimated, then treated as fixed.
  half_width <- qt((1 + level) / 2, df) * se
  data.frame(estimate = as.numeric(estimate), se = se,
             lower = estimate - half_width, upper = estimate + half_width)
}

pip_site_interval <- function(prediction, S, df) {
  stopifnot(length(prediction) == nrow(S), nrow(S) > 0,
            all(is.finite(prediction)))
  pip_interval(mean(prediction), sum(S) / nrow(S)^2, df)
}

pip_lambda_covariance <- function(V, lambda) {
  out <- V * lambda
  diag(out) <- diag(V)
  out
}

# Shared site discrepancy variance. Held-out site errors are dominated by a
# component common to all species at a site, which the phylogenetic
# covariance cannot represent. Method of moments on out-of-fold site
# residuals: E[r^2] = tau^2 + E[se^2]. Bias is retained inside tau^2.
pip_site_discrepancy_var <- function(residual, se) {
  stopifnot(length(residual) == length(se), length(residual) > 1,
            all(is.finite(residual)), all(is.finite(se)), all(se >= 0))
  max(0, mean(residual^2) - mean(se^2))
}

# tau^2 per target from the 10-fold site-grouped CV table written by 03.
pip_load_site_discrepancy <- function(path = "tables/pip_cv_site_uncertainty.csv") {
  if (!file.exists(path)) stop(path, " not found; run code/03_loso_cv.R first")
  cv <- read.csv(path, stringsAsFactors = FALSE)
  do.call(rbind, lapply(split(cv, cv$target), function(d) {
    d <- d[is.finite(d$observed), ]
    r <- d$observed - d$estimate
    data.frame(target = d$target[1], n_sites = nrow(d),
      rmse = sqrt(mean(r^2)), mean_se2 = mean(d$se^2),
      tau2 = pip_site_discrepancy_var(r, d$se))
  }))
}

# Site interval combining shared site discrepancy with the conditional
# phylogenetic SE of this assemblage.
pip_calibrated_interval <- function(estimate, se, tau2, level = 0.95) {
  stopifnot(length(tau2) == 1, is.finite(tau2), tau2 >= 0,
            all(is.finite(se)), all(se >= 0))
  total_se <- sqrt(tau2 + se^2)
  half_width <- qnorm((1 + level) / 2) * total_se
  data.frame(total_se = total_se, calibrated_lower = estimate - half_width,
             calibrated_upper = estimate + half_width)
}
