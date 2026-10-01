source("code/pip_uncertainty.R")

# Independent responses reduce to the usual linear-regression prediction
# variance. Site averaging retains shared coefficient uncertainty.
X <- cbind(intercept = 1, x = c(-2, -1, 0, 1, 2))
K <- diag(5)
r <- c(1, -2, 2, -2, 1)
fit <- pip_uncertainty_fit(X, K, r)
Xn <- cbind(intercept = 1, x = c(0, 1))
C <- matrix(0, 2, 5)
S <- pip_prediction_covariance(fit, Xn, C, diag(2))
expected <- fit$sigma2 * (diag(2) + Xn %*% solve(crossprod(X)) %*% t(Xn))
stopifnot(isTRUE(all.equal(unname(S), unname(expected))),
          sum(S) > sum(diag(S)),
          pip_site_interval(c(2, 4), S, fit$df)$estimate == 3)

# Observing an identical latent response with identical traits leaves no
# prediction error. This also checks the coefficient correction A, rather
# than incorrectly adding X_new Var(beta) X_new' on its own.
Si <- pip_prediction_covariance(fit, X, diag(5), diag(5))
stopifnot(max(abs(Si)) < 1e-10)

# Monte Carlo check against joint training/test errors, refitting GLS beta
# each time. This independently exercises cross-covariance and coefficients.
set.seed(92)
B <- matrix(rnorm(49), 7, 7)
V <- tcrossprod(B) + diag(7)
Ve <- V[1:5, 1:5]
C <- V[6:7, 1:5]
Vn <- V[6:7, 6:7]
Km <- solve(Ve)
fm <- pip_uncertainty_fit(X, Km, r)
fm$sigma2 <- 1
fm$beta_vcov <- solve(crossprod(X, Km %*% X))
Sm <- pip_prediction_covariance(fm, Xn, C, Vn)
errors <- t(chol(V)) %*% matrix(rnorm(7 * 60000), 7)
beta_hat <- solve(crossprod(X, Km %*% X), crossprod(X, Km %*% errors[1:5, ]))
pred <- Xn %*% beta_hat + C %*% Km %*% (errors[1:5, ] - X %*% beta_hat)
pe <- errors[6:7, ] - pred
stopifnot(max(abs(cov(t(pe)) - Sm)) / max(abs(Sm)) < 0.025,
          abs(var(colMeans(pe)) / (sum(Sm) / 4) - 1) < 0.025)

# Permuting training and prediction rows preserves the same joint errors.
p <- c(5, 3, 1, 4, 2)
fp <- pip_uncertainty_fit(X[p, ], Km[p, p], r[p])
f0 <- pip_uncertainty_fit(X, Km, r)
S0 <- pip_prediction_covariance(f0, Xn, C, Vn)
Sp <- pip_prediction_covariance(fp, Xn[2:1, ], C[2:1, p], Vn[2:1, 2:1])
stopifnot(isTRUE(all.equal(Sp, S0[2:1, 2:1], tolerance = 1e-8)))

# Lambda must leave the diagonal on the original tree scale.
Vl <- pip_lambda_covariance(V, 0.4)
stopifnot(identical(diag(Vl), diag(V)), Vl[1, 2] == 0.4 * V[1, 2])
bad <- try(pip_prediction_covariance(fit, Xn, matrix(100, 2, 5), diag(2)), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
# Site discrepancy: E[r^2] = tau^2 + E[se^2], floored at zero.
stopifnot(isTRUE(all.equal(pip_site_discrepancy_var(c(3, -4), c(1, 1)), 11.5)),
          pip_site_discrepancy_var(c(0.1, -0.1), c(1, 1)) == 0)
ci <- pip_calibrated_interval(c(10, 20), c(3, 0), 16)
stopifnot(isTRUE(all.equal(ci$total_se, c(5, 4))),
          isTRUE(all.equal(ci$calibrated_upper - c(10, 20), qnorm(0.975) * c(5, 4))),
          isTRUE(all.equal(c(10, 20) - ci$calibrated_lower, qnorm(0.975) * c(5, 4))))
cat("PIP covariance analytical and simulation checks passed.\n")
