# Joint trait-kernel and phylogenetic covariance model for log(MAP).
# All B blocks must come from the same rooted scaffold and use one root height.

map_kernel_dist2 <- function(x, z = x) {
  x <- as.matrix(x); z <- as.matrix(z)
  stopifnot(ncol(x) == ncol(z), ncol(x) > 0L)
  d <- outer(rowSums(x^2), rowSums(z^2), "+") - 2 * tcrossprod(x, z)
  pmax(d, 0)
}

map_kernel_reference <- function(x, dist2 = NULL) {
  d <- if (is.null(dist2)) map_kernel_dist2(x) else dist2
  v <- sqrt(d[upper.tri(d) & d > 1e-12])
  if (!length(v)) stop("No positive training trait distances")
  median(v)
}

map_kernel_matrix <- function(x, z, kind, ell = NA_real_) {
  if (kind == "linear") return(tcrossprod(x, z) / ncol(x))
  if (kind == "rbf") {
    if (!is.finite(ell) || ell <= 0) stop("RBF bandwidth must be positive")
    return(exp(-map_kernel_dist2(x, z) / (2 * ell^2)))
  }
  stop("Unknown trait kernel: ", kind)
}

map_kernel_design <- function(z, mean_kind, retained = NULL) {
  h <- if (mean_kind == "linear") cbind(`(Intercept)` = 1, z) else
    matrix(1, nrow(z), 1L, dimnames = list(NULL, "(Intercept)"))
  if (!is.null(retained)) return(h[, retained, drop = FALSE])
  keep <- 1L
  if (ncol(h) > 1L) for (j in 2:ncol(h)) {
    if (qr(h[, c(keep, j), drop = FALSE], tol = 1e-10)$rank > length(keep))
      keep <- c(keep, j)
  }
  h[, keep, drop = FALSE]
}

map_kernel_chol <- function(a) {
  scale <- max(1, max(abs(diag(a))))
  for (relative in c(0, 1e-12, 1e-10, 1e-8, 1e-6)) {
    jitter <- relative * scale
    u <- tryCatch(chol(a + diag(jitter, nrow(a))), error = function(e) NULL)
    if (!is.null(u)) return(list(u = u, jitter = jitter))
  }
  stop("Joint covariance is not positive definite within numerical jitter limit")
}

map_kernel_solve <- function(u, b) {
  backsolve(u, forwardsolve(t(u), b))
}

map_kernel_fit <- function(y, z, b, lambda_phy, eta, kind,
                           mean_kind = "intercept", ell = NA_real_,
                           k_train = NULL) {
  y <- as.numeric(y); z <- as.matrix(z); b <- as.matrix(b)
  stopifnot(length(y) == nrow(z), identical(dim(b), c(length(y), length(y))),
            all(is.finite(y)), all(is.finite(z)), all(is.finite(b)),
            lambda_phy >= 0, lambda_phy <= 1, eta >= 0)
  h <- map_kernel_design(z, mean_kind)
  k <- if (eta == 0) matrix(0, length(y), length(y)) else if
    (!is.null(k_train)) as.matrix(k_train) else
      map_kernel_matrix(z, z, kind, ell)
  stopifnot(identical(dim(k), c(length(y), length(y))))
  v <- lambda_phy * b
  diag(v) <- diag(b)
  ch <- map_kernel_chol(v + eta * k)
  ah <- map_kernel_solve(ch$u, h)
  ay <- map_kernel_solve(ch$u, y)
  beta <- as.numeric(solve(crossprod(h, ah), crossprod(h, ay)))
  names(beta) <- colnames(h)
  alpha <- as.numeric(map_kernel_solve(ch$u, y - as.numeric(h %*% beta)))
  list(beta = beta, alpha = alpha, z = z, lambda_phy = lambda_phy,
       eta = eta, kind = kind, ell = ell, mean_kind = mean_kind,
       retained = colnames(h), removed = setdiff(colnames(
         if (mean_kind == "linear") cbind(`(Intercept)` = 1, z) else h),
         colnames(h)), jitter = ch$jitter)
}

map_kernel_predict <- function(fit, new_z, b_cross, k_cross = NULL) {
  new_z <- as.matrix(new_z); b_cross <- as.matrix(b_cross)
  stopifnot(ncol(new_z) == ncol(fit$z), nrow(b_cross) == nrow(new_z),
            ncol(b_cross) == nrow(fit$z), all(is.finite(b_cross)))
  h <- map_kernel_design(new_z, fit$mean_kind, fit$retained)
  mean <- as.numeric(h %*% fit$beta)
  if (!is.null(k_cross))
    stopifnot(identical(dim(k_cross), c(nrow(new_z), nrow(fit$z))))
  trait <- if (fit$eta == 0) rep(0, nrow(new_z)) else
    as.numeric(fit$eta * (if (is.null(k_cross))
      map_kernel_matrix(new_z, fit$z, fit$kind, fit$ell) else
        k_cross) %*% fit$alpha)
  phylo <- as.numeric(fit$lambda_phy * b_cross %*% fit$alpha)
  data.frame(linear_mean = mean, trait_adjustment = trait,
             phylo_adjustment = phylo, prediction = mean + trait + phylo)
}

map_kernel_scale_fit <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x); spread <- apply(x, 2, sd)
  keep <- is.finite(spread) & spread > 1e-10
  if (!any(keep)) stop("No usable nonconstant training predictors")
  list(center = center[keep], spread = spread[keep],
       retained = colnames(x)[keep], removed = colnames(x)[!keep])
}

map_kernel_scale_apply <- function(x, scaling) {
  z <- sweep(sweep(as.matrix(x[, scaling$retained, drop = FALSE]), 2,
                   scaling$center, "-"), 2, scaling$spread, "/")
  if (any(!is.finite(z))) stop("Nonfinite completed/scaled traits")
  z
}
