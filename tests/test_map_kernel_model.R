source("code/map_kernel_model.R")
stopifnot(requireNamespace("ape", quietly = TRUE))
set.seed(510)
tree <- ape::rtree(48)
b_all <- ape::vcv(tree)
b_all <- b_all / max(ape::node.depth.edgelength(tree))
ids <- tree$tip.label
train <- ids[1:36]
query <- ids[37:48]
x <- setNames(seq(-2.3, 2.3, length.out = length(ids)), sample(ids))
z <- matrix(x[ids], ncol = 1, dimnames = list(ids, "leaf_trait"))
b <- b_all[train, train, drop = FALSE]
cross <- b_all[query, train, drop = FALSE]
noise <- setNames(as.numeric(t(chol(b_all + diag(1e-8, nrow(b_all)))) %*%
                             rnorm(length(ids))), ids) * .04
y_curve <- 1.5 * sin(2 * z[, 1]) + noise[ids]
y_line <- 1 + .8 * z[, 1] + noise[ids]

# eta=0 is the fixed-lambda PIP GLS equation, including its cross block.
fit <- map_kernel_fit(y_line[train], z[train, , drop = FALSE], b,
                      lambda_phy = .5, eta = 0, kind = "rbf",
                      mean_kind = "linear")
pred <- map_kernel_predict(fit, z[query, , drop = FALSE], cross)
h <- cbind(1, z[train, , drop = FALSE])
hn <- cbind(1, z[query, , drop = FALSE])
v <- .5 * b; diag(v) <- diag(b)
beta <- solve(t(h) %*% solve(v, h), t(h) %*% solve(v, y_line[train]))
manual <- as.numeric(hn %*% beta + .5 * cross %*%
                       solve(v, y_line[train] - h %*% beta))
stopifnot(max(abs(pred$prediction - manual)) < 1e-8,
          max(abs(pred$trait_adjustment)) == 0)

# A zero phylogenetic lambda has no prediction-time phylogenetic contribution.
zero <- map_kernel_fit(y_curve[train], z[train, , drop = FALSE], b,
                       lambda_phy = 0, eta = 1, kind = "rbf",
                       ell = 1, mean_kind = "intercept")
stopifnot(max(abs(map_kernel_predict(zero, z[query, , drop = FALSE], cross)$
                    phylo_adjustment)) == 0)

# Calibration and query order, and query batch composition, cannot change results.
permutation <- sample(seq_along(train))
reordered <- map_kernel_fit(y_line[train][permutation],
  z[train[permutation], , drop = FALSE], b[permutation, permutation],
  lambda_phy = .5, eta = 0, kind = "rbf", mean_kind = "linear")
stopifnot(max(abs(pred$prediction - map_kernel_predict(reordered,
  z[query, , drop = FALSE], cross[, permutation])$prediction)) < 1e-8)
query_order <- sample(seq_along(query))
stopifnot(max(abs(pred$prediction[query_order] - map_kernel_predict(fit,
  z[query[query_order], , drop = FALSE], cross[query_order, ])$prediction)) < 1e-8)
stopifnot(abs(pred$prediction[1] - map_kernel_predict(fit,
  z[query[1], , drop = FALSE], cross[1, , drop = FALSE])$prediction) < 1e-8)

# The numerical jitter is absent from the biological joint covariance.
v_all <- .7 * b_all; diag(v_all) <- diag(b_all)
k_all <- map_kernel_matrix(z, z, "rbf", ell = 1)
stopifnot(min(eigen(v_all + k_all, symmetric = TRUE,
                    only.values = TRUE)$values) > -1e-8)

# A curved truth should favour RBF over a linear trait kernel on this
# deterministic small example. A linear truth is recovered by the linear mean.
rbf <- map_kernel_fit(y_curve[train], z[train, , drop = FALSE], b,
                      .5, 1, "rbf", "intercept", 1)
cached <- map_kernel_fit(y_curve[train], z[train, , drop = FALSE], b,
                         .5, 1, "rbf", "intercept", 1,
                         k_train = map_kernel_matrix(z[train, , drop = FALSE],
                                                     z[train, , drop = FALSE],
                                                     "rbf", 1))
stopifnot(max(abs(map_kernel_predict(rbf, z[query, , drop = FALSE], cross)$prediction -
  map_kernel_predict(cached, z[query, , drop = FALSE], cross,
    k_cross = map_kernel_matrix(z[query, , drop = FALSE],
                                z[train, , drop = FALSE], "rbf", 1))$prediction)) < 1e-10)
linear_kernel <- map_kernel_fit(y_curve[train], z[train, , drop = FALSE], b,
                               .5, 1, "linear", "intercept")
e_rbf <- mean((map_kernel_predict(rbf, z[query, , drop = FALSE], cross)$prediction -
                 y_curve[query])^2)
e_linear <- mean((map_kernel_predict(linear_kernel, z[query, , drop = FALSE],
                                    cross)$prediction - y_curve[query])^2)
stopifnot(e_rbf < e_linear,
          mean((pred$prediction - y_line[query])^2) < .01)

# Duplicate mean columns are removed by a deterministic training-only rule.
duplicate <- cbind(z[train, , drop = FALSE], copy = z[train, 1])
rank_fit <- map_kernel_fit(y_line[train], duplicate, b, .5, 0,
                           "linear", "linear")
stopifnot(identical(rank_fit$removed, "copy"))
cat("map kernel model checks passed\n")
