# Independent reconstruction of every CV interval from saved fits and raw
# held-out traits. Does not call pip_uncertainty_fit/pip_prediction_covariance.
# Uses the production preprocessing, but independent GLS algebra for variance.
# Run from repository root; optionally set DILP_SOURCE to a local checkout.
if (nzchar(Sys.getenv("DILP_SOURCE"))) pkgload::load_all(Sys.getenv("DILP_SOURCE"), quiet = TRUE)
# Evaluate only setup/data preparation, before the first output-writing step.
# Stop by expression identity rather than fragile source line numbers.
expressions <- parse("code/03_loso_cv.R")
found_stop <- FALSE
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("dilp_cv"))) {
    found_stop <- TRUE
    break
  }
  eval(expr, envir = globalenv())
}
stopifnot(found_stop)
saved <- read.csv("tables/pip_cv_site_uncertainty.csv")
raw_truth <- aggregate(raw_dat[, c("mat", "map")], list(site = raw_dat$Site), mean, na.rm = TRUE)
raw_truth$log_map <- log(raw_truth$map)
rownames(raw_truth) <- raw_truth$site
rows <- list()
for (fold in seq_len(K_FOLDS)) {
  archive <- readRDS(sprintf("models/loso_cv_fold_%02d.rds", fold))
  stopifnot(setequal(archive$held_sites, names(fold_assignment)[fold_assignment == fold]),
            !length(intersect(archive$held_sites, archive$train_sites)))
  train <- agg_species(fill_tooth_traits(raw_dat[raw_dat$Site %in% archive$train_sites, ]))
  train <- train[intersect(rownames(train), rownames(full_vcv)), ]
  predictors <- active_pred_names(train, fossil_traits)
  stopifnot(identical(predictors, archive$pred_names))
  set.seed(SEED + fold)
  imputer <- preProcess(train[, predictors, drop = FALSE], method = "bagImpute")
  for (target in c("mat", "log_map")) {
    pc <- archive$pgls_res$impute[[target]]
    ids_train <- pc$sp_fit
    # Reconstruct covariance from scaffold and fitted lambda, independently
    # verify against the fitted model's matrix and saved inverse.
    V <- full_vcv[ids_train, ids_train] * pc$lambda
    diag(V) <- diag(full_vcv[ids_train, ids_train])
    stopifnot(max(abs(V - pc$pgls_fit$V[ids_train, ids_train])) < 1e-8)
    K <- chol2inv(chol(V))
    stopifnot(max(abs(K - pc$K_train)) < 1e-8)
    X <- pc$pgls_fit$x[ids_train, pc$common_vars, drop = FALSE]
    Y <- setNames(pc$pgls_fit$y, rownames(pc$pgls_fit$x))[ids_train]
    B <- solve(crossprod(X, K %*% X), crossprod(X, K))
    beta <- B %*% Y
    residual <- Y - as.vector(X %*% beta)
    sigma2 <- as.numeric(crossprod(residual, K %*% residual)) / (nrow(X) - ncol(X))
    stopifnot(abs(sigma2 - as.numeric(pc$pgls_fit$RMS)) < 1e-8)
    for (site in archive$held_sites) {
      leaves <- fill_tooth_traits(raw_dat[raw_dat$Site == site, ])
      local <- aggregate(leaves[, predictors, drop = FALSE],
                         list(species = leaves$genusSpecies), mean, na.rm = TRUE)
      rownames(local) <- local$species
      local <- nan_to_na(local)
      newdata <- predict(imputer, local[, predictors, drop = FALSE])
      rownames(newdata) <- local$species
      Xf <- model.matrix(delete.response(terms(pc$formula)), newdata)[, colnames(X), drop = FALSE]
      ids <- rownames(Xf)
      present <- intersect(ids, rownames(full_vcv))
      n <- length(ids)
      C <- matrix(0, n, length(ids_train), dimnames = list(ids, ids_train))
      C[present, ] <- pc$lambda * full_vcv[present, ids_train, drop = FALSE]
      Vf <- diag(max(diag(full_vcv)), n)
      dimnames(Vf) <- list(ids, ids)
      Vf[present, present] <- full_vcv[present, present, drop = FALSE] * pc$lambda
      diag(Vf)[match(present, ids)] <- diag(full_vcv[present, present, drop = FALSE])
      # Prediction is a linear function H Y of training responses. Therefore
      # prediction error variance is obtained directly from Var(Y), Var(Yf),
      # and Cov(Y,Yf), with no conditional-covariance helper or stored beta VC.
      H <- C %*% K + (Xf - C %*% K %*% X) %*% B
      stopifnot(max(abs(H %*% X - Xf)) < 1e-7)
      h <- colMeans(H)
      v_future <- mean(Vf)
      v_predictor <- as.numeric(crossprod(h, V %*% h))
      c_joint <- as.numeric(crossprod(h, colMeans(C)))
      variance <- sigma2 * (v_future + v_predictor - 2 * c_joint)
      stopifnot(variance > 0)
      estimate <- sum(h * Y)
      half <- qt(.975, nrow(X) - ncol(X)) * sqrt(variance)
      truth <- raw_truth[site, target]
      record <- saved[saved$site == site & saved$target == target, ]
      stopifnot(nrow(record) == 1L, record$n_species == n,
                abs(record$observed - truth) < 1e-7,
                abs(record$estimate - estimate) < 1e-7,
                abs(record$se - sqrt(variance)) < 1e-7,
                abs(record$lower - (estimate-half)) < 1e-7,
                abs(record$upper - (estimate+half)) < 1e-7)
      rows[[paste(fold, site, target)]] <- data.frame(site, fold, target,
        n_species = n, observed = truth, estimate, se = sqrt(variance),
        lower = estimate-half, upper = estimate+half,
        covered = abs(truth-estimate) <= half,
        se_difference = sqrt(variance)-record$se,
        estimate_difference = estimate-record$estimate)
    }
  }
  cat("Independently reconstructed fold", fold, "\n")
  rm(archive); gc()
}
audit <- do.call(rbind, rows)
stopifnot(nrow(audit) == 184L)
write.csv(audit, "doc/audits/cv_coverage_independent.csv", row.names = FALSE)
print(aggregate(covered ~ target, audit, function(x) c(covered=sum(x), total=length(x))))
cat("Maximum SE difference:", max(abs(audit$se_difference)), "\n")
cat("Maximum prediction difference:", max(abs(audit$estimate_difference)), "\n")
