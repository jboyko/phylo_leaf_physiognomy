# Refit the published DiLP equation structures on whole-site CV folds.
# site_traits must come from dilp()$processed_site_data: in particular,
# log tooth traits are transformed before site averaging, not afterwards.
# cv supplies the exact site folds and observed responses used by other models.
dilp_site_cv <- function(site_traits, cv) {
  stopifnot(!anyDuplicated(site_traits$site), !anyDuplicated(cv$site),
            all(cv$site %in% site_traits$site), !anyNA(cv$fold))
  d <- site_traits[match(cv$site, site_traits$site), , drop = FALSE]
  d$mat <- cv$obs_mat
  d$log_map <- cv$obs_log_map
  formulas <- list(mat = mat ~ margin + fdr + tc_ip,
                   log_map = log_map ~ ln_leaf_area + ln_tc_ip + ln_pr)
  predictions <- cv[, c("site", "fold", "obs_mat", "obs_log_map")]
  coefficients <- list()
  for (target in names(formulas)) {
    formula <- formulas[[target]]
    vars <- all.vars(formula)
    finite <- apply(d[, vars, drop = FALSE], 1, function(x) all(is.finite(x)))
    pred <- rep(NA_real_, nrow(d))
    for (fold in sort(unique(cv$fold))) {
      train <- cv$fold != fold & finite
      held <- cv$fold == fold
      predictors_ok <- apply(d[, vars[-1], drop = FALSE], 1,
                             function(x) all(is.finite(x)))
      stopifnot(sum(train) > length(vars))
      fit <- lm(formula, data = d[train, , drop = FALSE])
      stopifnot(!anyNA(coef(fit)))
      usable <- held & predictors_ok
      pred[usable] <- predict(fit, newdata = d[usable, , drop = FALSE])
      coefficients[[paste(target, fold)]] <- data.frame(
        target = target, fold = fold, n_train = sum(train),
        predictor = names(coef(fit)), estimate = unname(coef(fit)))
    }
    predictions[[paste0("dilp_cv_site_", target)]] <- pred
  }
  list(predictions = predictions, coefficients = do.call(rbind, coefficients))
}
