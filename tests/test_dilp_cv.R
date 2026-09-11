source("code/dilp_cv.R")
set.seed(42)
n <- 40
d <- data.frame(site = paste0("site", seq_len(n)),
                margin = runif(n), fdr = runif(n), tc_ip = runif(n),
                ln_leaf_area = rnorm(n), ln_tc_ip = rnorm(n), ln_pr = rnorm(n))
cv <- data.frame(site = d$site, fold = rep(1:10, 4),
                 obs_mat = rnorm(n), obs_log_map = rnorm(n))
base <- dilp_site_cv(d, cv)
# Held-out outcomes must have no influence on their own predictions.
altered <- cv
altered$obs_mat[cv$fold == 1] <- altered$obs_mat[cv$fold == 1] + 100
altered$obs_log_map[cv$fold == 1] <- altered$obs_log_map[cv$fold == 1] - 100
changed <- dilp_site_cv(d, altered)
cols <- c("dilp_cv_site_mat", "dilp_cv_site_log_map")
stopifnot(identical(base$predictions[cv$fold == 1, cols],
                    changed$predictions[cv$fold == 1, cols]))
# Independent fold fit agrees; row order of trait input does not matter.
fit <- lm(obs_mat ~ margin + fdr + tc_ip,
          data = cbind(cv, d[, -1])[cv$fold != 1, ])
stopifnot(isTRUE(all.equal(unname(predict(fit, d[cv$fold == 1, ])),
  base$predictions$dilp_cv_site_mat[cv$fold == 1])))
stopifnot(identical(base, dilp_site_cv(d[n:1, ], cv)))
# Non-finite logged tooth traits exclude MAP predictions, not MAT predictions.
d$ln_tc_ip[1] <- -Inf
missing <- dilp_site_cv(d, cv)
stopifnot(is.na(missing$predictions$dilp_cv_site_log_map[1]),
          is.finite(missing$predictions$dilp_cv_site_mat[1]))
cat("DiLP fold isolation, regression, alignment and missing-data checks passed.\n")
