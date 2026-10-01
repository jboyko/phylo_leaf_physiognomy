# End-to-end invariant: the extant covariance used to fit climate must equal
# the extant block on the scaffold used to graft fossils (apart from jitter).
library(ape)
source("code/pip_uncertainty.R")
p <- readRDS("models/pip_components.rds")
V <- vcv(read.tree("data/tre_scaffold.tre"))
for (cfg in names(p$configs)) for (target in c("mat", "log_map")) {
  fit <- p$configs[[cfg]][[target]]
  ids <- rownames(fit$X)
  expected <- V[ids, ids]
  diag(expected) <- diag(expected) + 1e-6
  expected <- pip_lambda_covariance(expected, fit$lambda)
  stopifnot(max(abs(expected - fit$V_lam[ids, ids])) < 1e-7)
}
cat("Climate fit and fossil scaffold covariance roots agree.\n")
