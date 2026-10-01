# Real-data eta=0 parity using the production PIP's response-inclusive
# training imputation and trait-only held-out imputation, at its fitted lambda.
if (nzchar(Sys.getenv("DILP_SOURCE")))
  pkgload::load_all(Sys.getenv("DILP_SOURCE"), quiet = TRUE)
Sys.setenv(PIP_CV_SCHEME = "ten_fold")
found_stop <- FALSE
for (expr in parse("code/03_loso_cv.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("dilp_cv"))) {
    found_stop <- TRUE
    break
  }
  eval(expr, envir = globalenv())
}
stopifnot(found_stop)
source("code/map_kernel_model.R")
archive <- readRDS("models/loso_cv_fold_01.rds")
pc <- archive$pgls_res$impute$log_map
pn <- pc$pred_names
filled <- fill_tooth_traits(raw_dat[raw_dat$Site %in% archive$train_sites, ])
sp <- agg_species(filled)
sp <- sp[intersect(rownames(sp), rownames(full_vcv)), ]
train_raw <- sp[, c("log_map", pn), drop = FALSE]
set.seed(43)
train_imp <- preProcess(train_raw, method = "bagImpute")
completed <- predict(train_imp, train_raw)[pc$sp_fit, ]
set.seed(43)
query_imp <- preProcess(sp[, pn, drop = FALSE], method = "bagImpute")
trait_only_train <- predict(query_imp, sp[, pn, drop = FALSE])
imputation_delta <- abs(as.matrix(trait_only_train[pc$sp_fit, pn]) -
                        as.matrix(completed[, pn]))
held <- fill_tooth_traits(raw_dat[raw_dat$Site %in% archive$held_sites, ])
occ <- aggregate(held[, pn, drop = FALSE],
                 list(site = held$Site, species = held$genusSpecies),
                 mean, na.rm = TRUE)
occ <- nan_to_na(occ)
new_x <- predict(query_imp, occ[, pn, drop = FALSE])
ids <- pc$sp_fit
root_height <- max(phytools::nodeHeights(read.tree("data/tre_scaffold.tre")))
b <- full_vcv[ids, ids, drop = FALSE]
diag(b) <- diag(b) - 1e-6
b <- b / root_height
cross <- matrix(0, nrow(occ), length(ids))
represented <- occ$species %in% rownames(full_vcv)
cross[represented, ] <- full_vcv[occ$species[represented], ids,
                                 drop = FALSE] / root_height
fit <- map_kernel_fit(completed$log_map, as.matrix(completed[, pn]), b,
                      lambda_phy = pc$lambda, eta = 0, kind = "linear",
                      mean_kind = "linear")
pred <- map_kernel_predict(fit, as.matrix(new_x), cross)$prediction
site <- tapply(pred, occ$site, mean)
published <- read.csv("tables/pip_cv_site_uncertainty.csv")
published <- published[published$fold == 1 & published$target == "log_map", ]
stopifnot(setequal(names(site), published$site),
          max(abs(site[published$site] - published$estimate)) < 1e-6)

# A real training/query joint block uses one root scale and is PSD before
# numerical jitter, including a held occurrence of a represented species.
pick_train <- ids[1:30]
pick_query <- which(represented)[1:5]
joint_ids <- c(pick_train, occ$species[pick_query])
b_base <- full_vcv
diag(b_base) <- diag(b_base) - 1e-6
b_joint <- b_base[joint_ids, joint_ids, drop = FALSE]
b_joint <- b_joint / root_height
v_joint <- .5 * b_joint
diag(v_joint) <- diag(b_joint)
z_joint <- rbind(as.matrix(completed[pick_train, pn]),
                 as.matrix(new_x[pick_query, pn]))
sc <- map_kernel_scale_fit(z_joint)
z_joint <- map_kernel_scale_apply(z_joint, sc)
k_joint <- map_kernel_matrix(z_joint, z_joint, "rbf", ell = 1)
stopifnot(min(eigen(v_joint + k_joint, symmetric = TRUE,
                    only.values = TRUE)$values) > -1e-8)

# Changing only held-out responses cannot affect fitting or predictions.
changed <- raw_dat
changed$map[changed$Site %in% archive$held_sites] <- 1e9
changed$mat[changed$Site %in% archive$held_sites] <- -1e9
held_changed <- fill_tooth_traits(changed[changed$Site %in% archive$held_sites, ])
occ_changed <- aggregate(held_changed[, pn, drop = FALSE],
                         list(site = held_changed$Site,
                              species = held_changed$genusSpecies),
                         mean, na.rm = TRUE)
occ_changed <- nan_to_na(occ_changed)
stopifnot(identical(occ, occ_changed))
cat("Real-data PIP limit and held-response exclusion checks passed\n")
cat("Training trait cells changed by trait-only imputation:",
    sum(imputation_delta > 1e-10), "/", length(imputation_delta),
    "maximum absolute change:", max(imputation_delta), "\n")
