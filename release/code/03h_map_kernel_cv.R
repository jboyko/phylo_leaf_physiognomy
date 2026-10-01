# Nested, site-grouped log(MAP) kernel/phylogeny point-prediction experiment.
# Keeps production climate models and fossil tables untouched.
if (nzchar(Sys.getenv("DILP_SOURCE")))
  pkgload::load_all(Sys.getenv("DILP_SOURCE"), quiet = TRUE)
Sys.setenv(PIP_CV_SCHEME = "ten_fold")
# Reuse the exact calibration preprocessing and fold definitions, not the CV loop.
found_stop <- FALSE
for (expr in parse("code/03_loso_cv.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("dilp_cv"))) {
    found_stop <- TRUE
    break
  }
  eval(expr, envir = globalenv())
}
stopifnot(found_stop, length(all_sites) == 92L)
source("code/map_kernel_model.R")
out_dir <- "models/map_kernel_experiment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 03_loso_cv adds 1e-6 for its own numerical solve. The biological Brownian
# covariance here excludes that jitter and uses one full-scaffold root scale.
scaffold <- read.tree("data/tre_scaffold.tre")
root_height <- max(phytools::nodeHeights(scaffold))
stopifnot(is.finite(root_height), root_height > 0,
          all(diag(full_vcv) > 1e-6))
b_scaffold <- full_vcv
diag(b_scaffold) <- diag(b_scaffold) - 1e-6
b_scaffold <- b_scaffold / root_height

prepare_kernel_split <- function(train_sites, test_sites, seed) {
  stopifnot(!length(intersect(train_sites, test_sites)),
            all(train_sites %in% all_sites), all(test_sites %in% all_sites))
  filled <- fill_tooth_traits(raw_dat[raw_dat$Site %in% train_sites, ])
  sp <- agg_species(filled)
  sp <- sp[intersect(rownames(sp), rownames(b_scaffold)), , drop = FALSE]
  pn <- active_pred_names(sp, fossil_traits)
  if (!length(pn)) stop("No eligible fossil-measurable predictors")
  held <- fill_tooth_traits(raw_dat[raw_dat$Site %in% test_sites, ])
  occ <- aggregate(held[, pn, drop = FALSE],
                   list(site = held$Site, species = held$genusSpecies),
                   mean, na.rm = TRUE)
  occ <- nan_to_na(occ)
  sp <- nan_to_na(sp)
  set.seed(seed)
  imp <- preProcess(sp[, pn, drop = FALSE], method = "bagImpute")
  x <- predict(imp, sp[, pn, drop = FALSE])
  new_x <- predict(imp, occ[, pn, drop = FALSE])
  assert_finite(x, "kernel training traits")
  assert_finite(new_x, "kernel query traits")
  scaling <- map_kernel_scale_fit(x)
  z <- map_kernel_scale_apply(x, scaling)
  new_z <- map_kernel_scale_apply(new_x, scaling)
  colnames(z) <- colnames(new_z) <- scaling$retained
  ids <- rownames(sp)
  b <- b_scaffold[ids, ids, drop = FALSE]
  cross <- matrix(0, nrow(occ), length(ids))
  represented <- occ$species %in% rownames(b_scaffold)
  cross[represented, ] <- b_scaffold[occ$species[represented], ids,
                                     drop = FALSE]
  dist_train <- map_kernel_dist2(z)
  dist_cross <- map_kernel_dist2(new_z, z)
  list(y = sp$log_map, z = z, new_z = new_z, b = b, cross = cross,
       occ = occ[, c("site", "species")], ids = ids, pn = pn,
       imputer = imp, scaling = scaling,
       dist_train = dist_train, dist_cross = dist_cross,
       linear_train = tcrossprod(z) / ncol(z),
       linear_cross = tcrossprod(new_z, z) / ncol(z),
       reference = map_kernel_reference(z, dist_train),
       train_sites = train_sites, test_sites = test_sites, seed = seed,
       root_height = root_height)
}

grid <- do.call(rbind, list(
  expand.grid(candidate = "linear_kernel", mean_kind = "intercept",
              kind = "linear", lambda_phy = c(0, .5, .9, .99),
              eta = c(.1, 1, 10), ell_factor = NA_real_),
  expand.grid(candidate = "rbf_kernel", mean_kind = "intercept",
              kind = "rbf", lambda_phy = c(0, .5, .9, .99),
              eta = c(.1, 1, 10), ell_factor = c(.5, 1, 2)),
  expand.grid(candidate = "linear_mean_rbf", mean_kind = "linear",
              kind = "rbf", lambda_phy = c(0, .5, .9, .99),
              eta = c(.1, 1, 10), ell_factor = c(.5, 1, 2)),
  expand.grid(candidate = "linear_mean_phylogeny", mean_kind = "linear",
              kind = "rbf", lambda_phy = c(0, .5, .9, .99),
              eta = 0, ell_factor = NA_real_)
))
grid$config_id <- seq_len(nrow(grid))
grid$boundary <- with(grid, lambda_phy %in% c(0, .99) |
                        eta %in% c(.1, 10) |
                        (!is.na(ell_factor) & ell_factor %in% c(.5, 2)))
candidate_names <- unique(grid$candidate)

predict_config <- function(d, cfg) {
  ell <- if (is.na(cfg$ell_factor)) NA_real_ else cfg$ell_factor * d$reference
  k_train <- k_cross <- NULL
  if (cfg$eta > 0) {
    if (cfg$kind == "linear") {
      k_train <- d$linear_train
      k_cross <- d$linear_cross
    } else {
      k_train <- exp(-d$dist_train / (2 * ell^2))
      k_cross <- exp(-d$dist_cross / (2 * ell^2))
    }
  }
  fit <- map_kernel_fit(d$y, d$z, d$b, cfg$lambda_phy, cfg$eta,
                        cfg$kind, cfg$mean_kind, ell, k_train)
  components <- map_kernel_predict(fit, d$new_z, d$cross, k_cross)
  cbind(d$occ, components, map_cm = exp(components$prediction))
}

run_outer <- function(outer) {
  started <- proc.time()["elapsed"]
  archive <- readRDS(sprintf("models/map_model_experiment/outer_%02d.rds", outer))
  train_sites <- archive$train_sites
  held_sites <- archive$held_sites
  stopifnot(setequal(held_sites, names(fold_assignment)[fold_assignment == outer]))
  inner_membership <- archive$inner_membership
  losses <- matrix(NA_real_, nrow(grid), length(inner_membership))
  inner_records <- vector("list", length(inner_membership))
  for (j in seq_along(inner_membership)) {
    m <- inner_membership[[j]]
    stopifnot(!length(intersect(c(m$train_sites, m$held_sites), held_sites)),
              setequal(c(m$train_sites, m$held_sites), train_sites))
    d <- prepare_kernel_split(m$train_sites, m$held_sites, m$seed)
    truth <- dat_site_obs[m$held_sites, "log_map"]
    names(truth) <- m$held_sites
    for (i in seq_len(nrow(grid))) {
      p <- predict_config(d, grid[i, ])
      site <- tapply(p$prediction, p$site, mean)
      losses[i, j] <- sum((site - truth[names(site)])^2)
    }
    inner_records[[j]] <- list(train_sites = m$train_sites,
                               held_sites = m$held_sites, seed = m$seed,
                               predictors = d$pn,
                               scaling = d$scaling,
                               n_species = length(d$ids))
    cat("outer", outer, "inner", j, "complete\n")
  }
  # Each inner site is held out exactly once: pooled site MSE weights sites equally.
  counts <- vapply(inner_membership, function(m) length(m$held_sites), integer(1))
  stopifnot(setequal(unlist(lapply(inner_membership, `[[`, "held_sites")), train_sites),
            !anyDuplicated(unlist(lapply(inner_membership, `[[`, "held_sites"))))
  inner_scores <- cbind(grid, mse = rowSums(losses) / sum(counts))
  selected <- do.call(rbind, lapply(candidate_names, function(nm) {
    subset <- inner_scores[inner_scores$candidate == nm, ]
    subset[which.min(subset$mse), , drop = FALSE]
  }))
  d <- prepare_kernel_split(train_sites, held_sites, 42 + outer)
  occurrence <- do.call(rbind, lapply(seq_len(nrow(selected)), function(k) {
    cfg <- selected[k, ]
    p <- predict_config(d, cfg)
    cbind(fold = outer, candidate = cfg$candidate, p)
  }))
  rownames(occurrence) <- NULL
  site <- aggregate(occurrence[, c("linear_mean", "trait_adjustment",
                                   "phylo_adjustment", "prediction")],
                    list(fold = occurrence$fold,
                         candidate = occurrence$candidate,
                         site = occurrence$site), mean)
  site$observed <- dat_site_obs[site$site, "log_map"]
  site$map_cm <- exp(site$prediction)
  old <- archive$predictions
  comparator <- rbind(
    data.frame(fold = outer, candidate = "pip", site = old$site,
               prediction = old$pip, observed = old$observed,
               map_cm = exp(old$pip)),
    data.frame(fold = outer, candidate = "pip_recalibrated", site = old$site,
               prediction = old$pip_recalibrated, observed = old$observed,
               map_cm = exp(old$pip_recalibrated)))
  stopifnot(setequal(site$site, held_sites),
            max(abs(site$prediction - site$linear_mean -
                    site$trait_adjustment - site$phylo_adjustment)) < 1e-9,
            all(is.finite(occurrence$prediction)),
            all(is.finite(site$prediction)))
  checkpoint <- list(outer = outer, occurrence = occurrence, site = site,
                     comparator = comparator, selected = selected,
                     inner_scores = inner_scores, inner_losses = losses,
                     inner_records = inner_records,
                     inner_membership = inner_membership,
                     train_sites = train_sites, held_sites = held_sites,
                     outer_preprocessing = list(predictors = d$pn,
                       scaling = d$scaling, n_species = length(d$ids),
                       seed = d$seed), grid = grid,
                     root_height = root_height, session = sessionInfo())
  saveRDS(checkpoint, file.path(out_dir, sprintf("outer_%02d.rds", outer)))
  cat("OUTER", outer, "DONE in",
      round(proc.time()["elapsed"] - started, 1), "seconds\n")
  invisible(outer)
}

if (Sys.getenv("PIP_KERNEL_SETUP_ONLY", "0") != "1") {
  fold_arg <- Sys.getenv("PIP_KERNEL_FOLDS", "")
  folds <- if (nzchar(fold_arg))
    as.integer(strsplit(fold_arg, ",", fixed = TRUE)[[1]]) else 1:10
  stopifnot(length(folds) > 0L, !anyNA(folds), !anyDuplicated(folds),
            all(folds %in% 1:10))
  workers <- as.integer(Sys.getenv("PIP_KERNEL_WORKERS", "1"))
  stopifnot(!is.na(workers), workers >= 1)
  if (workers > 1L && .Platform$OS.type != "windows") {
    ans <- parallel::mclapply(folds, run_outer,
                               mc.cores = min(workers, length(folds)),
                               mc.set.seed = FALSE)
    if (any(vapply(ans, inherits, logical(1), "try-error")))
      stop("One or more outer folds failed; inspect worker output")
  } else invisible(lapply(folds, run_outer))
}
