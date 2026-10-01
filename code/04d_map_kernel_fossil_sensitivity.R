# Exploratory fossil point estimates for one explicitly selected kernel family.
# Requires completed nested extant experiment. Does not replace production MAP.
candidate <- Sys.getenv("PIP_KERNEL_FOSSIL_CANDIDATE", "")
if (!candidate %in% c("linear_kernel", "rbf_kernel", "linear_mean_rbf",
                      "linear_mean_phylogeny"))
  stop("Set PIP_KERNEL_FOSSIL_CANDIDATE to one documented kernel family")
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
stopifnot(found_stop, length(all_sites) == 92L)
source("code/map_kernel_model.R")
source("code/fossil_taxonomy.R")
source("code/fossil_placement.R")
paths <- file.path("models/map_kernel_experiment",
                   sprintf("outer_%02d.rds", 1:10))
if (!all(file.exists(paths))) stop("Complete the nested experiment first")
folds <- lapply(paths, readRDS)
config <- folds[[1]]$grid
stopifnot(all(vapply(folds, function(f) identical(f$grid, config), logical(1))))
# Every extant site is an inner validation site nine times across the ten
# outer-training subsets. Pool squared errors, not fold-level RMSEs.
losses <- Reduce(`+`, lapply(folds, function(f) rowSums(f$inner_losses)))
totals <- sum(vapply(folds, function(f) length(f$train_sites), integer(1)))
config$pooled_inner_mse <- losses / totals
subset <- config[config$candidate == candidate, ]
chosen <- subset[which.min(subset$pooled_inner_mse), ]

filled <- fill_tooth_traits(raw_dat)
sp <- agg_species(filled)
scaffold <- ape::read.tree("data/tre_scaffold.tre")
scaffold_tips <- scaffold$tip.label
sp <- sp[intersect(rownames(sp), scaffold_tips), , drop = FALSE]
pn <- active_pred_names(sp, fossil_traits)
stopifnot(length(pn) > 0L)
set.seed(42)
imputer <- caret::preProcess(sp[, pn, drop = FALSE], method = "bagImpute")
x <- predict(imputer, sp[, pn, drop = FALSE])
assert_finite(x, "final extant traits")
scaling <- map_kernel_scale_fit(x)
z <- map_kernel_scale_apply(x, scaling)
colnames(z) <- scaling$retained
root_height <- max(phytools::nodeHeights(scaffold))
b <- ape::vcv(scaffold)[rownames(sp), rownames(sp), drop = FALSE] /
  root_height
reference <- map_kernel_reference(z)
ell <- if (is.na(chosen$ell_factor)) NA_real_ else
  chosen$ell_factor * reference
fit <- map_kernel_fit(sp$log_map, z, b, chosen$lambda_phy,
                      chosen$eta, chosen$kind, chosen$mean_kind, ell)
fit$training_ids <- rownames(sp)
fit$predictors <- pn
fit$scaling <- scaling
fit$root_height <- root_height
fit$selection <- chosen
fit$imputer <- imputer
saveRDS(fit, "models/map_kernel_experiment/final_fit.rds")

foss_base <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)
expected_ages <- fossil_site_ages()[foss_base$site]
stopifnot(nrow(foss_base) == 360L, length(unique(foss_base$site)) == 10L,
          !anyNA(expected_ages),
          max(abs(foss_base$age_ma - expected_ages)) < 1e-8,
          all(foss_base$site != "Palacio de los Loros PL1"),
          all(foss_base$site != "Palacio de los Loros PL2"))
pip <- readRDS("models/pip_components.rds")
name_table <- pip$name_table_full
old <- read.csv("tables/fossil_map_recalibration_trial.csv",
                stringsAsFactors = FALSE)
all_occ <- list(); all_site <- list(); all_placement <- list()
for (scenario in c("formal_only", "include_informal")) {
  foss <- taxonomy_for_scenario(foss_base, scenario)
  tree <- scaffold
  placements <- vector("list", nrow(foss))
  for (j in seq_len(nrow(foss))) {
    p <- graft_fossil_tip(tree, scaffold_tips, foss$fossil_name[j],
                          foss$age_ma[j], foss$genus[j], foss$family[j],
                          foss$order[j], placement_fallback =
                            "ancestral_branch", name_table = name_table)
    if (!p$placed) stop("Fossil placement failed for ", foss$fossil_name[j],
                         ": ", p$error)
    tree <- p$tree
    placements[[j]] <- data.frame(taxonomy_scenario = scenario,
      fossil_name = foss$fossil_name[j], site = foss$site[j],
      placement_level = p$placement_level,
      placement_target = p$placement_target,
      age_fallback = p$age_fallback)
  }
  vcv <- ape::vcv(tree)
  check_b <- vcv[fit$training_ids, fit$training_ids, drop = FALSE] /
    root_height
  if (max(abs(check_b - b)) > 1e-8)
    stop("Grafting changed calibration covariance/root geometry")
  raw_traits <- as.data.frame(matrix(NA_real_, nrow(foss), length(pn),
                                     dimnames = list(NULL, pn)))
  present <- intersect(pn, names(foss))
  raw_traits[, present] <- foss[, present, drop = FALSE]
  completed <- predict(imputer, raw_traits)
  assert_finite(completed, "fossil traits")
  new_z <- map_kernel_scale_apply(completed, scaling)
  colnames(new_z) <- scaling$retained
  cross <- vcv[foss$fossil_name, fit$training_ids, drop = FALSE] /
    root_height
  solo <- graft_fossil_tip(scaffold, scaffold_tips,
    foss$fossil_name[1], foss$age_ma[1], foss$genus[1],
    foss$family[1], foss$order[1], placement_fallback =
      "ancestral_branch", name_table = name_table)
  stopifnot(solo$placed)
  solo_cross <- ape::vcv(solo$tree)[foss$fossil_name[1],
                                      fit$training_ids] / root_height
  stopifnot(max(abs(cross[1, ] - solo_cross)) < 1e-8)
  components <- map_kernel_predict(fit, new_z, cross)
  distance <- map_kernel_dist2(new_z, z)
  nearest <- sqrt(apply(distance, 1, min))
  similarity <- if (chosen$kind == "rbf")
    exp(-nearest^2 / (2 * fit$ell^2)) else rep(NA_real_, nrow(foss))
  occ <- cbind(foss[, c("fossil_name", "species", "site", "age_ma")],
               taxonomy_scenario = scenario, candidate = candidate,
               components, map_cm = exp(components$prediction),
               nearest_scaled_trait_distance = nearest,
               max_kernel_similarity = similarity,
               weak_kernel_similarity = similarity < 0.5,
               outside_training_trait_range = apply(new_z, 1, function(v)
                 any(v < apply(z, 2, min) | v > apply(z, 2, max))))
  site <- aggregate(occ[, c("linear_mean", "trait_adjustment",
                            "phylo_adjustment", "prediction")],
                    list(taxonomy_scenario = rep(scenario, nrow(occ)),
                         site = occ$site,
                         age_ma = occ$age_ma), mean)
  site$map_cm <- exp(site$prediction)
  comparator <- old[old$taxonomy_scenario == scenario,
                    c("site", "original_map_cm", "recalibrated_map_cm")]
  site <- merge(site, comparator, by = "site", sort = FALSE)
  stopifnot(nrow(site) == 10L,
            all(is.finite(site$prediction)),
            max(abs(site$map_cm - exp(site$prediction))) < 1e-10)
  all_occ[[scenario]] <- occ
  all_site[[scenario]] <- site
  all_placement[[scenario]] <- do.call(rbind, placements)
}
out <- "tables/map_kernel_experiment"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, all_occ), file.path(out,
  "fossil_occurrence_sensitivity.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_site), file.path(out,
  "fossil_site_sensitivity.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_placement), file.path(out,
  "fossil_placement_sensitivity.csv"), row.names = FALSE)
write.csv(chosen, file.path(out, "fossil_selected_parameters.csv"),
          row.names = FALSE)
cat("Exploratory fossil point estimates saved for", candidate, "\n")
