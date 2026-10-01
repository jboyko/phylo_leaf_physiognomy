# Extant held-site transfer check with genus/family/order-only placement.
# All placement anchors belong to the outer training species.
Sys.setenv(PIP_KERNEL_SETUP_ONLY = "1")
source("code/03h_map_kernel_cv.R")
source("code/fossil_placement.R")
paths <- file.path(out_dir, sprintf("outer_%02d.rds", 1:10))
fold_arg <- Sys.getenv("PIP_KERNEL_TRANSFER_FOLDS", "")
target_folds <- if (nzchar(fold_arg))
  as.integer(strsplit(fold_arg, ",", fixed = TRUE)[[1]]) else 1:10
stopifnot(length(target_folds) > 0L, !anyNA(target_folds),
          !anyDuplicated(target_folds), all(target_folds %in% 1:10))
if (!all(file.exists(paths[target_folds])))
  stop("Complete the requested extant kernel folds first")
name_table <- read.csv("data/name_table_full.csv", stringsAsFactors = FALSE)
records <- list()
placement_records <- list()
scaffold_tips <- scaffold$tip.label
for (outer in target_folds) {
  checkpoint <- readRDS(paths[outer])
  d <- prepare_kernel_split(checkpoint$train_sites, checkpoint$held_sites,
                            42 + outer)
  held <- raw_dat[raw_dat$Site %in% checkpoint$held_sites,
                  c("Site", "genusSpecies", "genus", "Family", "Order")]
  key_raw <- paste(held$Site, held$genusSpecies, sep = "\r")
  collapse_rank <- function(v) {
    known <- sort(unique(v[!is.na(v) & nzchar(v) & v != "unknown"]))
    # Operational genus-only labels can pool unrelated morphotypes. A rank
    # with contradictory evidence is unavailable for this occurrence.
    if (length(known) == 1L) known else "unknown"
  }
  taxonomy <- do.call(rbind, lapply(split(held, key_raw), function(a)
    data.frame(Site = a$Site[1], genusSpecies = a$genusSpecies[1],
               genus = collapse_rank(a$genus),
               Family = collapse_rank(a$Family),
               Order = collapse_rank(a$Order))))
  key <- paste(taxonomy$Site, taxonomy$genusSpecies, sep = "\r")
  lookup <- match(paste(d$occ$site, d$occ$species, sep = "\r"), key)
  stopifnot(!anyNA(lookup))
  taxonomy <- taxonomy[lookup, ]
  for (level in c("genus", "family", "order")) {
    names_new <- sprintf("cv_%02d_%s_%04d", outer, level,
                         seq_len(nrow(d$occ)))
    placements <- character(nrow(d$occ))
    b_cross <- matrix(NA_real_, nrow(d$occ), length(d$ids))
    train_genera <- sub("_.*$", "", d$ids)
    verified_geometry <- character()
    # Place each query against the unchanged scaffold. Extant tips differ in
    # depth by a few machine-scale branch units; 1e-5 Ma (10 years) avoids an
    # edge-rounding failure at age exactly zero and is negligible here.
    for (j in seq_len(nrow(d$occ))) {
      tree <- scaffold
      genus <- if (level == "genus") taxonomy$genus[j] else "unknown"
      family <- if (level != "order") taxonomy$Family[j] else "unknown"
      p <- graft_fossil_tip(tree, d$ids, names_new[j], 1e-5,
                            genus, family, taxonomy$Order[j],
                            placement_fallback = "ancestral_branch",
                            name_table = name_table)
      if (!p$placed) stop("CV placement failed: outer=", outer,
                           " level=", level, " occurrence=", j,
                           " site=", d$occ$site[j],
                           " species=", d$occ$species[j],
                           " genus=", genus, " family=", family,
                           " order=", taxonomy$Order[j], ": ", p$error)
      placements[j] <- p$placement_level
      # The new tip branches from one scaffold lineage. Its shared path with
      # each training tip is the shorter of that lineage's shared path and
      # the insertion depth. This avoids materializing a full VCV per query.
      anchor_genus <- switch(p$placement_level,
        genus = p$placement_target,
        family = name_table$genus[name_table$family == p$placement_target],
        order = name_table$genus[name_table$order == p$placement_target],
        root = train_genera)
      anchor <- d$ids[train_genera %in% anchor_genus][1]
      if (is.na(anchor)) stop("No training anchor after placement")
      new_tip <- match(names_new[j], p$tree$tip.label)
      parent <- p$tree$edge[p$tree$edge[, 2] == new_tip, 1]
      insertion_depth <- phytools::nodeheight(p$tree, parent) / root_height
      b_cross[j, ] <- pmin(b_scaffold[anchor, d$ids], insertion_depth)
      geometry <- paste(p$placement_level, p$age_fallback)
      if (!geometry %in% verified_geometry) {
        vcv <- ape::vcv(p$tree)
        stopifnot(max(abs(vcv[d$ids, d$ids, drop = FALSE] /
                              root_height - d$b)) < 1e-8,
                  max(abs(b_cross[j, ] -
                            vcv[names_new[j], d$ids] / root_height)) < 1e-8)
        verified_geometry <- c(verified_geometry, geometry)
      }
    }
    placement_records[[length(placement_records) + 1L]] <- data.frame(
      fold = outer, requested_resolution = level,
      site = d$occ$site, species = d$occ$species,
      actual_placement_level = placements)
    stopifnot(all(is.finite(b_cross)))
    for (k in seq_len(nrow(checkpoint$selected))) {
      cfg <- checkpoint$selected[k, ]
      ell <- if (is.na(cfg$ell_factor)) NA_real_ else
        cfg$ell_factor * d$reference
      fit <- map_kernel_fit(d$y, d$z, d$b, cfg$lambda_phy, cfg$eta,
                            cfg$kind, cfg$mean_kind, ell)
      p <- map_kernel_predict(fit, d$new_z, b_cross)
      site_pred <- tapply(p$prediction, d$occ$site, mean)
      records[[length(records) + 1L]] <- data.frame(
        fold = outer, candidate = cfg$candidate, placement_resolution = level,
        site = names(site_pred), prediction = as.numeric(site_pred),
        observed = dat_site_obs[names(site_pred), "log_map"],
        map_cm = exp(as.numeric(site_pred)))
    }
    cat("outer", outer, "coarsened", level, "complete\n")
  }
}
out <- do.call(rbind, records)
expected_sites <- sum(vapply(target_folds, function(i)
  length(readRDS(paths[i])$held_sites), integer(1)))
stopifnot(nrow(out) == expected_sites * 4L * 3L,
          all(is.finite(out$prediction)),
          !anyDuplicated(out[, c("candidate", "placement_resolution", "site")]))
dir.create("tables/map_kernel_experiment", recursive = TRUE,
           showWarnings = FALSE)
write.csv(out, "tables/map_kernel_experiment/coarsened_site_predictions.csv",
          row.names = FALSE)
write.csv(do.call(rbind, placement_records),
          "tables/map_kernel_experiment/coarsened_placement_log.csv",
          row.names = FALSE)
metrics <- do.call(rbind, lapply(split(out,
  interaction(out$candidate, out$placement_resolution, drop = TRUE)),
  function(d) data.frame(candidate = d$candidate[1],
    placement_resolution = d$placement_resolution[1],
    rmse = sqrt(mean((d$prediction - d$observed)^2)),
    mae = mean(abs(d$prediction - d$observed)),
    mean_signed_error = mean(d$prediction - d$observed))))
write.csv(metrics, "tables/map_kernel_experiment/coarsened_comparison.csv",
          row.names = FALSE)
print(metrics, row.names = FALSE)
