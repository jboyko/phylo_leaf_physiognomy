# Shared fossil placement helpers.
#
# `scaffold_tip_labels` is deliberately an immutable snapshot of the extant
# scaffold.  Fossils already grafted into `tree` must never become taxonomic
# evidence for later fossil placements.

fossil_placement_failure <- function(tree, message) {
  list(
    tree = tree,
    placed = FALSE,
    placement_level = NA_character_,
    placement_target = NA_character_,
    age_fallback = NA_character_,
    error = message
  )
}

graft_fossil_tip <- function(tree, scaffold_tip_labels, tip_label, age_ma,
                             genus, family, order,
                             placement_fallback = "ancestral_branch",
                             name_table) {
  if (!is.character(scaffold_tip_labels) || !length(scaffold_tip_labels)) {
    return(fossil_placement_failure(tree, "scaffold_tip_labels must be non-empty"))
  }
  if (!all(c("genus", "family", "order") %in% names(name_table))) {
    return(fossil_placement_failure(
      tree, "name_table must contain genus, family, and order columns"
    ))
  }
  if (length(tip_label) != 1L || is.na(tip_label) || !nzchar(tip_label)) {
    return(fossil_placement_failure(tree, "tip_label must be one non-empty string"))
  }
  if (length(age_ma) != 1L || !is.finite(age_ma) || age_ma < 0) {
    return(fossil_placement_failure(tree, "age_ma must be one finite, non-negative value"))
  }
  if (!placement_fallback %in% c("ancestral_branch", "node")) {
    return(fossil_placement_failure(
      tree, "placement_fallback must be 'ancestral_branch' or 'node'"
    ))
  }
  if (tip_label %in% tree$tip.label) {
    return(list(
      tree = tree, placed = TRUE, placement_level = "existing_tip",
      placement_target = tip_label, age_fallback = "none", error = NA_character_
    ))
  }

  # Restrict every taxonomic lookup to the original scaffold.  `tree` grows
  # after each call, but this candidate set does not.
  anchor_tips <- intersect(scaffold_tip_labels, tree$tip.label)
  if (!length(anchor_tips)) {
    return(fossil_placement_failure(tree, "no scaffold tips remain in tree"))
  }
  tip_genera <- sub("_.*$", "", anchor_tips)

  known <- function(x) {
    length(x) == 1L && !is.na(x) && nzchar(x) && x != "unknown"
  }
  target_node <- NULL
  placement_level <- "root"
  placement_target <- "root"
  age_fallback <- "none"

  lookup_target <- function(candidate_tips) {
    candidate_tips <- intersect(anchor_tips, candidate_tips)
    if (length(candidate_tips) >= 2L) return(ape::getMRCA(tree, candidate_tips))
    if (length(candidate_tips) == 1L) return(match(candidate_tips, tree$tip.label))
    NULL
  }

  if (known(genus)) {
    target_node <- lookup_target(anchor_tips[tip_genera == genus])
    if (!is.null(target_node)) {
      placement_level <- "genus"
      placement_target <- genus
    }
  }
  if (is.null(target_node) && known(family)) {
    family_genera <- unique(name_table$genus[name_table$family == family])
    target_node <- lookup_target(anchor_tips[tip_genera %in% family_genera])
    if (!is.null(target_node)) {
      placement_level <- "family"
      placement_target <- family
    }
  }
  if (is.null(target_node) && known(order)) {
    order_genera <- unique(name_table$genus[name_table$order == order])
    target_node <- lookup_target(anchor_tips[tip_genera %in% order_genera])
    if (!is.null(target_node)) {
      placement_level <- "order"
      placement_target <- order
    }
  }
  if (is.null(target_node)) target_node <- length(tree$tip.label) + 1L

  tree_before <- tree
  tryCatch({
    tree_height <- max(phytools::nodeHeights(tree))
    fossil_height <- tree_height - age_ma
    node_height <- phytools::nodeheight(tree, target_node)
    edge_length <- fossil_height - node_height

    if (edge_length >= 0) {
      tree <- phytools::bind.tip(
        tree, tip.label = tip_label, where = target_node, edge.length = edge_length
      )
    } else if (placement_fallback == "ancestral_branch") {
      age_fallback <- "ancestral_branch"
      node <- target_node
      found <- FALSE
      repeat {
        edge_idx <- which(tree$edge[, 2] == node)
        if (!length(edge_idx)) break
        parent_node <- tree$edge[edge_idx, 1]
        parent_height <- phytools::nodeheight(tree, parent_node)
        if (parent_height <= fossil_height) {
          # bind.tip expects a position along this exact parent edge.  Permit
          # only machine-scale roundoff at an endpoint; a material mismatch
          # signals a broken time geometry and must not be silently altered.
          edge_limit <- tree$edge.length[edge_idx[1]]
          position <- phytools::nodeheight(tree, node) - fossil_height
          tolerance <- sqrt(.Machine$double.eps) * max(1, abs(edge_limit))
          if (position < -tolerance || position > edge_limit + tolerance) {
            stop("computed fossil position lies outside its parent edge")
          }
          position <- max(0, min(position, edge_limit))
          tree <- phytools::bind.tip(
            tree, tip.label = tip_label, where = node, position = position,
            edge.length = 0
          )
          found <- TRUE
          break
        }
        node <- parent_node
      }
      if (!found) {
        stop(
          "fossil age predates the scaffold root; no branch exists at its age"
        )
      }
    } else {
      age_fallback <- "node"
      tree <- phytools::bind.tip(
        tree, tip.label = tip_label, where = target_node, edge.length = 0.001
      )
    }
    placed_tip <- match(tip_label, tree$tip.label)
    placed_height <- phytools::nodeheight(tree, placed_tip)
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(fossil_height))
    if (age_fallback %in% c("none", "ancestral_branch") &&
        abs(placed_height - fossil_height) > tolerance) {
      stop("grafted tip depth does not match scaffold height minus fossil age")
    }
    list(
      tree = tree, placed = TRUE, placement_level = placement_level,
      placement_target = placement_target, age_fallback = age_fallback,
      error = NA_character_
    )
  }, error = function(e) fossil_placement_failure(tree_before, conditionMessage(e)))
}
