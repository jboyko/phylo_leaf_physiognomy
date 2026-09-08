# Shared helpers for time-aware fossil-tip placement.
#
# Placement targets are always resolved from the original extant scaffold tips.
# Fossils grafted earlier in a loop must never become evidence for the placement
# of a later fossil: doing so makes results depend on input row order and can
# leak informal fossil labels back into a formal-only taxonomy scenario.

fossil_tip_genera <- function(tip_labels) {
  vapply(strsplit(tip_labels, "_", fixed = TRUE), `[`, character(1), 1L)
}

fossil_target_node <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  if (length(tips) >= 2L) return(ape::getMRCA(tree, tips))
  if (length(tips) == 1L) return(match(tips, tree$tip.label))
  NULL
}

resolve_fossil_target <- function(tree, genus, family, order, name_table,
                                  reference_tip_labels) {
  reference_tips <- intersect(reference_tip_labels, tree$tip.label)
  reference_genera <- fossil_tip_genera(reference_tips)

  if (is_known_taxon(genus)) {
    genus_tips <- reference_tips[reference_genera == genus]
    target <- fossil_target_node(tree, genus_tips)
    if (!is.null(target)) {
      return(list(node = target, level = "genus", target = genus))
    }
  }

  if (is_known_taxon(family)) {
    family_genera <- unique(name_table$genus[name_table$family == family])
    family_tips <- reference_tips[reference_genera %in% family_genera]
    target <- fossil_target_node(tree, family_tips)
    if (!is.null(target)) {
      return(list(node = target, level = "family", target = family))
    }
  }

  if (is_known_taxon(order)) {
    order_genera <- unique(name_table$genus[name_table$order == order])
    order_tips <- reference_tips[reference_genera %in% order_genera]
    target <- fossil_target_node(tree, order_tips)
    if (!is.null(target)) {
      return(list(node = target, level = "order", target = order))
    }
  }

  list(
    node = ape::Ntip(tree) + 1L,
    level = "root",
    target = "root"
  )
}

graft_fossil_tip <- function(tree, tip_label, age_ma, genus, family, order,
                             name_table, reference_tip_labels,
                             tree_height = max(phytools::nodeHeights(tree)),
                             placement_fallback = "ancestral_branch") {
  if (tip_label %in% tree$tip.label) {
    return(list(
      tree = tree,
      placed = TRUE,
      placement_level = "existing_tip",
      placement_target = tip_label,
      age_fallback = "none",
      error = NA_character_
    ))
  }

  tryCatch({
    if (!is.finite(age_ma)) stop("age_ma must be finite")

    resolved <- resolve_fossil_target(
      tree, genus, family, order, name_table, reference_tip_labels
    )
    target_node <- resolved$node
    placement_level <- resolved$level
    placement_target <- resolved$target
    age_fallback <- "none"

    if (placement_level == "root") {
      warning("No taxonomy match for '", tip_label, "'. Placing at root.")
    }

    fossil_height <- tree_height - age_ma
    node_height <- phytools::nodeheight(tree, target_node)
    terminal_length <- fossil_height - node_height
    tolerance <- sqrt(.Machine$double.eps) * max(1, tree_height)

    if (terminal_length >= -tolerance) {
      tree <- phytools::bind.tip(
        tree,
        tip.label = tip_label,
        where = target_node,
        edge.length = max(terminal_length, 0)
      )
    } else if (placement_fallback == "ancestral_branch") {
      age_fallback <- "ancestral_branch"
      node <- target_node
      found <- FALSE

      repeat {
        edge_idx <- which(tree$edge[, 2] == node)
        if (length(edge_idx) == 0L) break

        parent_node <- tree$edge[edge_idx, 1]
        parent_height <- phytools::nodeheight(tree, parent_node)
        child_height <- phytools::nodeheight(tree, node)
        branch_length <- tree$edge.length[edge_idx]

        spans_fossil <-
          parent_height <= fossil_height + tolerance &&
          child_height >= fossil_height - tolerance

        if (spans_fossil) {
          position <- child_height - fossil_height
          if (position < -tolerance || position > branch_length + tolerance) {
            stop("computed fossil position falls outside the spanning branch")
          }
          position <- min(max(position, 0), branch_length)
          tree <- phytools::bind.tip(
            tree,
            tip.label = tip_label,
            where = node,
            position = position,
            edge.length = 0
          )
          found <- TRUE
          break
        }
        node <- parent_node
      }

      if (!found) {
        warning("No spanning branch for '", tip_label, "'. Placing at root.")
        placement_level <- "root"
        placement_target <- "root"
        age_fallback <- "root"
        tree <- phytools::bind.tip(
          tree,
          tip.label = tip_label,
          where = ape::Ntip(tree) + 1L,
          edge.length = 0.001
        )
      }
    } else if (placement_fallback == "node") {
      age_fallback <- "node"
      tree <- phytools::bind.tip(
        tree,
        tip.label = tip_label,
        where = target_node,
        edge.length = 0.001
      )
    } else {
      stop("Unknown placement_fallback: ", placement_fallback)
    }

    list(
      tree = tree,
      placed = TRUE,
      placement_level = placement_level,
      placement_target = placement_target,
      age_fallback = age_fallback,
      error = NA_character_
    )
  }, error = function(error) {
    list(
      tree = tree,
      placed = FALSE,
      placement_level = NA_character_,
      placement_target = NA_character_,
      age_fallback = NA_character_,
      error = conditionMessage(error)
    )
  })
}
