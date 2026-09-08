library(ape)
library(phytools)

source("code/fossil_taxonomy.R")
source("code/fossil_placement.R")

base_tree <- read.tree(
  text = "(((Cornus_a:20,Cornus_b:20):60,Crataegus_a:80):20,Other_x:100);"
)
tree_height <- max(nodeHeights(base_tree))
reference_tips <- base_tree$tip.label
name_table <- data.frame(
  genus = c("Cornus", "Crataegus", "Other"),
  family = c("Cornaceae", "Rosaceae", "Nyssaceae"),
  order = c("Cornales", "Rosales", "Cornales"),
  stringsAsFactors = FALSE
)

occurrences <- data.frame(
  fossil_name = c(
    "Cornus_informal_Old_Site",
    "Cornus_repeated_Old_Site",
    "Cornus_repeated_Young_Site",
    "Family_only_Site",
    "Order_only_Site"
  ),
  age_ma = c(60, 60, 10, 30, 30),
  genus = c("unknown", "Cornus", "Cornus", "unknown", "unknown"),
  family = c("unknown", "Cornaceae", "Cornaceae", "Cornaceae", "unknown"),
  order = c("unknown", "Cornales", "Cornales", "Cornales", "Cornales"),
  stringsAsFactors = FALSE
)

place_occurrences <- function(rows) {
  tree <- base_tree
  log <- vector("list", nrow(rows))
  for (index in seq_len(nrow(rows))) {
    result <- graft_fossil_tip(
      tree,
      tip_label = rows$fossil_name[index],
      age_ma = rows$age_ma[index],
      genus = rows$genus[index],
      family = rows$family[index],
      order = rows$order[index],
      name_table = name_table,
      reference_tip_labels = reference_tips,
      tree_height = tree_height,
      placement_fallback = "ancestral_branch"
    )
    tree <- result$tree
    log[[index]] <- data.frame(
      fossil_name = rows$fossil_name[index],
      placed = result$placed,
      level = result$placement_level,
      target = result$placement_target,
      error = result$error,
      stringsAsFactors = FALSE
    )
  }
  list(tree = tree, log = do.call(rbind, log))
}

forward <- place_occurrences(occurrences)
reverse <- place_occurrences(occurrences[nrow(occurrences):1, ])

stopifnot(all(forward$log$placed), all(reverse$log$placed))

# The same named taxon at different sites must produce distinct time-correct
# tips rather than a single mean-age species tip.
for (index in seq_len(nrow(occurrences))) {
  tip <- occurrences$fossil_name[index]
  tip_index <- match(tip, forward$tree$tip.label)
  inferred_age <- tree_height - nodeheight(forward$tree, tip_index)
  stopifnot(isTRUE(all.equal(inferred_age, occurrences$age_ma[index], tolerance = 1e-8)))
}

# A fossil label beginning with Cornus but carrying censored taxonomy must not
# become genus evidence for later occurrences. Placement metadata is invariant
# to input row order and respects the explicit formal ranks.
forward_log <- forward$log[match(occurrences$fossil_name, forward$log$fossil_name), ]
reverse_log <- reverse$log[match(occurrences$fossil_name, reverse$log$fossil_name), ]
row.names(forward_log) <- NULL
row.names(reverse_log) <- NULL
stopifnot(
  identical(forward_log[, c("fossil_name", "level", "target")],
            reverse_log[, c("fossil_name", "level", "target")]),
  forward_log$level[forward_log$fossil_name == "Cornus_informal_Old_Site"] == "root",
  all(forward_log$level[grepl("Cornus_repeated", forward_log$fossil_name)] == "genus"),
  forward_log$level[forward_log$fossil_name == "Family_only_Site"] == "family",
  forward_log$level[forward_log$fossil_name == "Order_only_Site"] == "order"
)

# The extant-to-fossil covariance block—the quantity used by PIP—must not
# depend on fossil input order.
tip_order <- c(reference_tips, occurrences$fossil_name)
vcv_forward <- vcv(keep.tip(forward$tree, tip_order))
vcv_reverse <- vcv(keep.tip(reverse$tree, tip_order))
stopifnot(isTRUE(all.equal(
  vcv_forward[reference_tips, occurrences$fossil_name],
  vcv_reverse[reference_tips, occurrences$fossil_name],
  tolerance = 1e-10
)))

cat("Fossil placement unit checks passed.\n")

