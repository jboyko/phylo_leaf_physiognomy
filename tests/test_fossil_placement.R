source("code/fossil_placement.R")
source("code/fossil_taxonomy.R")

library(ape)
library(phytools)

fixture_tree <- read.tree(text = "((A_one:5,A_two:5):5,B_one:10);")
fixture_scaffold <- fixture_tree$tip.label
fixture_taxonomy <- data.frame(
  genus = c("A", "A", "B"),
  family = c("AFam", "AFam", "BFam"),
  order = c("AOrd", "AOrd", "BOrd"),
  stringsAsFactors = FALSE
)

place_fixture_rows <- function(rows) {
  tree <- fixture_tree
  for (i in seq_len(nrow(rows))) {
    res <- graft_fossil_tip(
      tree, fixture_scaffold, rows$tip_label[i], rows$age_ma[i],
      rows$genus[i], rows$family[i], rows$order[i],
      name_table = fixture_taxonomy
    )
    stopifnot(isTRUE(res$placed))
    tree <- res$tree
  }
  tree
}

# A fossil older than the A crown is grafted on the ancestral edge at the
# depth defined by its own age.
age_rows <- data.frame(
  tip_label = c("A_old", "A_young"),
  age_ma = c(7, 2), genus = "A", family = "AFam", order = "AOrd",
  stringsAsFactors = FALSE
)
age_tree <- place_fixture_rows(age_rows)
expected_heights <- 10 - age_rows$age_ma
actual_heights <- vapply(
  age_rows$tip_label,
  function(label) nodeheight(age_tree, match(label, age_tree$tip.label)),
  numeric(1)
)
stopifnot(isTRUE(all.equal(unname(actual_heights), expected_heights, tolerance = 1e-10)))

# The fixed scaffold prevents an earlier A fossil from becoming an A anchor.
# Reversing input order preserves every fossil-to-extant covariance.
reverse_tree <- place_fixture_rows(age_rows[2:1, ])
extant <- fixture_scaffold
forward_vcv <- vcv(age_tree)[age_rows$tip_label, extant, drop = FALSE]
reverse_vcv <- vcv(reverse_tree)[age_rows$tip_label, extant, drop = FALSE]
stopifnot(isTRUE(all.equal(forward_vcv, reverse_vcv, tolerance = 1e-10)))

# The two previously dropped primary occurrences now both graft successfully
# and retain their individual site ages.
pipeline_tree <- read.tree("data/tre_scaffold.tre")
pipeline_scaffold <- pipeline_tree$tip.label
pipeline_taxonomy <- readRDS("models/pip_components.rds")$name_table_full
fossils <- read.csv("data/fossil_traits.csv", stringsAsFactors = FALSE)
problem_rows <- fossils[fossils$species %in%
  c("Cornus nebrascensis", "Crataegus sp. (rp42)"), ]
stopifnot(nrow(problem_rows) == 2L)
for (i in seq_len(nrow(problem_rows))) {
  row <- problem_rows[i, ]
  label <- paste0("placement_regression_", i)
  res <- graft_fossil_tip(
    pipeline_tree, pipeline_scaffold, label, row$age_ma,
    row$genus, row$family, row$order, name_table = pipeline_taxonomy
  )
  stopifnot(isTRUE(res$placed))
  expected_height <- max(nodeHeights(pipeline_tree)) - row$age_ma
  actual_height <- nodeheight(res$tree, match(label, res$tree$tip.label))
  stopifnot(isTRUE(all.equal(actual_height, expected_height, tolerance = 1e-8)))
  pipeline_tree <- res$tree
}

# Replay every formal-only occurrence in original, reversed, and shuffled
# order.  Each final tip must retain its own age and the fossil-to-extant VCV
# block must be invariant to insertion order.
formal_rows <- taxonomy_for_scenario(fossils, "formal_only")
order_rows <- formal_rows
stopifnot(nrow(order_rows) == 361L)
fossil_labels <- order_rows$fossil_name
extant_labels <- rownames(readRDS("models/pip_components.rds")$dat_imputed_mat)
place_all_occurrences <- function(rows) {
  tree <- read.tree("data/tre_scaffold.tre")
  scaffold <- tree$tip.label
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    result <- graft_fossil_tip(
      tree, scaffold, row$fossil_name, row$age_ma,
      row$genus, row$family, row$order, name_table = pipeline_taxonomy
    )
    stopifnot(isTRUE(result$placed))
    tree <- result$tree
    expected <- max(nodeHeights(tree)) - row$age_ma
    observed <- nodeheight(tree, match(row$fossil_name, tree$tip.label))
    stopifnot(isTRUE(all.equal(observed, expected, tolerance = 1e-8)))
  }
  final_vcv <- vcv(tree)
  expected_heights <- max(nodeHeights(tree)) - order_rows$age_ma
  stopifnot(isTRUE(all.equal(unname(diag(final_vcv)[fossil_labels]),
                            expected_heights, tolerance = 1e-8)))
  final_vcv[fossil_labels, extant_labels, drop = FALSE]
}

original_block <- place_all_occurrences(order_rows)
reversed_block <- place_all_occurrences(order_rows[nrow(order_rows):1, ])
set.seed(20260907)
shuffled_block <- place_all_occurrences(order_rows[sample(nrow(order_rows)), ])
stopifnot(
  isTRUE(all.equal(original_block, reversed_block, tolerance = 1e-10)),
  isTRUE(all.equal(original_block, shuffled_block, tolerance = 1e-10))
)

# An age older than the root has no branch at which to attach.  Preserve the
# input tree rather than silently inventing a different age.
too_old <- graft_fossil_tip(
  fixture_tree, fixture_scaffold, "too_old", 11, "A", "AFam", "AOrd",
  name_table = fixture_taxonomy
)
stopifnot(!too_old$placed, identical(too_old$tree, fixture_tree))

cat("Fossil placement regression checks passed.\n")
