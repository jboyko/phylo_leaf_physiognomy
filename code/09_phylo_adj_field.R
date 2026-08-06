# ==============================================================================
# 09_phylo_adj_field.R
# Visualize the PIP phylogenetic-adjustment field over `tre_pruned`.
#
# For every tip and internal node X of the pruned tree, compute the adjustment
# a hypothetical fossil placed at X would receive:
#
#   phylo_adj(X) = lambda * v_X' %*% V_inv %*% resid
#
# where v_X[i] = depth of MRCA(X, training-tip i). No leave-one-out: a
# hypothetical fossil is genuinely outside the training data.
#
# Output: separate trees for MAT and log(MAP), rendered as fan phylogenies
# with order-level bands around the perimeter.
# ==============================================================================

source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(ape)
library(phytools)

# ==============================================================================
# 1. LOAD
# ==============================================================================

pip <- readRDS("models/pip_components.rds")
phy <- pip$tree_pruned

Ntip  <- length(phy$tip.label)
Nnode <- phy$Nnode

# Depth of every node (root = 0, tips = h for ultrametric tree)
depths   <- node.depth.edgelength(phy)
mrca_mat <- mrca(phy, full = TRUE)  # (Ntip+Nnode) x (Ntip+Nnode) of node IDs

# ==============================================================================
# 2. FIELD COMPUTATION
# ==============================================================================

phylo_adj_field <- function(lambda, V_lam, resid) {
  # V_lam is keyed by training tip labels; align to phy tip order
  tip_order <- phy$tip.label
  V <- V_lam[tip_order, tip_order, drop = FALSE]
  r <- resid[tip_order]

  w <- solve(V, r)  # V_inv %*% resid, length Ntip, aligned to phy$tip.label

  # MRCA depths between every (tip+node) and every training tip
  # mrca_mat[X, 1:Ntip] gives node IDs of MRCAs; look up their depths
  M <- matrix(depths[mrca_mat[, seq_len(Ntip)]],
              nrow = Ntip + Nnode, ncol = Ntip)

  as.numeric(lambda * (M %*% w))
}

adj_mat <- phylo_adj_field(pip$lambda_mat, pip$V_lam_mat, pip$resid_mat)
adj_map <- phylo_adj_field(pip$lambda_map, pip$V_lam_map, pip$resid_map)

# Split into tip values and node values per contMap convention
tips_mat <- setNames(adj_mat[seq_len(Ntip)], phy$tip.label)
node_mat <- setNames(adj_mat[(Ntip + 1):(Ntip + Nnode)],
                     as.character((Ntip + 1):(Ntip + Nnode)))

tips_map <- setNames(adj_map[seq_len(Ntip)], phy$tip.label)
node_map <- setNames(adj_map[(Ntip + 1):(Ntip + Nnode)],
                     as.character((Ntip + 1):(Ntip + Nnode)))

cat("Field computed.\n")
cat("  MAT     adj range: [", round(min(adj_mat), 2), ",",
    round(max(adj_mat), 2), "]\n")
cat("  log MAP adj range: [", round(min(adj_map), 3), ",",
    round(max(adj_map), 3), "]\n")

# ==============================================================================
# 3. ORDER-LEVEL MRCAs FOR LABELING
# ==============================================================================
# Parse genus from each tip label, look up order in name_table_full, take MRCA.

name_tbl <- pip$name_table_full
h_tree   <- max(node.depth.edgelength(phy))

tip_genus_phy <- sapply(strsplit(phy$tip.label, "_"), `[`, 1)
gen2order     <- setNames(name_tbl$order, name_tbl$genus)
tip_order_phy <- gen2order[tip_genus_phy]

order_tbl <- data.frame(tip = phy$tip.label, order = tip_order_phy,
                        stringsAsFactors = FALSE)
order_tbl <- order_tbl[!is.na(order_tbl$order) & nzchar(order_tbl$order), ]

ord_counts <- sort(table(order_tbl$order), decreasing = TRUE)
ord_keep   <- names(ord_counts)[ord_counts >= 5]
cat("Labeling", length(ord_keep), "orders (>=5 tips each)\n")

ord_tip_idx <- lapply(ord_keep, function(o)
  which(phy$tip.label %in% order_tbl$tip[order_tbl$order == o]))
names(ord_tip_idx) <- ord_keep

# ==============================================================================
# 4. RENDER
# ==============================================================================

# Draw a filled annular wedge between radii (r1, r2) and angles (a1, a2),
# centred on (x0, y0), then place a label at the mid-angle just outside.
draw_order_band <- function(x0, y0, r1, r2, a1, a2, fill, label,
                            r_label, cex_label) {
  n <- max(8, ceiling((a2 - a1) / (pi / 180)))  # 1-degree resolution
  outer_a <- seq(a1, a2, length.out = n)
  inner_a <- rev(outer_a)
  xs <- c(x0 + r2 * cos(outer_a), x0 + r1 * cos(inner_a))
  ys <- c(y0 + r2 * sin(outer_a), y0 + r1 * sin(inner_a))
  polygon(xs, ys, col = fill, border = "white", lwd = 0.8)
  am <- (a1 + a2) / 2
  deg <- (am * 180 / pi) %% 360
  srt <- if (deg > 90 && deg < 270) deg - 180 else deg
  text(x0 + r_label * cos(am), y0 + r_label * sin(am),
       labels = label, srt = srt, cex = cex_label, adj = c(0.5, 0.5))
}

render_field <- function(tips_v, node_v, file_base, title) {
  v_abs <- max(abs(c(tips_v, node_v)))
  cm <- contMap(phy, x = tips_v, method = "user",
                anc.states = node_v, plot = FALSE,
                lims = c(-v_abs, v_abs))
  cm <- setMap(cm, c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))

  pdf_path <- paste0("plots/", file_base, ".pdf")
  png_path <- paste0("plots/", file_base, ".png")

  for (dev_call in list(
    function() pdf(pdf_path, width = 12, height = 12),
    function() png(png_path, width = 2400, height = 2400, res = 200)
  )) {
    dev_call()
    par(mar = c(1, 1, 2, 1))
    plot(cm, type = "fan", ftype = "off",
         legend = 0.3 * h_tree, outline = FALSE, lwd = 1.2,
         xlim = c(-1.4 * h_tree, 1.4 * h_tree),
         ylim = c(-1.4 * h_tree, 1.4 * h_tree))
    title(main = title, line = 0, cex.main = 1.0)

    obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    root_idx <- Ntip + 1L
    x0 <- obj$xx[root_idx]; y0 <- obj$yy[root_idx]

    # Tip angle = atan2 from root, then the per-tip "slot" half-width
    tip_x  <- obj$xx[seq_len(Ntip)] - x0
    tip_y  <- obj$yy[seq_len(Ntip)] - y0
    tip_ang <- atan2(tip_y, tip_x)
    # Phytools fan: tips evenly spaced — half-slot width
    half_slot <- pi / Ntip

    r_inner <- 1.03 * h_tree
    r_outer <- 1.13 * h_tree
    r_label <- 1.22 * h_tree

    pal <- c("#d9d9d9", "#bdbdbd")  # alternating shades for visual separation

    # Sort orders by mean tip angle so alternation tracks the perimeter
    ord_meanang <- sapply(ord_tip_idx, function(idx) {
      a <- tip_ang[idx]
      # Use circular mean to handle wrap
      atan2(mean(sin(a)), mean(cos(a)))
    })
    ord_keep_sorted <- ord_keep[order(ord_meanang)]

    for (i in seq_along(ord_keep_sorted)) {
      o   <- ord_keep_sorted[i]
      idx <- ord_tip_idx[[o]]
      a   <- sort(tip_ang[idx])
      a1  <- min(a) - half_slot
      a2  <- max(a) + half_slot
      draw_order_band(x0, y0, r_inner, r_outer, a1, a2,
                      fill = pal[(i %% 2) + 1], label = o,
                      r_label = r_label, cex_label = 0.6)
    }
    dev.off()
  }
  cat("Wrote", pdf_path, "and", png_path, "\n")
}

render_field(tips_mat, node_mat,
             file_base = "fig9_phylo_adj_field_mat",
             title = "Phylogenetic adjustment field - MAT")

render_field(tips_map, node_map,
             file_base = "fig9_phylo_adj_field_map",
             title = "Phylogenetic adjustment field - log(MAP)")

cat("Done.\n")
