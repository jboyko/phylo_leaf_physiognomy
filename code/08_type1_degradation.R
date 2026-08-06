source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")

library(ape)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# 1. LOAD PIP COMPONENTS
# ==============================================================================

pip <- readRDS("models/pip_components.rds")

phylomat <- vcv(pip$tree_pruned)
diag(phylomat) <- diag(phylomat) + 1e-6

tree    <- pip$tree_pruned
tax     <- pip$taxonomy          # species, genus, family, order
species <- rownames(pip$V_lam_mat)
n       <- length(species)

cat("n species:", n, "\n")
cat("n unique genera:", length(unique(tax$genus)), "\n")
cat("n unique families:", length(unique(tax$family)), "\n")
cat("n unique orders:", length(unique(tax$order)), "\n")

# Branching times (internal node ages, Ma from root)
bt <- branching.times(tree)

# Per-target components
targets <- list(
  mat = list(
    V_lam  = pip$V_lam_mat,
    resid  = pip$resid_mat[species],
    lambda = pip$lambda_mat,
    label  = "MAT (°C)"
  ),
  map = list(
    V_lam  = pip$V_lam_map,
    resid  = pip$resid_map[species],
    lambda = pip$lambda_map,
    label  = "log(MAP)"
  )
)

# ==============================================================================
# 2. HELPER — maps lambda-scaled cross-covariance to PIP adjustment
#    (reused from 07_type2_degradation.R)
# ==============================================================================

pip_adj_from_vcov <- function(v_cross, V_inv, resid) {
  as.numeric(t(v_cross) %*% V_inv %*% resid)
}

# ==============================================================================
# 3. BUILD CLADE TABLES — for each rank, find the MRCA node and its depth
#    for every multi-species clade present in the training set.
#
#    Placement covariance for a fossil at MRCA node M (depth d_M):
#      - For training tip T descending from M: cov(fossil, T) = d_M
#      - For training tip T outside M's clade: cov(fossil, T) = phylomat[rep, T]
#        where rep is any member of the focal clade (same value for all members
#        because MRCA(clade, T) is the same regardless of which member we pick).
# ==============================================================================

# Lookup: species index in the species vector
sp_idx <- setNames(seq_len(n), species)

# Map from species name to its row in taxonomy
tax_sp <- setNames(seq_len(nrow(tax)), tax$species)

build_clade_table <- function(rank_col) {
  clade_members <- split(tax$species, tax[[rank_col]])
  rows <- list()
  skipped_singleton <- 0L

  for (clade_name in names(clade_members)) {
    # Training species in this clade
    in_training <- intersect(clade_members[[clade_name]], species)
    if (length(in_training) < 2) {
      skipped_singleton <- skipped_singleton + length(in_training)
      next
    }

    # MRCA node and its depth
    mrca_node  <- getMRCA(tree, in_training)
    depth_M    <- bt[as.character(mrca_node)]

    # Descendants of MRCA that are in the training set
    desc_all   <- extract.clade(tree, mrca_node)$tip.label
    desc_train <- intersect(desc_all, species)

    rows[[clade_name]] <- list(
      clade      = clade_name,
      mrca_node  = mrca_node,
      depth_M    = depth_M,
      in_training = in_training,   # focal species for this clade
      desc_train  = desc_train,    # all descendants in training set
      rep_sp      = in_training[1] # representative for non-descendant covs
    )
  }

  cat("  Rank:", rank_col, "| clades with 2+ training species:",
      length(rows), "| singleton focal species skipped:", skipped_singleton, "\n")
  rows
}

cat("Building clade tables...\n")
clade_tables <- list(
  genus  = build_clade_table("genus"),
  family = build_clade_table("family"),
  order  = build_clade_table("order")
)

# ==============================================================================
# 4. COMPUTE LOO ADJUSTMENTS FOR EACH TARGET × RANK × FOCAL
#
#  For a focal species i placed at its rank-level MRCA (node M):
#
#    adj_loo_i(M) = adj_full(M) - w[i] * c_M[i]
#
#  where:
#    v_M      = lambda * cov_M          (lambda-scaled cross-cov vector)
#    adj_full = t(v_M) %*% p            (full, non-LOO adjustment)
#    c_M      = V_inv %*% v_M           (n-vector; need only the focal's element)
#    w[i]     = p[i] / V_inv[i,i]       (leverage weight)
#    p        = V_inv %*% resid
#
#  Error (b): marginal — adj_loo_i(MRCA) minus adj_correct_i
#  Error (a): absolute — adj_loo_i(MRCA) minus resid[i]
#
#  adj_correct_i = adj_full(tip_i) - w[i] * c_tip_i[i]
#               = a[i]            - lambda * w[i] * E[i,i]
#  (same formula as Type 2 diagonal, recomputed here for self-containment)
# ==============================================================================

cat("Computing LOO adjustments...\n")

all_results <- list()

for (tname in names(targets)) {
  cat(" Processing target:", tname, "\n")
  comp   <- targets[[tname]]
  V_lam  <- comp$V_lam[species, species]
  resid  <- comp$resid
  lambda <- comp$lambda

  V_inv  <- solve(V_lam)
  p      <- as.numeric(V_inv %*% resid)
  names(p) <- species

  # Leverage weights
  diag_Vinv <- diag(V_inv)
  w <- p / diag_Vinv

  # Full adjustment vector (placement at own tip, no LOO) — reused from type2
  E           <- phylomat[species, species] %*% V_inv   # n×n
  a           <- as.numeric(lambda * phylomat[species, species] %*% p)
  names(a)    <- species

  # Correct LOO adjustment for each focal (diagonal of Type 2 adj_all)
  adj_correct <- a - lambda * w * diag(E)

  rows <- list()

  for (rank in names(clade_tables)) {
    cat("  Rank:", rank, "\n")
    clades <- clade_tables[[rank]]

    for (clade_name in names(clades)) {
      cl     <- clades[[clade_name]]
      depth_M <- cl$depth_M

      # Build cross-cov vector cov_M (raw, not lambda-scaled)
      rep_sp   <- cl$rep_sp
      cov_M    <- phylomat[rep_sp, species]     # non-descendant covariances
      cov_M[cl$desc_train] <- depth_M           # descendants share history to M

      v_M_lam  <- lambda * cov_M               # lambda-scaled

      # Full adjustment (no LOO)
      adj_full <- as.numeric(t(v_M_lam) %*% p)

      # LOO correction vector: c_M = V_inv %*% v_M_lam
      # We only need element [i] for each focal i in this clade
      c_M_i <- as.numeric(V_inv %*% v_M_lam)

      # One row per focal species in this clade
      for (sp in cl$in_training) {
        i           <- sp_idx[sp]
        adj_loo_i   <- adj_full - w[i] * c_M_i[i]
        error_b_i   <- adj_loo_i - adj_correct[i]
        error_a_i   <- adj_loo_i - resid[i]

        rows[[length(rows) + 1L]] <- data.frame(
          focal      = sp,
          rank       = rank,
          clade      = clade_name,
          depth_M    = depth_M,
          adj_loo    = adj_loo_i,
          adj_correct = adj_correct[i],
          error_b    = error_b_i,
          error_a    = error_a_i,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  result_df <- do.call(rbind, rows)
  result_df$target <- tname
  all_results[[tname]] <- result_df

  cat("   Rows:", nrow(result_df), "\n")
}

# ==============================================================================
# 5. MERGE TARGETS AND SAVE TABLE
# ==============================================================================

mat_df <- all_results$mat[, c("focal","rank","clade","depth_M","error_b","error_a")]
names(mat_df)[names(mat_df) == "error_b"] <- "error_b_mat"
names(mat_df)[names(mat_df) == "error_a"] <- "error_a_mat"

map_df <- all_results$map[, c("focal","rank","clade","depth_M","error_b","error_a")]
names(map_df)[names(map_df) == "error_b"] <- "error_b_map"
names(map_df)[names(map_df) == "error_a"] <- "error_a_map"

result_df <- merge(mat_df, map_df,
                   by = c("focal","rank","clade","depth_M"),
                   all = TRUE)

# Rank as ordered factor for plotting
result_df$rank <- factor(result_df$rank,
                         levels = c("genus","family","order"),
                         ordered = TRUE)

write.csv(result_df, "tables/type1_degradation.csv", row.names = FALSE)
cat("Saved tables/type1_degradation.csv (", nrow(result_df), "rows)\n")

# ==============================================================================
# 6. SUMMARY TABLE — mean errors by rank
# ==============================================================================

summary_df <- do.call(rbind, lapply(split(result_df, result_df$rank), function(d) {
  data.frame(
    rank             = as.character(d$rank[1]),
    n_focal          = nrow(d),
    mean_depth_M     = mean(d$depth_M),
    mean_abs_eb_mat  = mean(abs(d$error_b_mat), na.rm = TRUE),
    sd_eb_mat        = sd(d$error_b_mat, na.rm = TRUE),
    mean_abs_ea_mat  = mean(abs(d$error_a_mat), na.rm = TRUE),
    mean_abs_eb_map  = mean(abs(d$error_b_map), na.rm = TRUE),
    sd_eb_map        = sd(d$error_b_map, na.rm = TRUE),
    mean_abs_ea_map  = mean(abs(d$error_a_map), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_df, "tables/type1_degradation_summary.csv", row.names = FALSE)
cat("Saved tables/type1_degradation_summary.csv\n")
print(summary_df)

# ==============================================================================
# 7. FIGURES
# ==============================================================================

rank_colours <- c(genus = "#1b9e77", family = "#d95f02", order = "#7570b3")

# sd(resid) asymptotes for reference lines
sd_mat <- sd(targets$mat$resid)
sd_map <- sd(targets$map$resid)

# Thin scatter for background points
set.seed(42)
plot_idx <- sample(nrow(result_df), min(20000, nrow(result_df)))
plot_df  <- result_df[plot_idx, ]

clip99 <- function(x) {
  q <- quantile(abs(x), 0.99, na.rm = TRUE)
  pmax(pmin(x, q), -q)
}

plot_df$error_b_mat_c <- clip99(plot_df$error_b_mat)
plot_df$error_b_map_c <- clip99(plot_df$error_b_map)
plot_df$error_a_mat_c <- clip99(plot_df$error_a_mat)
plot_df$error_a_map_c <- clip99(plot_df$error_a_map)

# Per-rank summary bands (median ± IQR) for overlay
make_rank_bands <- function(df, y_col) {
  do.call(rbind, lapply(split(df, df$rank), function(d) {
    y <- d[[y_col]]
    data.frame(
      rank   = as.character(d$rank[1]),
      median = median(y, na.rm = TRUE),
      q25    = quantile(y, 0.25, na.rm = TRUE),
      q75    = quantile(y, 0.75, na.rm = TRUE),
      depth  = mean(d$depth_M),
      stringsAsFactors = FALSE
    )
  }))
}

make_error_b_panel <- function(scatter_df, scatter_col, label, sd_val) {
  bands <- make_rank_bands(result_df, sub("_c$", "", scatter_col))
  ggplot() +
    geom_jitter(data = scatter_df,
                aes(x = depth_M, y = .data[[scatter_col]], colour = rank),
                alpha = 0.06, size = 0.4, width = 0.3) +
    geom_pointrange(data = bands,
                    aes(x = depth, y = median, ymin = q25, ymax = q75,
                        colour = rank),
                    size = 0.7, linewidth = 1.2) +
    geom_hline(yintercept =  sd_val, linetype = "dashed",
               colour = "darkorange", linewidth = 0.7) +
    geom_hline(yintercept = -sd_val, linetype = "dashed",
               colour = "darkorange", linewidth = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = rank_colours, name = "Rank") +
    labs(x = "Attachment depth (Ma)",
         y = paste("Marginal error (b) —", label),
         caption = "Point: median | bar: IQR | dashed: ±sd(resid)") +
    theme_bw(base_size = 11)
}

make_error_a_panel <- function(scatter_df, scatter_col, label) {
  bands <- make_rank_bands(result_df, sub("_c$", "", scatter_col))
  ggplot() +
    geom_jitter(data = scatter_df,
                aes(x = depth_M, y = .data[[scatter_col]], colour = rank),
                alpha = 0.06, size = 0.4, width = 0.3) +
    geom_pointrange(data = bands,
                    aes(x = depth, y = median, ymin = q25, ymax = q75,
                        colour = rank),
                    size = 0.7, linewidth = 1.2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = rank_colours, name = "Rank") +
    labs(x = "Attachment depth (Ma)",
         y = paste("Absolute prediction error (a) —", label),
         caption = "Point: median | bar: IQR") +
    theme_bw(base_size = 11)
}

cat("Generating figures...\n")

p1 <- make_error_b_panel(plot_df, "error_b_mat_c", "MAT (°C)", sd_mat)
p2 <- make_error_b_panel(plot_df, "error_b_map_c", "log(MAP)", sd_map)
p3 <- make_error_a_panel(plot_df, "error_a_mat_c", "MAT (°C)")
p4 <- make_error_a_panel(plot_df, "error_a_map_c", "log(MAP)")

fig1 <- arrangeGrob(p1, p2, ncol = 2,
           top = "Type 1 degradation: marginal PIP adjustment error vs attachment depth")
ggsave("plots/type1_error_b.pdf", fig1, width = 12, height = 5)
ggsave("plots/type1_error_b.png", fig1, width = 12, height = 5, dpi = 150)
cat("Saved plots/type1_error_b.{pdf,png}\n")

fig2 <- arrangeGrob(p3, p4, ncol = 2,
           top = "Type 1 degradation: absolute climate prediction error vs attachment depth")
ggsave("plots/type1_error_a.pdf", fig2, width = 12, height = 5)
ggsave("plots/type1_error_a.png", fig2, width = 12, height = 5, dpi = 150)
cat("Saved plots/type1_error_a.{pdf,png}\n")

# ==============================================================================
# 8. SUMMARY DIAGNOSTICS
# ==============================================================================

cat("\n=== SUMMARY ===\n")
cat("MAT sd(resid):", round(sd_mat, 4), "\n")
cat("MAP sd(resid):", round(sd_map, 4), "\n")
cat("\nMean |error_b| by rank:\n")
print(summary_df[, c("rank","n_focal","mean_depth_M",
                      "mean_abs_eb_mat","mean_abs_eb_map")])

cat("\nDone.\n")
