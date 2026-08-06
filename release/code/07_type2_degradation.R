source("code/_setup.R")

library(ape)
library(ggplot2)
library(gridExtra)

# 1. LOAD PIP COMPONENTS (impute variant, backward-compatible top-level keys)

pip <- readRDS("models/pip_components.rds")

phylomat   <- vcv(pip$tree_pruned)
diag(phylomat) <- diag(phylomat) + 1e-6

species    <- rownames(pip$V_lam_mat)    # same order for mat and map
n          <- length(species)

cat("n species:", n, "\n")

# Per-target components (impute variant)
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

# 2. REUSABLE HELPER — maps a covariance vector to its PIP adjustment
#    Consumed by the Type 1 (taxonomic imprecision) analysis in 08_.

# Phylogenetic adjustment for one new placement. v_cross is the lambda-scaled
# cross-covariance to every training species; V_inv and resid share its order.
pip_adj_from_vcov <- function(v_cross, V_inv, resid) {
  as.numeric(t(v_cross) %*% V_inv %*% resid)
}

# 3. Vectorised leave-one-out adjustment: adj(i, R) = a[R] - lambda * w[i] * E[R,i],
#    with p = V_inv %*% resid, a = lambda * phylomat %*% p, E = phylomat %*% V_inv,
#    w = p / diag(V_inv).

cat("Computing LOO adjustment matrices...\n")

results <- list()

for (tname in names(targets)) {
  cat(" ", tname, "...\n")
  comp   <- targets[[tname]]
  V_lam  <- comp$V_lam[species, species]
  resid  <- comp$resid
  lambda <- comp$lambda

  # O(n^3) LAPACK — once per target
  V_inv  <- solve(V_lam)

  # O(n^2) matrix-vector products
  p      <- as.numeric(V_inv %*% resid)
  names(p) <- species

  # O(n^3) BLAS — once per target; E[R,i] = (phylomat %*% V_inv)[R,i]
  E      <- phylomat[species, species] %*% V_inv

  # Full adjustment vector (no LOO) — a[R] is PIP adj for a fossil placed at R
  a      <- as.numeric(lambda * phylomat[species, species] %*% p)
  names(a) <- species

  # Per-species leverage weight
  w      <- p / diag(V_inv)

  # n×n LOO adjustment matrix: adj_all[i, R] = adj for focal i placed at R
  # adj_all[i, R] = a[R] - lambda * w[i] * E[R, i]
  adj_all <- matrix(a, nrow = n, ncol = n, byrow = TRUE) -
             lambda * outer(w, rep(1, n)) * t(E)
  rownames(adj_all) <- species   # focal species
  colnames(adj_all) <- species   # recipient (placement) species

  # Correct adjustment for each focal: diagonal of adj_all
  adj_correct <- diag(adj_all)

  # Error (b): marginal — how much wrong placement shifts the PIP adjustment
  # error_b[i, R] = adj_all[i, R] - adj_all[i, i]
  error_b <- sweep(adj_all, 1, adj_correct, "-")

  # Error (a): absolute error vs true climate, which reduces to
  # adj_all[i, R] - resid[i] because the Xbeta term cancels.
  error_a <- sweep(adj_all, 1, resid, "-")

  results[[tname]] <- list(
    adj_all     = adj_all,
    adj_correct = adj_correct,
    error_b     = error_b,
    error_a     = error_a,
    resid       = resid,
    sd_resid    = sd(resid),
    V_inv       = V_inv,
    label       = comp$label
  )

  cat("    sd(resid):", round(sd(resid), 4),
      " | adj_correct range:", round(range(adj_correct), 3), "\n")
}

# 4. COPHENETIC DISTANCE MATRIX

cat("Computing cophenetic distances...\n")
coph <- cophenetic(pip$tree_pruned)[species, species]

# 5. BUILD PAIRWISE RESULTS TABLE

cat("Assembling pairwise table (n*(n-1) =", n*(n-1), "rows)...\n")

# Use lower + upper triangle (directed: focal → recipient)
idx <- which(row(results$mat$error_b) != col(results$mat$error_b), arr.ind = TRUE)

pair_df <- data.frame(
  focal       = species[idx[, 1]],
  recipient   = species[idx[, 2]],
  coph_dist   = coph[idx],
  error_b_mat = results$mat$error_b[idx],
  error_a_mat = results$mat$error_a[idx],
  error_b_map = results$map$error_b[idx],
  error_a_map = results$map$error_a[idx],
  stringsAsFactors = FALSE
)

cat("Table rows:", nrow(pair_df), "\n")
write.csv(pair_df, "tables/type2_degradation_pairwise.csv", row.names = FALSE)
cat("Saved tables/type2_degradation_pairwise.csv\n")

# 6. SUMMARY TABLE — binned by cophenetic distance

pair_df$dist_bin <- cut(pair_df$coph_dist,
                        breaks = unique(quantile(pair_df$coph_dist,
                                                 probs = seq(0, 1, by = 0.05))),
                        include.lowest = TRUE)

summary_df <- do.call(rbind, lapply(split(pair_df, pair_df$dist_bin), function(d) {
  data.frame(
    dist_bin        = as.character(d$dist_bin[1]),
    n_pairs         = nrow(d),
    mean_coph_dist  = mean(d$coph_dist),
    mean_abs_eb_mat = mean(abs(d$error_b_mat)),
    sd_eb_mat       = sd(d$error_b_mat),
    mean_abs_ea_mat = mean(abs(d$error_a_mat)),
    mean_abs_eb_map = mean(abs(d$error_b_map)),
    sd_eb_map       = sd(d$error_b_map),
    mean_abs_ea_map = mean(abs(d$error_a_map))
  )
}))

write.csv(summary_df, "tables/type2_degradation_summary.csv", row.names = FALSE)
cat("Saved tables/type2_degradation_summary.csv\n")

# 7. FIGURES

# Pre-compute per-bin quantiles (10th, 50th, 90th) from all 3.4M pairs.
# log_x = TRUE bins on log10 scale for equal visual spacing at short distances.
make_quantile_bands <- function(df, y_col, n_bins = 40, log_x = FALSE) {
  x <- if (log_x) log10(pmax(df$coph_dist, 1e-3)) else df$coph_dist
  breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1)))
  df$bin <- cut(x, breaks = breaks, include.lowest = TRUE)
  do.call(rbind, lapply(split(df, df$bin), function(d) {
    y <- d[[y_col]]
    x_d <- if (log_x) log10(pmax(d$coph_dist, 1e-3)) else d$coph_dist
    data.frame(
      mid  = mean(x_d),
      q10  = quantile(y, 0.10),
      q50  = quantile(y, 0.50),
      q90  = quantile(y, 0.90),
      stringsAsFactors = FALSE
    )
  }))
}

cat("Computing quantile bands...\n")
qbands <- list(
  error_b_mat = make_quantile_bands(pair_df, "error_b_mat"),
  error_b_map = make_quantile_bands(pair_df, "error_b_map"),
  error_a_mat = make_quantile_bands(pair_df, "error_a_mat"),
  error_a_map = make_quantile_bands(pair_df, "error_a_map")
)
qbands_log <- list(
  error_b_mat = make_quantile_bands(pair_df, "error_b_mat", log_x = TRUE),
  error_b_map = make_quantile_bands(pair_df, "error_b_map", log_x = TRUE),
  error_a_mat = make_quantile_bands(pair_df, "error_a_mat", log_x = TRUE),
  error_a_map = make_quantile_bands(pair_df, "error_a_map", log_x = TRUE)
)

# Thin scatter for background points only
set.seed(42)
plot_idx <- sample(nrow(pair_df), min(50000, nrow(pair_df)))
plot_df  <- pair_df[plot_idx, ]

# Clip extreme outliers for scatter layer (keep 99th percentile)
clip99 <- function(x) {
  q <- quantile(abs(x), 0.99, na.rm = TRUE)
  pmax(pmin(x, q), -q)
}
plot_df$error_b_mat_c <- clip99(plot_df$error_b_mat)
plot_df$error_b_map_c <- clip99(plot_df$error_b_map)
plot_df$error_a_mat_c <- clip99(plot_df$error_a_mat)
plot_df$error_a_map_c <- clip99(plot_df$error_a_map)

make_error_b_plot <- function(scatter_df, scatter_col, bands, label, sd_resid,
                             log_x = FALSE) {
  x_col <- if (log_x) "log10_coph_dist" else "coph_dist"
  x_lab <- if (log_x) "log10 cophenetic distance (Ma)" else "Cophenetic distance (Ma)"
  ggplot() +
    geom_point(data = scatter_df,
               aes(x = .data[[x_col]], y = .data[[scatter_col]]),
               alpha = 0.03, size = 0.3, colour = "steelblue") +
    geom_ribbon(data = bands,
                aes(x = mid, ymin = q10, ymax = q90),
                fill = "firebrick", alpha = 0.20) +
    geom_line(data = bands,
              aes(x = mid, y = q50),
              colour = "firebrick", linewidth = 0.9) +
    geom_hline(yintercept =  sd_resid, linetype = "dashed",
               colour = "darkorange", linewidth = 0.7) +
    geom_hline(yintercept = -sd_resid, linetype = "dashed",
               colour = "darkorange", linewidth = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    labs(x = x_lab,
         y = paste("Marginal error (b) —", label),
         caption = "Red line: median  |  band: 10th–90th percentile  |  dashed: ±sd(resid)") +
    theme_bw(base_size = 11)
}

make_error_a_plot <- function(scatter_df, scatter_col, bands, label,
                              log_x = FALSE) {
  x_col <- if (log_x) "log10_coph_dist" else "coph_dist"
  x_lab <- if (log_x) "log10 cophenetic distance (Ma)" else "Cophenetic distance (Ma)"
  ggplot() +
    geom_point(data = scatter_df,
               aes(x = .data[[x_col]], y = .data[[scatter_col]]),
               alpha = 0.03, size = 0.3, colour = "mediumpurple") +
    geom_ribbon(data = bands,
                aes(x = mid, ymin = q10, ymax = q90),
                fill = "firebrick", alpha = 0.20) +
    geom_line(data = bands,
              aes(x = mid, y = q50),
              colour = "firebrick", linewidth = 0.9) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    labs(x = x_lab,
         y = paste("Absolute prediction error (a) —", label),
         caption = "Red line: median  |  band: 10th–90th percentile") +
    theme_bw(base_size = 11)
}

cat("Generating figures...\n")

plot_df$log10_coph_dist <- log10(pmax(plot_df$coph_dist, 1e-3))

p1 <- make_error_b_plot(plot_df, "error_b_mat_c", qbands$error_b_mat,
                        "MAT (°C)", results$mat$sd_resid)
p2 <- make_error_b_plot(plot_df, "error_b_map_c", qbands$error_b_map,
                        "log(MAP)", results$map$sd_resid)
p3 <- make_error_a_plot(plot_df, "error_a_mat_c", qbands$error_a_mat, "MAT (°C)")
p4 <- make_error_a_plot(plot_df, "error_a_map_c", qbands$error_a_map, "log(MAP)")

p1l <- make_error_b_plot(plot_df, "error_b_mat_c", qbands_log$error_b_mat,
                         "MAT (°C)", results$mat$sd_resid, log_x = TRUE)
p2l <- make_error_b_plot(plot_df, "error_b_map_c", qbands_log$error_b_map,
                         "log(MAP)", results$map$sd_resid, log_x = TRUE)
p3l <- make_error_a_plot(plot_df, "error_a_mat_c", qbands_log$error_a_mat,
                         "MAT (°C)", log_x = TRUE)
p4l <- make_error_a_plot(plot_df, "error_a_map_c", qbands_log$error_a_map,
                         "log(MAP)", log_x = TRUE)

# Figure 1: marginal error (b) — linear distance
fig1 <- arrangeGrob(p1, p2, ncol = 2,
                    top = "Type 2 degradation: marginal PIP adjustment error vs cophenetic distance")
ggsave("plots/type2_error_b.pdf", fig1, width = 12, height = 5)
ggsave("plots/type2_error_b.png", fig1, width = 12, height = 5, dpi = 150)
cat("Saved plots/type2_error_b.{pdf,png}\n")

# Figure 1 (log x): marginal error (b) — log10 distance
fig1l <- arrangeGrob(p1l, p2l, ncol = 2,
                     top = "Type 2 degradation: marginal error vs log10 cophenetic distance")
ggsave("plots/type2_error_b_log10.pdf", fig1l, width = 12, height = 5)
ggsave("plots/type2_error_b_log10.png", fig1l, width = 12, height = 5, dpi = 150)
cat("Saved plots/type2_error_b_log10.{pdf,png}\n")

# Figure 2: absolute prediction error (a) — linear distance
fig2 <- arrangeGrob(p3, p4, ncol = 2,
                    top = "Type 2 degradation: absolute climate prediction error vs cophenetic distance")
ggsave("plots/type2_error_a.pdf", fig2, width = 12, height = 5)
ggsave("plots/type2_error_a.png", fig2, width = 12, height = 5, dpi = 150)
cat("Saved plots/type2_error_a.{pdf,png}\n")

# Figure 2 (log x): absolute prediction error (a) — log10 distance
fig2l <- arrangeGrob(p3l, p4l, ncol = 2,
                     top = "Type 2 degradation: absolute error vs log10 cophenetic distance")
ggsave("plots/type2_error_a_log10.pdf", fig2l, width = 12, height = 5)
ggsave("plots/type2_error_a_log10.png", fig2l, width = 12, height = 5, dpi = 150)
cat("Saved plots/type2_error_a_log10.{pdf,png}\n")

# 8. SUMMARY DIAGNOSTICS

cat("\n=== SUMMARY ===\n")
cat("MAT sd(resid):", round(results$mat$sd_resid, 4), "\n")
cat("MAP sd(resid):", round(results$map$sd_resid, 4), "\n")
q10 <- quantile(pair_df$coph_dist, 0.10)
q90 <- quantile(pair_df$coph_dist, 0.90)
cat("Cophenetic dist range:", round(range(pair_df$coph_dist), 2),
    "| 10th pct:", round(q10, 2), "| 90th pct:", round(q90, 2), "\n")
cat("  Closest 10% (dist <", round(q10, 1), "): mean|error_b_mat| =",
    round(mean(abs(pair_df$error_b_mat[pair_df$coph_dist <= q10])), 4), "\n")
cat("  Farthest 10% (dist >", round(q90, 1), "): mean|error_b_mat| =",
    round(mean(abs(pair_df$error_b_mat[pair_df$coph_dist >= q90])), 4), "\n")
cat("  Closest 10%: mean|error_b_map| =",
    round(mean(abs(pair_df$error_b_map[pair_df$coph_dist <= q10])), 4), "\n")
cat("  Farthest 10%: mean|error_b_map| =",
    round(mean(abs(pair_df$error_b_map[pair_df$coph_dist >= q90])), 4), "\n")

cat("\nDone.\n")
