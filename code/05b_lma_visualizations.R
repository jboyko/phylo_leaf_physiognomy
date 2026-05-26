setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ggplot2)
library(dplyr)
library(tidyr)

# ── 1. Load CSVs ───────────────────────────────────────────────────────────────
rmse  <- read.csv("tables/lma_cv_rmse.csv",                  stringsAsFactors = FALSE)
preds <- read.csv("tables/lma_cv_site_predictions.csv",      stringsAsFactors = FALSE)
coefs <- read.csv("tables/lma_cv_model_coefs.csv",           stringsAsFactors = FALSE)
fit   <- read.csv("tables/lma_cv_model_fit.csv",             stringsAsFactors = FALSE)

# ── 2. Shared aesthetics ───────────────────────────────────────────────────────
model_labels <- c(PIP = "PIP", LM = "LM (Peppe)", PGLS = "PGLS")
model_colours <- c(PIP = "#1b7837", LM = "#d6604d", PGLS = "#9970ab")

# ── 3. Figure 1: RMSE dot plot ────────────────────────────────────────────────
rmse_plot <- rmse |>
  mutate(
    label = factor(model_labels[model], levels = rev(model_labels)),
    colour = model_colours[model]
  )

p1 <- ggplot(rmse_plot, aes(x = rmse, y = label, colour = model)) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.4f", rmse)), hjust = -0.25, size = 3.5) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_x_continuous(
    limits = c(min(rmse$rmse) * 0.97, max(rmse$rmse) * 1.05),
    name   = "10-fold CV RMSE (log10 LMA)"
  ) +
  labs(
    title = "LMA model comparison — 10-fold CV RMSE",
    y     = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank()
  )

ggsave("plots/lma_fig1_rmse.pdf", p1, width = 6, height = 3)
ggsave("plots/lma_fig1_rmse.png", p1, width = 6, height = 3, dpi = 150)
message("Saved lma_fig1_rmse")

# ── 4. Figure 2: Observed vs. predicted ───────────────────────────────────────
preds_long <- preds |>
  pivot_longer(
    cols      = c(pred_lm, pred_pgls, pred_pip),
    names_to  = "model_raw",
    values_to = "predicted"
  ) |>
  mutate(
    model = case_when(
      model_raw == "pred_lm"   ~ "LM",
      model_raw == "pred_pgls" ~ "PGLS",
      model_raw == "pred_pip"  ~ "PIP"
    ),
    label = factor(model_labels[model], levels = model_labels)
  ) |>
  filter(!is.na(predicted))

lims <- range(c(preds_long$obs, preds_long$predicted), na.rm = TRUE)

p2 <- ggplot(preds_long, aes(x = obs, y = predicted, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  facet_wrap(~ label, nrow = 1) +
  scale_colour_manual(values = model_colours, guide = "none") +
  coord_fixed(xlim = lims, ylim = lims) +
  labs(
    title = "Observed vs. predicted — 10-fold CV (site level)",
    x     = "Observed log10(LMA)",
    y     = "Predicted log10(LMA)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("plots/lma_fig2_obs_vs_pred.pdf", p2, width = 10, height = 4)
ggsave("plots/lma_fig2_obs_vs_pred.png", p2, width = 10, height = 4, dpi = 150)
message("Saved lma_fig2_obs_vs_pred")

# ── 5. Figure 3: Slope stability across folds ─────────────────────────────────
coef_slope <- coefs |>
  filter(predictor == "log10_petiole_metric") |>
  mutate(
    label = factor(model_labels[method], levels = model_labels[c("LM", "PGLS")])
  )

p3 <- ggplot(coef_slope, aes(x = label, y = estimate, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(aes(fill = method), alpha = 0.4, width = 0.4, outlier.size = 1.5) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_fill_manual(values   = model_colours, guide = "none") +
  labs(
    title = "Slope stability across CV folds — log10(PW2/A) coefficient",
    x     = NULL,
    y     = "Coefficient estimate"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("plots/lma_fig3_slope_stability.pdf", p3, width = 5, height = 4.5)
ggsave("plots/lma_fig3_slope_stability.png", p3, width = 5, height = 4.5, dpi = 150)
message("Saved lma_fig3_slope_stability")

# ── 6. Figure 4: Lambda across folds (PGLS) ───────────────────────────────────
lambda_df <- fit |> filter(method == "PGLS", !is.na(lambda))

p4 <- ggplot(lambda_df, aes(x = "PGLS", y = lambda)) +
  geom_boxplot(fill = model_colours["PGLS"], alpha = 0.5, width = 0.3, outlier.size = 1.5) +
  geom_jitter(colour = model_colours["PGLS"], width = 0.06, size = 2.5, alpha = 0.85) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = expression(paste("Pagel's ", lambda, " across CV folds — PGLS")),
    x     = NULL,
    y     = expression(paste("Pagel's ", lambda))
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave("plots/lma_fig4_lambda.pdf", p4, width = 4, height = 4.5)
ggsave("plots/lma_fig4_lambda.png", p4, width = 4, height = 4.5, dpi = 150)
message("Saved lma_fig4_lambda")

# ── 7. Diagnostics table ───────────────────────────────────────────────────────
diag_table <- preds_long |>
  filter(!is.na(obs), !is.na(predicted)) |>
  group_by(model, label) |>
  summarise(
    n            = n(),
    mean_bias    = mean(predicted - obs),
    bias_sq      = mean_bias^2,
    rmse         = sqrt(mean((predicted - obs)^2)),
    residual_var = rmse^2 - bias_sq,
    slope        = coef(lm(predicted ~ obs))[2],
    intercept    = coef(lm(predicted ~ obs))[1],
    .groups      = "drop"
  ) |>
  arrange(rmse)

write.csv(diag_table, "tables/lma_prediction_diagnostics.csv", row.names = FALSE)
message("Saved tables/lma_prediction_diagnostics.csv")
print(as.data.frame(diag_table), digits = 3)

message("\nAll LMA figures saved to plots/")
