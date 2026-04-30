setwd("/Users/jboyko/phylo_leaf_physiognomy")

library(ggplot2)
library(dplyr)
library(tidyr)

# ── 1. Load CSVs ──────────────────────────────────────────────────────────────
rmse    <- read.csv("tables/loso_cv_rmse.csv", stringsAsFactors = FALSE)
preds   <- read.csv("tables/loso_cv_site_predictions.csv", stringsAsFactors = FALSE)
coefs   <- read.csv("tables/loso_cv_model_coefs.csv", stringsAsFactors = FALSE)
fit     <- read.csv("tables/loso_cv_model_fit.csv", stringsAsFactors = FALSE)

# ── 2. Shared label lookup ────────────────────────────────────────────────────
model_labels <- c(
  pip_sp_site_impute              = "PIP (impute)",
  pip_sp_site_cc                  = "PIP (CC)",
  lm_site_site_peppe_impute       = "LM site, Peppe (impute)",
  lm_site_site_peppe_cc           = "LM site, Peppe (CC)",
  lm_site_site_sp_zero_impute     = "LM site, sp+zero (impute)",
  lm_site_site_sp_zero_cc         = "LM site, sp+zero (CC)",
  lm_site_site_specimen_impute    = "LM site, specimen (impute)",
  lm_site_site_specimen_cc        = "LM site, specimen (CC)",
  lm_sp_site_impute               = "LM species (impute)",
  lm_sp_site_cc                   = "LM species (CC)",
  pgls_sp_site_impute             = "PGLS (impute)",
  pgls_sp_site_cc                 = "PGLS (CC)"
)

model_family <- c(
  pip_sp_site_impute              = "PIP",
  pip_sp_site_cc                  = "PIP",
  lm_site_site_peppe_impute       = "LM site",
  lm_site_site_peppe_cc           = "LM site",
  lm_site_site_sp_zero_impute     = "LM site",
  lm_site_site_sp_zero_cc         = "LM site",
  lm_site_site_specimen_impute    = "LM site",
  lm_site_site_specimen_cc        = "LM site",
  lm_sp_site_impute               = "LM species",
  lm_sp_site_cc                   = "LM species",
  pgls_sp_site_impute             = "PGLS",
  pgls_sp_site_cc                 = "PGLS"
)

family_colours <- c(
  PIP        = "#1b7837",
  "LM site"  = "#4393c3",
  "LM species" = "#d6604d",
  PGLS       = "#9970ab"
)

# ── 3. Figure 1: RMSE dot plot ────────────────────────────────────────────────
# Strip the trailing _mat / _log_map suffix to get the model key
rmse <- rmse |>
  mutate(
    model_key = sub("_(mat|log_map)$", "", column),
    label     = model_labels[model_key],
    family    = model_family[model_key],
    target_lab = ifelse(target == "mat", "MAT (°C)", "log(MAP) (log cm)")
  )

# Order models by mean RMSE across targets (worst → best so best is at top)
order_key <- rmse |>
  group_by(model_key) |>
  summarise(mean_rmse = mean(rmse), .groups = "drop") |>
  arrange(desc(mean_rmse)) |>
  pull(model_key)

rmse$label <- factor(rmse$label, levels = model_labels[order_key])

p1 <- ggplot(rmse, aes(x = rmse, y = label, colour = family)) +
  geom_point(size = 3.5) +
  facet_wrap(~ target_lab, scales = "free_x") +
  scale_colour_manual(values = family_colours, name = "Model family") +
  labs(
    title = "LOSO-CV RMSE across model configurations",
    x     = "RMSE",
    y     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank()
  )

ggsave("plots/fig1_rmse_comparison.pdf", p1, width = 10, height = 6)
ggsave("plots/fig1_rmse_comparison.png", p1, width = 10, height = 6, dpi = 150)
message("Saved fig1_rmse_comparison")

# ── 4. Figure 2: Observed vs. Predicted (impute only) ─────────────────────────
focal_models <- c(
  "pip_sp_site_impute",
  "lm_site_site_peppe_impute",
  "lm_sp_site_impute"
)

preds_long <- preds |>
  pivot_longer(
    cols      = -c(site, fold, obs_mat, obs_log_map),
    names_to  = "col",
    values_to = "predicted"
  ) |>
  mutate(
    target    = ifelse(grepl("_mat$", col), "mat", "log_map"),
    model_key = sub("_(mat|log_map)$", "", col),
    observed  = ifelse(target == "mat", obs_mat, obs_log_map)
  ) |>
  filter(model_key %in% focal_models) |>
  mutate(
    label      = factor(model_labels[model_key], levels = model_labels[focal_models]),
    target_lab = ifelse(target == "mat", "MAT (°C)", "log(MAP) (log cm)"),
    family     = model_family[model_key]
  )

# Force consistent axes within each variable column across model rows.
# geom_blank points at (lo, lo) and (hi, hi) per target × model combination
# make facet_grid(scales = "free") lock all MAT panels to the same range and
# all log(MAP) panels to their own shared range.
axis_lims <- preds_long |>
  group_by(target_lab) |>
  summarise(lo = min(c(observed, predicted), na.rm = TRUE),
            hi = max(c(observed, predicted), na.rm = TRUE),
            .groups = "drop")

blank_df <- preds_long |>
  distinct(label, target_lab, family) |>
  left_join(axis_lims, by = "target_lab") |>
  pivot_longer(cols = c(lo, hi), names_to = "lim", values_to = "val") |>
  mutate(observed = val, predicted = val)

p2 <- ggplot(preds_long, aes(x = observed, y = predicted, colour = family)) +
  geom_blank(data = blank_df) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  facet_grid(label ~ target_lab, scales = "free") +
  scale_colour_manual(values = family_colours, name = "Model family") +
  scale_fill_manual(values = family_colours, guide = "none") +
  labs(
    title = "Observed vs. predicted - LOSO-CV site-level predictions",
    x     = "Observed",
    y     = "Predicted"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text.y     = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

ggsave("plots/fig2_obs_vs_pred.pdf", p2, width = 8, height = 8)
ggsave("plots/fig2_obs_vs_pred.png", p2, width = 8, height = 8, dpi = 150)
message("Saved fig2_obs_vs_pred")

# ── 5. Figure 3: Coefficient stability across LOSO folds (LM species impute) ──
coef_lm_sp <- coefs |>
  filter(
    method      == "LM",
    train_level == "sp",
    config      == "grand_mean",
    na_handling == "impute",
    predictor   != "(Intercept)"
  ) |>
  mutate(target_lab = ifelse(target == "mat", "MAT (°C)", "log(MAP) (log cm)"))

# Rank predictors by median absolute estimate across targets
pred_order <- coef_lm_sp |>
  group_by(predictor) |>
  summarise(med_abs = median(abs(estimate)), .groups = "drop") |>
  arrange(med_abs) |>
  pull(predictor)

coef_lm_sp$predictor <- factor(coef_lm_sp$predictor, levels = pred_order)

p3 <- ggplot(coef_lm_sp, aes(x = estimate, y = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(aes(fill = target_lab), outlier.size = 1.2, width = 0.5,
               position = position_dodge(width = 0.6), alpha = 0.7) +
  scale_fill_manual(values = c("MAT (°C)" = "#d6604d", "log(MAP) (log cm)" = "#4393c3"),
                    name = NULL) +
  labs(
    title = "Coefficient stability across LOSO folds - LM species (impute)",
    x     = "Coefficient estimate",
    y     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("plots/fig3_coef_stability.pdf", p3, width = 8, height = 6)
ggsave("plots/fig3_coef_stability.png", p3, width = 8, height = 6, dpi = 150)
message("Saved fig3_coef_stability")

# ── 6. Figure 4: Lambda across folds (PGLS) ──────────────────────────────────
lambda_df <- fit |>
  filter(method == "PGLS", !is.na(lambda)) |>
  mutate(
    target_lab = ifelse(target == "mat", "MAT (°C)", "log(MAP) (log cm)"),
    na_label   = ifelse(na_handling == "impute", "Impute", "Complete-case")
  )

p4 <- ggplot(lambda_df, aes(x = na_label, y = lambda, fill = target_lab)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.75, width = 0.5,
               position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("MAT (°C)" = "#d6604d", "log(MAP) (log cm)" = "#4393c3"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Pagel's lambda across LOSO folds - PGLS models",
    x     = NULL,
    y     = expression(paste("Pagel's ", lambda))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("plots/fig4_lambda.pdf", p4, width = 5, height = 5)
ggsave("plots/fig4_lambda.png", p4, width = 5, height = 5, dpi = 150)
message("Saved fig4_lambda")

# ── 7. Table: Prediction diagnostics (slope, bias, RMSE decomposition) ────────
all_impute_models <- c(
  "pip_sp_site_impute",
  "lm_site_site_peppe_impute",
  "lm_site_site_sp_zero_impute",
  "lm_site_site_specimen_impute",
  "lm_sp_site_impute",
  "pgls_sp_site_impute"
)

diag_table <- preds |>
  pivot_longer(
    cols      = -c(site, fold, obs_mat, obs_log_map),
    names_to  = "col",
    values_to = "predicted"
  ) |>
  mutate(
    target    = ifelse(grepl("_mat$", col), "mat", "log_map"),
    model_key = sub("_(mat|log_map)$", "", col),
    observed  = ifelse(target == "mat", obs_mat, obs_log_map)
  ) |>
  filter(model_key %in% all_impute_models, !is.na(observed), !is.na(predicted)) |>
  group_by(model_key, target) |>
  summarise(
    n            = n(),
    mean_bias    = mean(predicted - observed),
    bias_sq      = mean_bias^2,
    rmse         = sqrt(mean((predicted - observed)^2)),
    residual_var = rmse^2 - bias_sq,
    slope        = coef(lm(predicted ~ observed))[2],
    intercept    = coef(lm(predicted ~ observed))[1],
    .groups = "drop"
  ) |>
  mutate(
    label      = model_labels[model_key],
    target_lab = ifelse(target == "mat", "MAT", "log(MAP)")
  ) |>
  arrange(target, rmse) |>
  select(label, target_lab, n, slope, intercept, mean_bias, bias_sq, residual_var, rmse)

write.csv(diag_table, "tables/prediction_diagnostics.csv", row.names = FALSE)
message("Saved tables/prediction_diagnostics.csv")
print(as.data.frame(diag_table), digits = 3)

message("\nAll figures saved to plots/")
