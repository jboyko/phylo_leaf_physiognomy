# Generate the extant validation displays from the repository root.
library(dplyr)
library(tidyr)
library(ggplot2)

submission <- "submission/submission_1"
climate <- read.csv("tables/loso_cv_site_predictions.csv")
lma <- read.csv("tables/lma_cv_site_predictions.csv")
paired <- read.csv("tables/cv_paired_rmse.csv")

# Use the zero-filled site LM that Dana requested for the fossil comparison.
climate_models <- c(
  pip_sp_site_impute = "PIP (imputed)",
  pip_sp_site_cc = "PIP (complete case)",
  dilp_cv_site = "Refitted DiLP",
  lm_site_site_sp_zero_cc = "Site LM (complete case)",
  lm_site_site_sp_zero_impute = "Site LM (imputed)",
  lm_sp_site_cc = "Species LM (complete case)",
  lm_sp_site_impute = "Species LM (imputed)",
  pgls_sp_site_impute = "PGLS (imputed)",
  pgls_sp_site_cc = "PGLS (complete case)"
)
supplementary_site_models <- c(
  lm_site_site_untoothed_excl_impute = "Untoothed-excluded site LM (imputed)",
  lm_site_site_untoothed_excl_cc = "Untoothed-excluded site LM (complete case)"
)
validation_models <- c(climate_models, supplementary_site_models)
targets <- c(mat = "MAT (°C)", log_map = "ln(MAP in cm)")

climate_predictions <- bind_rows(lapply(names(targets), function(target) {
  columns <- paste0(names(validation_models), "_", target)
  climate |>
    select(site, fold, all_of(columns)) |>
    mutate(observed = climate[[paste0("obs_", target)]]) |>
    pivot_longer(all_of(columns), names_to = "column", values_to = "predicted") |>
    mutate(
      model = unname(validation_models[match(column, columns)]),
      target = target
    ) |>
    select(site, fold, target, model, observed, predicted)
})) |>
  filter(is.finite(observed), is.finite(predicted))

lma_predictions <- lma |>
  select(site, observed = obs, LM = pred_uni_lm, PGLS = pred_uni_pgls,
         PIP = pred_uni_pip) |>
  drop_na() |>
  pivot_longer(c(LM, PGLS, PIP), names_to = "model", values_to = "predicted")

diagnostics <- function(predictions) {
  predictions |>
    mutate(error = predicted - observed) |>
    summarise(
      n_sites = n(),
      rmse = sqrt(mean(error^2)),
      slope = unname(coef(lm(predicted ~ observed))[2]),
      mean_bias = mean(error),
      bias_squared = mean_bias^2,
      residual_variance = mean((error - mean_bias)^2),
      .groups = "drop"
    )
}

climate_stats <- climate_predictions |> group_by(target, model) |> diagnostics()
main_climate_stats <- filter(climate_stats, model %in% unname(climate_models))
lma_stats <- lma_predictions |> group_by(model) |> diagnostics() |> arrange(rmse)

climate_rmse <- main_climate_stats |>
  select(model, target, rmse, n_sites) |>
  pivot_wider(names_from = target, values_from = c(rmse, n_sites),
              names_glue = "{target}_{.value}") |>
  arrange(mat_rmse) |>
  select(model, mat_rmse, mat_n_sites, log_map_rmse, log_map_n_sites)

climate_diagnostics <- main_climate_stats |>
  filter(!grepl("complete case", model, fixed = TRUE)) |>
  arrange(rmse) |>
  select(target, model, n_sites, slope, mean_bias, bias_squared, residual_variance, rmse)

write_table <- function(data, number) {
  write.csv(data, file.path(submission, "tables", paste0("Table_", number, ".csv")),
            row.names = FALSE)
}

write_table(climate_rmse, "1")
write_table(filter(climate_diagnostics, target == "mat") |> select(-target), "2")
write_table(filter(climate_diagnostics, target == "log_map") |> select(-target), "3")
write_table(select(lma_stats, model, rmse, n_sites), "4")
write_table(select(lma_stats, model, n_sites, slope, mean_bias, bias_squared,
                   residual_variance, rmse), "5")
write_table(paired, "S1")
write_table(
  climate_stats |>
    filter(model %in% unname(supplementary_site_models)) |>
    arrange(match(target, names(targets)), match(model, unname(supplementary_site_models))) |>
    select(target, model, n_sites, rmse, slope, mean_bias, bias_squared, residual_variance),
  "S2"
)

accent <- "#176B58"
baseline <- "#546673"
ink <- "#172832"
figure_theme <- theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E1E7E9", linewidth = 0.3),
    panel.spacing = grid::unit(15, "pt"),
    strip.text = element_text(face = "bold", size = 10, margin = margin(b = 8)),
    legend.position = "none",
    plot.margin = margin(9, 12, 9, 9)
  )

figure_models <- unname(climate_models[c(
  "pip_sp_site_impute", "lm_site_site_sp_zero_impute",
  "lm_sp_site_impute", "dilp_cv_site"
)])
scatter <- climate_predictions |>
  filter(model %in% figure_models) |>
  group_by(site, target) |>
  filter(n() == length(figure_models)) |>
  ungroup() |>
  mutate(
    model = factor(model, levels = figure_models),
    target = factor(targets[target], levels = unname(targets))
  )

scatter_labels <- setNames(c("PIP\n(imputed)", "Site LM\n(12 traits)",
                             "Species LM\n(imputed)", "Refitted DiLP"), figure_models)

scatter_row <- function(target_label) {
  points <- filter(scatter, target == target_label)
  limits <- range(c(points$observed, points$predicted))
  n_sites <- n_distinct(points$site)

  ggplot(points, aes(observed, predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "#87969E", linewidth = 0.4) +
    geom_point(aes(colour = model == figure_models[1]), alpha = 0.67, size = 1.2) +
    geom_smooth(aes(colour = model == figure_models[1]), method = "lm",
                formula = y ~ x, fill = "#A1B6BE", alpha = 0.25, linewidth = 0.65) +
    facet_wrap(vars(model), nrow = 1, labeller = as_labeller(scatter_labels)) +
    scale_colour_manual(values = c(baseline, accent)) +
    coord_equal(xlim = limits, ylim = limits) +
    labs(title = sprintf("%s  (%d sites)", target_label, n_sites),
         x = "Observed", y = "Predicted") +
    figure_theme +
    theme(plot.title = element_text(face = "bold", size = 11, margin = margin(b = 7)))
}

figure_s1 <- gridExtra::arrangeGrob(grobs = lapply(unname(targets), scatter_row), ncol = 1)

ggsave(file.path(submission, "figures/Figure_S1.pdf"), figure_s1,
       device = cairo_pdf, width = 8.5, height = 6)
write.csv(scatter, file.path(submission, "figures/Figure_S1_data.csv"), row.names = FALSE)
