# From the repository root, run after producing the climate and LMA analyses.
source("code/setup.R")
library(caret)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
source("code/fossil_taxonomy.R")

submission <- "submission/submission_1"
figures <- file.path(submission, "figures")
tables <- file.path(submission, "tables")
climate_fit <- readRDS("models/pip_components.rds")$configs$impute
lma_fit <- readRDS("models/lma_pip_components.rds")$univariate
traits <- read.csv("data/fossil_traits.csv")
climate_predictions <- read.csv("tables/fossil_predictions_formal_only.csv")
climate_intervals <- read.csv("tables/fossil_site_uncertainty_formal_only.csv")
lma_predictions <- read.csv("tables/lma_fossil_species_predictions.csv")
lma_reference <- read.csv("tables/lma_fossil_site_predictions.csv")

climate_traits <- climate_predictions |>
  select(fossil_name, site, age_ma) |>
  left_join(traits, by = c("fossil_name", "site", "age_ma"),
            relationship = "one-to-one")
stopifnot(!anyNA(climate_traits$morphotype))

mean_on_scale <- function(values, transform, inverse) {
  inverse(mean(transform(values)))
}

site_points <- function(occurrences, mean_x, mean_y) {
  occurrences |>
    group_by(site, age_ma) |>
    summarise(x = mean_x(x), y = mean_y(y), n_occurrences = n(), .groups = "drop")
}

climate_panel <- function(fit, target, reference_target, prediction_col, transform,
                          inverse) {
  # With twelve predictors, the horizontal coordinate is the fitted linear
  # predictor. An observation on the identity line equals its PGLS fit.
  predictor_terms <- delete.response(terms(fit$formula))
  predictors <- all.vars(predictor_terms)
  fossil_imputed <- predict(fit$impute_model_traits, climate_traits[, predictors])
  fossil_design <- model.matrix(predictor_terms, fossil_imputed)
  fossil_design <- fossil_design[, colnames(fit$X), drop = FALSE]
  fit_rows <- match(rownames(fit$X), rownames(fit$dat_fit))
  observed <- fit$dat_fit[[all.vars(fit$formula)[1]]][fit_rows]

  calibration <- data.frame(
    id = rownames(fit$X), x = inverse(drop(fit$X %*% fit$beta)),
    y = inverse(observed), record = "extant_species"
  )
  occurrences <- climate_predictions |>
    transmute(id = fossil_name, site, age_ma, x = inverse(drop(fossil_design %*% fit$beta)),
              y = .data[[prediction_col]], record = "fossil_occurrence")

  mean_response <- function(values) mean_on_scale(values, transform, inverse)
  sites <- site_points(occurrences, mean_response, mean_response) |>
    left_join(
      climate_intervals |>
        filter(.data$target == reference_target) |>
        select(site, age_ma, n_species, estimate_response,
               lower_95 = calibrated_lower_response, upper_95 = calibrated_upper_response),
      by = c("site", "age_ma"), relationship = "one-to-one"
    )
  stopifnot(
    all(sites$n_occurrences == sites$n_species),
    isTRUE(all.equal(sites$y, sites$estimate_response, tolerance = 1e-11))
  )
  sites <- sites |> select(-n_species, -estimate_response) |> mutate(record = "fossil_site")

  fit_range <- range(calibration$x, occurrences$x)
  line <- data.frame(x = fit_range, y = fit_range, record = "pgls_fit")
  bind_rows(calibration, occurrences, sites, line) |> mutate(target = target)
}

mat <- climate_panel(climate_fit$mat, "MAT", "mat", "mat_pip_site", identity, identity)
map <- climate_panel(climate_fit$log_map, "MAP", "log_map", "map_pip_site", log, exp)

# The LMA model uses one measured predictor, so its original predictor can be
# shown directly. Site LMA is the geometric mean of occurrence predictions.
lma_calibration <- lma_fit$dat_fit |>
  transmute(id = rownames(lma_fit$dat_fit), x = log10_pw2a_ratio,
            y = 10^log10_lma, record = "extant_species")
lma_occurrences <- lma_predictions |>
  left_join(traits |> select(fossil_name, site, age_ma, pw2.a.ratio),
            by = c("fossil_name", "site", "age_ma"), relationship = "one-to-one") |>
  transmute(id = fossil_name, site, age_ma, x = log10(pw2.a.ratio),
            y = lma_pip, record = "fossil_occurrence")
eligible_lma <- traits$fossil_name[is.finite(log10(traits$pw2.a.ratio))]
stopifnot(setequal(eligible_lma, lma_occurrences$id))
lma_sites <- site_points(lma_occurrences, mean, function(values) {
  mean_on_scale(values, log10, function(value) 10^value)
}) |>
  mutate(record = "fossil_site")
lma_match <- match(lma_sites$site, lma_reference$site)
stopifnot(
  nrow(lma_sites) == nrow(lma_reference),
  all(lma_sites$age_ma == lma_reference$age_ma[lma_match]),
  all(round(lma_sites$y, 1) == lma_reference$lma_pip[lma_match])
)
lma_x <- seq(min(lma_calibration$x, lma_occurrences$x),
             max(lma_calibration$x, lma_occurrences$x), length.out = 80)
lma_beta <- drop(lma_fit$beta)
lma_line <- data.frame(x = lma_x, y = 10^(lma_beta[1] + lma_beta[2] * lma_x),
                       record = "pgls_fit")
lma <- bind_rows(lma_calibration, lma_occurrences, lma_sites, lma_line) |>
  mutate(target = "LMA")

figure_data <- bind_rows(mat, map, lma) |>
  select(target, record, id, site, age_ma, x, y, n_occurrences, lower_95, upper_95)
stopifnot(all(is.finite(figure_data$x)), all(is.finite(figure_data$y)))
write.csv(figure_data, file.path(figures, "Figure_2_data.csv"), row.names = FALSE)

ages <- sort(fossil_site_ages(), decreasing = TRUE)
site_order <- names(ages)
fossil_sites <- data.frame(site = site_order, age_ma = unname(ages))
site_comparison <- read.csv("tables/fossil_site_comparison_formal_only.csv")
published_dilp <- read.csv("data/published_dilp_fossil.csv")

climate_sites <- figure_data |>
  filter(record == "fossil_site", target %in% c("MAT", "MAP")) |>
  select(site, age_ma, target, n_occurrences, pgls = x, phylodilp = y, lower_95, upper_95) |>
  left_join(
    site_comparison |> select(site, age_ma, target = variable, site_lm = lm_site, pip_site),
    by = c("site", "age_ma", "target"), relationship = "one-to-one"
  ) |>
  left_join(
    published_dilp |>
      rename(MAT = published_dilp_mat_c, MAP = published_dilp_map_cm) |>
      pivot_longer(c(MAT, MAP), names_to = "target", values_to = "published_dilp"),
    by = c("site", "target"), relationship = "one-to-one"
  )
stopifnot(
  nrow(climate_sites) == 2 * nrow(fossil_sites),
  setequal(climate_sites$site, site_order),
  !anyNA(climate_sites$site_lm),
  all(round(climate_sites$phylodilp, ifelse(climate_sites$target == "MAT", 1, 0)) ==
        climate_sites$pip_site),
  identical(unique(climate_sites$site[is.na(climate_sites$published_dilp)]), "Bonanza")
)

climate_table <- climate_sites |>
  select(-pgls, -pip_site) |>
  pivot_wider(names_from = target,
              values_from = c(n_occurrences, phylodilp, lower_95, upper_95,
                              site_lm, published_dilp),
              names_glue = "{tolower(target)}_{.value}") |>
  arrange(desc(age_ma))
stopifnot(identical(climate_table$site, site_order),
          identical(climate_table$age_ma, unname(ages)),
          all(climate_table$mat_n_occurrences == climate_table$map_n_occurrences))
write.csv(climate_table, file.path(tables, "Table_6.csv"), row.names = FALSE)

lma_table <- fossil_sites |>
  left_join(lma_sites |> select(site, age_ma, n_occurrences, lma_phylodilp_g_m2 = y),
            by = c("site", "age_ma"), relationship = "one-to-one") |>
  left_join(lma_reference |> select(site, age_ma, lma_lm_g_m2 = lma_lm),
            by = c("site", "age_ma"), relationship = "one-to-one") |>
  mutate(n_occurrences = coalesce(n_occurrences, 0L))
stopifnot(identical(lma_table$site[lma_table$n_occurrences == 0], "Cerrejon"),
          identical(is.na(lma_table$lma_phylodilp_g_m2), is.na(lma_table$lma_lm_g_m2)))
write.csv(lma_table, file.path(tables, "Table_7.csv"), row.names = FALSE)

# Separate the difference between estimators from the change induced by the
# prediction-time correction within the same fitted PGLS model.
climate_shifts <- climate_sites |>
  select(target, site, age_ma, pgls, phylodilp, site_lm, published_dilp) |>
  mutate(
    correction = phylodilp - pgls,
    phylodilp_minus_site_lm = phylodilp - site_lm,
    phylodilp_minus_published_dilp = phylodilp - published_dilp
  ) |>
  arrange(target, desc(age_ma))
write.csv(climate_shifts, file.path("tables", "fossil_climate_shifts.csv"), row.names = FALSE)

site_colours <- setNames(c("#00729B", "#A2470D", "#764695", "#16735F", "#AE3859",
                           "#716021", "#3C5B9C", "#AB6000", "#4B3C94", "#436B10"),
                         site_order)
site_shapes <- setNames(rep(c(21, 24, 22, 23, 25), 2), site_order)

panel <- function(data, title, x_label, y_label) {
  ggplot(mapping = aes(x, y)) +
    geom_point(data = filter(data, record == "extant_species"),
               colour = "#656B72", alpha = 0.21, size = 0.88, shape = 16) +
    geom_line(data = filter(data, record == "pgls_fit"),
              colour = "#424B52", linewidth = 0.38, linetype = "longdash") +
    geom_point(data = filter(data, record == "fossil_occurrence"),
               aes(colour = site, shape = site), alpha = 0.62, fill = NA,
               size = 1.65, stroke = 0.46, show.legend = FALSE) +
    geom_point(data = filter(data, record == "fossil_site"),
               aes(fill = site, shape = site), colour = "#172326", size = 3.35,
               stroke = 0.52, show.legend = TRUE) +
    scale_colour_manual(values = site_colours, guide = "none") +
    scale_fill_manual(values = site_colours, limits = site_order, name = NULL,
                      drop = FALSE) +
    scale_shape_manual(values = site_shapes, limits = site_order, name = NULL,
                       drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE),
           shape = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", margin = margin(b = 9)),
      axis.text = element_text(colour = "#323D44"),
      axis.title = element_text(colour = "#172326"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E7EBED", linewidth = 0.35),
      plot.margin = margin(8, 12, 8, 9),
      legend.key.width = grid::unit(18, "pt"),
      legend.key.height = grid::unit(14, "pt"),
      legend.text = element_text(size = 9)
    )
}

# For climate panels, identical per-axis scales make the 1:1 reference
# comparable on the horizontal and vertical axes.
mat_limits <- range(mat$x, mat$y)
map_limits <- range(map$x, map$y)
map_ticks <- c(25, 50, 100, 200, 400)
mat_plot <- panel(mat, "Mean annual temperature",
                 "PGLS estimate (°C)", "Observed or phyloDiLP (°C)") +
  scale_x_continuous(limits = mat_limits) +
  scale_y_continuous(limits = mat_limits) +
  coord_fixed()
map_plot <- panel(map, "Mean annual precipitation",
                 "PGLS estimate (cm)", "Observed or phyloDiLP (cm)") +
  scale_x_log10(breaks = map_ticks, limits = map_limits) +
  scale_y_log10(breaks = map_ticks, limits = map_limits) +
  coord_fixed()
lma_plot <- panel(lma, "Leaf mass per area",
                 expression(log[10]("petiole width"^2 / "blade area")),
                 expression("Observed or phyloDiLP (g " * m^{-2} * ")")) +
  scale_y_log10(breaks = c(25, 50, 100, 200, 400))

figure <- (mat_plot | map_plot | lma_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(8, 17, 8, 17))) &
  theme(legend.position = "bottom")

ggsave(file.path(figures, "Figure_2.pdf"), figure, device = cairo_pdf,
       width = 11.1, height = 4.4, bg = "white")
ggsave(file.path(figures, "Figure_2.png"), figure, width = 11.1, height = 4.4,
       dpi = 300, bg = "white")
