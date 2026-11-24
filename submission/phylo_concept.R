# A three-species illustration of the PhyloDiLP prediction.
library(grid)

fossil <- "#913D29"
gray <- "#777777"

label <- function(text, x, y, size = 10, just = "left", colour = "black", rot = 0) {
  grid.text(text, x, y, just = just, rot = rot,
            gp = gpar(fontfamily = "serif", fontsize = size, col = colour))
}

segment <- function(x, y, xend, yend, colour = "black", type = 1, width = 0.85) {
  grid.segments(x, y, xend, yend, gp = gpar(col = colour, lwd = width, lty = type))
}

point <- function(x, y, pch, colour) {
  grid.points(x, y, pch = pch, default.units = "npc",
              size = unit(2, "mm"), gp = gpar(col = colour, fill = "white"))
}

# Fix the fitted line and use lambda = 1. A and B have equal weighted residuals.
fit <- function(tooth_fraction) 25 - 15 * tooth_fraction
extant_traits <- c(A = 0.8, B = 0.2, C = 0.4)
residuals <- c(A = -2.8, B = -2.8, C = 3)
extant_cov <- matrix(c(1, 0.4, 0, 0.4, 1, 0, 0, 0, 1), nrow = 3)
shared <- c(A = 0.7, B = 0.4, C = 0)
weighted_residuals <- solve(extant_cov, residuals)
contribution <- shared * weighted_residuals
pgls_mat <- fit(0.5)
pip_mat <- pgls_mat + sum(contribution)

grDevices::cairo_pdf("submission/submission_1/figures/Figure_1.pdf",
                     width = 7.6, height = 2.75)
grid.newpage()

tooth_x <- function(fraction) 0.08 + 0.28 * fraction
mat_y <- function(mat) 0.23 + (mat - 8) / 18 * 0.62
fossil_x <- tooth_x(0.5)
x_axis <- mat_y(8)

label("(a)", 0.02, 0.93, 11)
segment(tooth_x(0), x_axis, tooth_x(1), x_axis)
segment(tooth_x(0), x_axis, tooth_x(0), mat_y(26))
for (tick in c(0, 1)) {
  segment(tooth_x(tick), x_axis, tooth_x(tick), x_axis - 0.024)
  label(tick, tooth_x(tick), 0.16, 9, just = "centre")
}
for (tick in c(10, 20, 25)) {
  segment(tooth_x(0) - 0.008, mat_y(tick), tooth_x(0), mat_y(tick))
  label(tick, tooth_x(0) - 0.014, mat_y(tick), 9, just = "right")
}
label("Fraction of toothed leaves", tooth_x(0.5), 0.073, 10, just = "centre")
label("MAT (°C)", 0.018, 0.52, 10, just = "centre", rot = 90)

segment(tooth_x(0), mat_y(fit(0)), tooth_x(1), mat_y(fit(1)))
observed_y <- mat_y(fit(extant_traits) + residuals)
extant_x <- tooth_x(extant_traits)
segment(extant_x, observed_y, extant_x, mat_y(fit(extant_traits)), gray, 2)
point(extant_x, observed_y, 16, "black")
label(names(extant_traits), extant_x - 0.012, observed_y, 10, just = "right")

segment(fossil_x, x_axis, fossil_x, mat_y(pip_mat), gray, 3)
segment(fossil_x, mat_y(pip_mat), fossil_x, mat_y(pgls_mat), fossil, width = 1.3)
point(fossil_x, x_axis, 18, fossil)
label("F (0.5)", fossil_x, 0.16, 9, just = "centre", colour = fossil)

predictions <- c(PGLS = pgls_mat, PIP = pip_mat)
prediction_y <- mat_y(predictions)
segment(tooth_x(0), prediction_y, fossil_x, prediction_y, fossil, 2)
point(rep(fossil_x, 2), prediction_y, 21, fossil)
label(names(predictions), fossil_x + 0.015, prediction_y + c(0.034, -0.041),
      10, colour = fossil)
label(sprintf("%.1f", predictions), tooth_x(0) - 0.014, prediction_y, 9,
      just = "right", colour = fossil)

# Time increases to the right. C diverges at the root.
tree_x <- function(time) 0.455 + 0.25 * time
tips <- c(A = 0.76, B = 0.52, C = 0.28)
ab <- mean(tips[c("A", "B")])
f <- ab
fossil_time <- 0.82

label("(b)", 0.42, 0.93, 11)
segment(tree_x(0), tips["C"], tree_x(0), ab)
segment(tree_x(0), tips["C"], tree_x(1), tips["C"])
segment(tree_x(0), ab, tree_x(0.4), ab, gray, width = 2.5)
segment(tree_x(0.4), tips["B"], tree_x(0.4), tips["A"])
segment(tree_x(0.4), ab, tree_x(0.4), tips["A"], gray, width = 2.5)
segment(tree_x(0.4), tips["A"], tree_x(0.7), tips["A"], gray, width = 2.5)
segment(tree_x(0.7), tips["A"], tree_x(1), tips["A"])
segment(tree_x(0.4), tips["B"], tree_x(1), tips["B"])
segment(tree_x(0.7), tips["A"], tree_x(0.7), f, fossil)
segment(tree_x(0.7), f, tree_x(fossil_time), f, fossil)
point(rep(tree_x(1), 3), tips, 16, "black")
point(tree_x(fossil_time), f, 18, fossil)
label("F", tree_x(fossil_time) + 0.013, f, 10, colour = fossil)
label(names(tips), tree_x(1) + 0.015, tips, 10)

time_ticks <- tree_x(c(0, fossil_time, 1))
segment(tree_x(0), 0.165, tree_x(1), 0.165, gray, width = 0.6)
segment(time_ticks, 0.165, time_ticks, 0.19, gray, width = 0.6)
label("Root", time_ticks[1], 0.116, 9, just = "centre")
label("Fossil age", time_ticks[2] - 0.006, 0.116, 9, just = "right")
label("Present", time_ticks[3], 0.116, 9, just = "centre")

label("Shared branch\nlength", 0.775, 0.92, 8.6, just = "centre")
label("Weighted\nresidual", 0.861, 0.92, 8.6, just = "centre")
label("Contribution\n(°C)", 0.945, 0.92, 8.6, just = "centre")
label(sprintf("%.1f", shared), 0.775, tips, 10, just = "centre")
label(sprintf("%.1f", weighted_residuals), 0.861, tips, 10, just = "centre")
label(sprintf("%.1f", contribution), 0.945, tips, 10, just = "centre")
invisible(dev.off())
