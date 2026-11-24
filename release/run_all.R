# Runs the full pipeline in order from the release root: Rscript run_all.R
# Scripts are also runnable individually, e.g. Rscript code/02_phy_regression.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--paper")) stop("Expected --paper or no arguments.")

scripts <- c(
  # Climate pipeline
  "code/00_data_cleaning.R",
  "code/00c_fossil_data_cleaning.R",
  "code/01_nophy_regression.R",
  "code/02_phy_regression.R",
  "code/03_loso_cv.R",
  "code/03c_dilp_cv.R",
  "code/03d_uncertainty_diagnostics.R",
  "code/04_fossil_predictions.R",
  "code/05_visualizations.R",
  # Degradation and adjustment-field analyses
  "code/07_type2_degradation.R",
  "code/08_type1_degradation.R",
  "code/09_phylo_adj_field.R",
  # LMA pipeline
  "code/00b_data_cleaning.R",
  "code/01b_nophy_regression.R",
  "code/02b_phy_regression.R",
  "code/03b_loso_cv.R",
  "code/04b_lma_fossil_predictions.R",
  "code/05b_lma_visualizations.R",
  "code/05c_paired_site_comparison.R"
)

if (identical(args, "--paper")) {
  Sys.setenv(PIP_CV_SCHEME = "ten_fold")
  Sys.unsetenv("PIP_CV_FOLDS")
  supplementary_diagnostics <- c(
    "code/03c_dilp_cv.R", "code/05_visualizations.R", "code/07_type2_degradation.R",
    "code/08_type1_degradation.R", "code/09_phylo_adj_field.R",
    "code/05b_lma_visualizations.R"
  )
  scripts <- c(setdiff(scripts, supplementary_diagnostics),
               "paper/phylo_concept.R", "paper/validation_displays.R",
               "paper/fossil_calibration_figure.R")
}

for (s in scripts) {
  cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")
  started <- proc.time()["elapsed"]
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", s))
  if (status != 0) stop("Failed: ", s)
  cat(sprintf("-- %s finished in %.1f s\n", s, proc.time()["elapsed"] - started))
}

cat("\nPipeline complete. Outputs in models/, tables/ and plots/.\n")
