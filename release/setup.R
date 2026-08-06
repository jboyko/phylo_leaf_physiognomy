# Installs every package the pipeline needs, at the versions it was developed against.
# Run once from the release root: Rscript setup.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "ape", "phytools", "caper", "geiger", "caret", "dplyr", "tidyr", "tibble",
  "ggplot2", "gridExtra", "mvtnorm", "MASS", "remotes", "ipred", "e1071"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

# dilp is pinned to a specific commit of the fork used for this analysis; the
# CRAN release carries different regression constants and will not reproduce it.
DILP_COMMIT <- "b29be909355f7edb7809bf718f189b7a53b789d2"

remotes::install_github("jboyko/dilp", ref = DILP_COMMIT,
                        upgrade = "never", force = FALSE)

stopifnot(requireNamespace("dilp", quietly = TRUE))
cat("Setup complete. dilp version:",
    as.character(utils::packageVersion("dilp")), "\n")
cat("\nSession info:\n")
print(sessionInfo())
