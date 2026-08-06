# Common setup sourced by every analysis script.
# Resolves the project root, creates output directories, and fixes the RNG seed.

if (!file.exists("code/_setup.R")) {
  stop("Run scripts from the release root, e.g. Rscript code/00_data_cleaning.R")
}

SEED <- 42L
set.seed(SEED)

for (d in c("models", "tables", "plots", "dilp_update")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# The pipeline requires the pinned dilp fork installed by setup.R.
if (!requireNamespace("dilp", quietly = TRUE)) {
  stop("Package 'dilp' not found. Run: Rscript setup.R")
}

DILP_COMMIT_REQUIRED <- "b29be909355f7edb7809bf718f189b7a53b789d2"
