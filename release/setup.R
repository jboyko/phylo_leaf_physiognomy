# Installs the required CRAN packages and the pinned dilp commit.
# Run once from the release root: Rscript setup.R

source("code/setup.R")
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "ape", "phytools", "caret", "tidyverse",
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
remotes::install_github("jboyko/dilp", ref = PIP_DILP_REF,
                        upgrade = "never", force = FALSE)

stopifnot(requireNamespace("dilp", quietly = TRUE),
          identical(packageDescription("dilp", fields = "RemoteSha"), PIP_DILP_REF))
cat("Setup complete. dilp version:",
    as.character(utils::packageVersion("dilp")), "\n")
cat("\nSession info:\n")
print(sessionInfo())
