# ==============================================================================
# setup.R
# Shared bootstrap sourced at the top of every pipeline script, in place of a
# hardcoded setwd() to an absolute path (issue #7).
#
# 1. Locates the repository root by walking up from the current working
#    directory looking for README.md alongside a code/ directory, and sets it
#    as the working directory. This works whether the script is launched as
#    `Rscript code/00_data_cleaning.R` from the repo root, or from an R/
#    RStudio session opened anywhere inside the repository (including
#    inside code/).
# 2. Creates the gitignored output directories (models/, tables/, plots/) so
#    a fresh clone does not fail on its first write (issue #21).
#
# Every script sources this file with:
#   source(if (file.exists("code/setup.R")) "code/setup.R" else "setup.R")
# ==============================================================================

.pip_find_root <- function(start = getwd(), marker = "README.md") {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, marker)) && dir.exists(file.path(d, "code"))) {
      return(d)
    }
    parent <- dirname(d)
    if (parent == d) return(NULL)
    d <- parent
  }
}

.pip_root <- .pip_find_root()
if (is.null(.pip_root)) {
  stop(
    "setup.R: could not locate the phylo_leaf_physiognomy project root ",
    "(looked for README.md alongside a code/ directory) starting from '",
    getwd(), "'. Launch scripts with the working directory set inside the ",
    "repository, e.g. run `Rscript code/00_data_cleaning.R` from the repo root."
  )
}
setwd(.pip_root)

for (.pip_dir in c("models", "tables", "plots")) {
  dir.create(.pip_dir, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# dilp package loader (issue #22)
#
# The pipeline depends on the unpublished fork jboyko/dilp, pinned to a
# specific commit because the CRAN release carries different regression
# constants and will not reproduce these results. Scripts that need dilp
# call pip_require_dilp() instead of library(dilp) or devtools::load_all(),
# so there is exactly one loading path and one hard error message.
# ------------------------------------------------------------------------------

PIP_DILP_REF <- "b29be909355f7edb7809bf718f189b7a53b789d2"

pip_require_dilp <- function() {
  if (!requireNamespace("dilp", quietly = TRUE)) {
    stop(
      "The 'dilp' package is required but not installed.\n",
      "Install the pinned fork used for this analysis:\n\n",
      '  remotes::install_github("jboyko/dilp", ref = "', PIP_DILP_REF, '")\n\n',
      "The CRAN release of dilp carries different regression constants and ",
      "will NOT reproduce these results."
    )
  }
  library(dilp, character.only = FALSE)
  invisible(TRUE)
}

rm(.pip_root, .pip_dir, .pip_find_root)
