# Install the ordinary R dependencies used by the release scripts.
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "ape", "phytools", "caret", "tidyverse", "magrittr", "stringr",
  "ggplot2", "gridExtra", "mvtnorm", "MASS", "ipred", "e1071", "patchwork"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
unavailable <- cran_packages[!vapply(cran_packages, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(unavailable)) stop("Missing packages: ", paste(unavailable, collapse = ", "))
if (!capabilities("cairo")) stop("The paper figures require R with Cairo support.")

versions <- vapply(cran_packages, function(package) {
  as.character(utils::packageVersion(package))
}, character(1))
print(versions)
print(sessionInfo())
