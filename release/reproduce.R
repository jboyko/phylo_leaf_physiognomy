args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript reproduce.R /path/to/Data_S1.zip /path/to/new-output",
       call. = FALSE)
}
archive <- normalizePath(args[1], mustWork = TRUE)
output <- path.expand(args[2])
if (file.exists(output)) {
  stop("Choose a new output directory: ", output, call. = FALSE)
}

inputs <- c(
  "Data_S1_README.md",
  "data/Peppe_2011_calibration_data_leaf_level_clean.csv",
  "data/Peppe_2011_fossil_data_June_2026_leaf_level_clean.csv",
  "data/extra_calibration_data_for_LMA.csv",
  "data/best_wcvp.tre_dated",
  "data/fossil_site_ages.csv",
  "data/published_dilp_fossil.csv"
)
archive_members <- utils::unzip(archive, list = TRUE)$Name
if (!all(inputs %in% archive_members) || anyDuplicated(archive_members)) {
  stop("The archive must contain one copy of each Data S1 input.", call. = FALSE)
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
release <- dirname(normalizePath(script, mustWork = TRUE))
stopifnot(dir.create(output, recursive = TRUE))
source_files <- c("README.md", "code", "licenses", "paper", "tests",
                  "setup.R", "run_all.R")
stopifnot(all(file.copy(file.path(release, source_files), output, recursive = TRUE)))
utils::unzip(archive, files = inputs, exdir = output)
stopifnot(all(file.exists(file.path(output, inputs))))
setwd(output)
dir.create("submission/submission_1/figures", recursive = TRUE)
dir.create("submission/submission_1/tables")
writeLines(capture.output(sessionInfo()), "R-session.txt")
cat("Building the paper in", getwd(), "\n")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", "run_all.R", "--paper"))
if (status != 0L) stop("Paper rebuild failed.", call. = FALSE)
