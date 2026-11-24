args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript verify_reproduction.R /path/to/rebuild /path/to/submission_1",
       call. = FALSE)
}
rebuild <- normalizePath(args[1], mustWork = TRUE)
reference <- normalizePath(args[2], mustWork = TRUE)
produced <- file.path(rebuild, "submission/submission_1")

# Keep plotting order where Figure 2 fit lines share a key.
row_keys <- function(data, columns) {
  keys <- do.call(paste, c(data[columns], sep = "\r"))
  order_within_key <- ave(seq_along(keys), keys, FUN = seq_along)
  paste(keys, order_within_key, sep = "\r")
}

compare_csv <- function(path, keys) {
  actual <- read.csv(file.path(produced, path), stringsAsFactors = FALSE)
  expected <- read.csv(file.path(reference, path), stringsAsFactors = FALSE)
  if (!identical(names(actual), names(expected))) stop(path, " columns differ.")
  actual_key <- row_keys(actual, keys)
  expected_key <- row_keys(expected, keys)
  if (length(actual_key) != length(expected_key) || !setequal(actual_key, expected_key)) {
    stop(path, " row identifiers differ.")
  }
  actual <- actual[match(expected_key, actual_key), , drop = FALSE]

  max_difference <- 0
  for (column in names(expected)) {
    x <- actual[[column]]
    y <- expected[[column]]
    exact <- column %in% keys || !is.numeric(y) ||
      grepl("(^|_)n_(sites|occurrences)$|^age_ma$|^fold$", column)
    if (exact) {
      if (!identical(x, y)) stop(path, " differs in ", column, ".")
    } else {
      if (!identical(is.na(x), is.na(y))) stop(path, " missing values differ in ", column, ".")
      present <- !is.na(y)
      difference <- abs(x[present] - y[present])
      tolerance <- 1e-7 * pmax(1, abs(y[present]))
      if (any(!is.finite(difference) | difference > tolerance)) {
        stop(path, " exceeds numerical tolerance in ", column,
             ". Maximum absolute difference: ", max(difference), ".")
      }
      max_difference <- max(max_difference, difference)
    }
  }
  cat(path, "matched", nrow(expected), "rows.",
      "Maximum absolute numerical difference:", format(max_difference, digits = 4), "\n")
}

identifiers <- list(
  "tables/Table_1.csv" = "model",
  "tables/Table_2.csv" = "model",
  "tables/Table_3.csv" = "model",
  "tables/Table_4.csv" = "model",
  "tables/Table_5.csv" = "model",
  "tables/Table_6.csv" = "site",
  "tables/Table_7.csv" = "site",
  "tables/Table_S1.csv" = c("target", "baseline"),
  "tables/Table_S2.csv" = c("target", "model"),
  "figures/Figure_2_data.csv" = c("target", "record", "id", "site", "age_ma"),
  "figures/Figure_S1_data.csv" = c("site", "fold", "target", "model")
)
for (path in names(identifiers)) compare_csv(path, identifiers[[path]])

figures <- file.path(produced, "figures",
                     c("Figure_1.pdf", "Figure_2.pdf", "Figure_2.png", "Figure_S1.pdf"))
stopifnot(all(file.exists(figures)), all(file.info(figures)$size > 0))
cat("Submitted figure and table source data reproduced. All figure files were generated.\n")
