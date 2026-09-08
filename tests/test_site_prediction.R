source("code/site_prediction.R")

# Average on the fitted log scale: 100 and 400 cm give a 200 cm site estimate.
stopifnot(
  isTRUE(all.equal(exp(site_prediction(log(c(100, 400)), "log_map")), 200)),
  isTRUE(all.equal(site_prediction(c(10, 20), "mat"), 15)),
  isTRUE(all.equal(site_prediction(c(log(100), NA), "log_map"), log(100))),
  is.na(site_prediction(c(NA_real_, NA_real_), "log_map")),
  isTRUE(all.equal(site_prediction(c(1000, 1000), "log_map"), 1000))
)
cat("Site prediction aggregation tests passed\n")
