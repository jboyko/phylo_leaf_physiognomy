source("code/site_grouping.R")
d <- read.csv("data/Peppe_2011_calibration_data_leaf_level_clean.csv",
              fileEncoding = "latin1", stringsAsFactors = FALSE)
sites <- normalise_calibration_site(d$site)
stopifnot(length(unique(sites)) == 92L,
          "Yasuni-valley bottom" %in% sites,
          !any(c("Yasuni-ridgetop", "Yasuni-upper slope") %in% sites),
          identical(normalise_calibration_site(sites), sites))
# Both source collections must enter one fold, while unrelated site labels
# retain their identities.
pair <- d$site %in% c("Yasuni-ridgetop", "Yasuni-upper slope")
stopifnot(length(unique(sites[pair])) == 1L,
          identical(sites[!pair], trimws(d$site[!pair])))
source("code/fossil_taxonomy.R")
expected <- c("Fox Hills" = 67, "Williston Basin I" = 65.07,
  "Palacio de los Loros" = 64.08, "Williston Basin II" = 63.8,
  "Williston Basin III" = 60.4, "Cerrejon" = 59,
  "Hubble Bubble" = 55.8, "Laguna del Hunco" = 52,
  "Republic" = 51.18, "Bonanza" = 47.3)
stopifnot(identical(fossil_site_ages()[names(expected)], expected))
cat("92-site calibration grouping and revised fossil ages passed.\n")
