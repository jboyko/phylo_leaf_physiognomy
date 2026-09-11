# Dana Royer: ridgetop and upper slope form one extant calibration site,
# following the original collector's advice. Valley bottom remains separate.
normalise_calibration_site <- function(site) {
  site <- trimws(as.character(site))
  site[site %in% c("Yasuni-ridgetop", "Yasuni-upper slope")] <-
    "Yasuni-ridgetop and upper slope"
  site
}
