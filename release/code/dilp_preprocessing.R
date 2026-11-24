# SPDX-License-Identifier: GPL-3.0-or-later
# DiLP leaf processing by Matthew J. Butrim and DiLP contributors.
# Source, changes and license are documented in licenses/DiLP.md.
#
# The paper uses the processed traits. The package's diagnostics and climate
# regression are not part of this preprocessing entry point.
pip_preprocess_dilp <- local({
  `%>%` <- magrittr::`%>%`

  colnameClean <- function(data) {
    names <- colnames(data) %>%
      stringr::str_trim() %>%
      stringr::str_to_lower() %>%
      stringr::str_replace_all("[.]", " ") %>%
      stringr::str_replace_all("[ ]", "_") %>%
      stringr::str_replace_all("__", "_")
    names <- colnameSynonym(names, "petiole_metric", c("pwa", "pw^2/a"))
    names <- colnameSynonym(names, "leaf_area", c("la"))
    names <- colnameSynonym(names, "petiole_width", c("pw"))
    names <- colnameSynonym(names, "internal_raw_blade_area", c("raw_internal_blade_area"))
    names <- colnameSynonym(names, "internal_raw_blade_perimeter", c("raw_internal_blade_perimeter"))
    names <- colnameSynonym(names, "no_of_primary_teeth", c("no_primary_teeth", "number_primary_teeth", "number_of_primary_teeth", "#primary_teeth", "#_of_primary_teeth"))
    names <- colnameSynonym(names, "no_of_subsidiary_teeth", c(
      "no_subsidiary_teeth", "number_subsidiary_teeth", "number_of_subsidiary_teeth", "#subsidiary_teeth", "#_of_subsidiary_teeth",
      "no_secondary_teeth", "number_secondary_teeth", "number_of_secondary_teeth", "#secondary_teeth", "#_of_secondary_teeth",
      "no_of_secondary_teeth"
    ))
    names <- colnameSynonym(names, "margin", c("Margin (0 - Toothed; 1 - Untoothed)"))
    return(names)
  }

  colnameSynonym <- function(data, name, synonyms) {
    for (i in 1:length(synonyms)) {
      index <- which(data == tolower(synonyms[i]))
      data[index] <- name
    }
    return(data)
  }

  dilp_processing <- function(specimen_data) {
    specimen_data <- dplyr::as_tibble(specimen_data)
    colnames(specimen_data) <- colnameClean(specimen_data)

    required_columns <- c(
      "site", "specimen_number", "morphotype", "margin", "feret", "blade_area",
      "raw_blade_perimeter", "internal_raw_blade_perimeter", "length_of_cut_perimeter",
      "no_of_primary_teeth", "no_of_subsidiary_teeth"
    )

    recommended_columns <- c(
      "petiole_width", "petiole_area", "blade_perimeter",
      "minimum_feret", "raw_blade_area", "internal_raw_blade_area"
    )
    missing_columns <- required_columns[!required_columns %in% colnames(specimen_data)]
    if (length(missing_columns) > 0) {
      stop(paste("specimen_data is missing required columns:", stringr::str_flatten(missing_columns, collapse = ", ")))
    } else {
      other_columns <- recommended_columns[!recommended_columns %in% colnames(specimen_data)]
      if (length(other_columns) > 0) {
        warning(paste("specimen_data is missing recommended columns and has filled them with NAs:", stringr::str_flatten(other_columns, collapse = ", ")), call. = FALSE)
        for (column in other_columns) {
          specimen_data[[column]] <- NA
        }
      }

      # Identify 0.5 margin species
      mixed_margins <- specimen_data %>%
        dplyr::group_by(.data$site, .data$morphotype) %>%
        dplyr::filter(1 %in% .data$margin & 0 %in% .data$margin) %>%
        dplyr::filter(.data$margin == 1) %>%
        dplyr::select("site", "morphotype", "specimen_number")


      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$length_of_cut_perimeter <- 0
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$no_of_primary_teeth <- 0
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$no_of_subsidiary_teeth <- 0
        }
      }


      ## For specimens where blade area was measured, but there is an NA for petiole area, convert NA to 0 to permit the subsequent addition
      temp3 <- specimen_data$petiole_area
      specimen_data$petiole_area <- ifelse(specimen_data$blade_area > 0 & is.na(specimen_data$petiole_area), 0, specimen_data$petiole_area)
      # Sum petiole area and blade area
      specimen_data$leaf_area <- specimen_data$blade_area + specimen_data$petiole_area # if petiole area was measured but there is a NA for blade area, this addition will appropriately return a NA for leaf area

      # Feret Diameter Ratio
      specimen_data$feret_diameter <- 2 * sqrt(specimen_data$leaf_area / pi)
      specimen_data$fdr <- specimen_data$feret_diameter / specimen_data$feret

      # change NA's to zeros for length.of.cut.perimeter so that derived variables can be calculated
      temp1 <- specimen_data$length_of_cut_perimeter
      specimen_data$length_of_cut_perimeter <- ifelse(is.na(specimen_data$length_of_cut_perimeter), 0, specimen_data$length_of_cut_perimeter)

      # Corrected perimeters
      specimen_data$raw_blade_perimeter <- ifelse(is.na(specimen_data$raw_blade_perimeter) & !is.na(specimen_data$blade_perimeter) & specimen_data$margin == 0 & specimen_data$length_of_cut_perimeter == 0,
        specimen_data$blade_perimeter, specimen_data$raw_blade_perimeter
      )
      specimen_data$raw_blade_perimeter_corrected <- specimen_data$raw_blade_perimeter - specimen_data$length_of_cut_perimeter
      specimen_data$internal_raw_blade_perimeter_corrected <- specimen_data$internal_raw_blade_perimeter - specimen_data$length_of_cut_perimeter

      # Corrected raw blade area if there is no cut perimeter on a toothed leaf.
      specimen_data$raw_blade_area <- ifelse(is.na(specimen_data$raw_blade_area) & specimen_data$margin == 0 & specimen_data$length_of_cut_perimeter == 0,
        specimen_data$blade_area, specimen_data$raw_blade_area
      )

      ## Toothed variables
      # Total tooth count
      temp2 <- specimen_data$no_of_subsidiary_teeth
      specimen_data$no_of_subsidiary_teeth[is.na(specimen_data$no_of_subsidiary_teeth)] <- 0 # if it is left as NA then total tooth count will also return NA
      specimen_data$total_tooth_count <- ifelse(specimen_data$margin == 1, NA, specimen_data$no_of_primary_teeth + specimen_data$no_of_subsidiary_teeth)
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$total_tooth_count <- 0
        }
      }

      # Total tooth count : internal perimeter
      specimen_data$tc_ip <- specimen_data$total_tooth_count / specimen_data$internal_raw_blade_perimeter_corrected
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_ip <- 0
        }
      }

      # Perimeter ratio
      specimen_data$perimeter_ratio <- specimen_data$raw_blade_perimeter_corrected / specimen_data$internal_raw_blade_perimeter_corrected
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$perimeter_ratio <- 1
        }
      }

      #### Apply required natural logs to appropriate variables
      specimen_data$ln_leaf_area <- log(100 * specimen_data$leaf_area) # Leaf area is expressed in mm2 in DiLP MLR and SLR models, so here also converting from cm2 to mm2
      specimen_data$ln_pr <- log(specimen_data$perimeter_ratio)
      specimen_data$ln_tc_ip <- log(specimen_data$tc_ip)
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ln_tc_ip <- NA
        }
      }

      # Petiole metric for reconstructing leaf mass per area
      specimen_data$petiole_metric <- (specimen_data$petiole_width^2) / specimen_data$leaf_area
      # Aspect ratio
      specimen_data$aspect_ratio <- specimen_data$feret / specimen_data$minimum_feret
      # Shape factor
      specimen_data$shape_factor <- 4 * pi * (specimen_data$blade_area / (specimen_data$blade_perimeter^2))
      # Compactness
      specimen_data$compactness <- (specimen_data$blade_perimeter^2) / specimen_data$blade_area
      # Tooth area
      specimen_data$tooth_area <- specimen_data$raw_blade_area - specimen_data$internal_raw_blade_area
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tooth_area <- 0
        }
      }

      # Tooth area : perimeter
      specimen_data$ta_p <- specimen_data$tooth_area / specimen_data$raw_blade_perimeter_corrected
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_p <- 0
        }
      }

      # Tooth area : internal perimeter
      specimen_data$ta_ip <- specimen_data$tooth_area / specimen_data$internal_raw_blade_perimeter_corrected
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_ip <- 0
        }
      }

      # Tooth area : blade area
      specimen_data$ta_ba <- specimen_data$tooth_area / specimen_data$raw_blade_area
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$ta_ba <- 0
        }
      }

      # Average primary tooth area
      specimen_data$avg_ta <- specimen_data$tooth_area / specimen_data$no_of_primary_teeth # double check
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$avg_ta <- 0
        }
      }

      # Tooth count : blade area
      specimen_data$tc_ba <- specimen_data$total_tooth_count / specimen_data$raw_blade_area
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_ba <- 0
        }
      }

      # tooth count : perimeter
      specimen_data$tc_p <- specimen_data$total_tooth_count / specimen_data$raw_blade_perimeter_corrected
      if (length(mixed_margins$specimen_number) > 0) {
        for (i in 1:length(mixed_margins$specimen_number)) {
          specimen_data[which(specimen_data$specimen_number == mixed_margins$specimen_number[i]), ]$tc_p <- 0
        }
      }


      # Return length of cut perimeter and no subsidiary teeth back to original values to keep averages intact
      specimen_data$length_of_cut_perimeter <- temp1
      specimen_data$no_of_subsidiary_teeth <- temp2
      specimen_data$petiole_area <- temp3




      return(specimen_data)
    }
  }

  .dilp_morphotype_means <- function(processed_specimen_data) {
    processed_specimen_data %>%
      dplyr::select(-dplyr::any_of("specimen_number")) %>%
      dplyr::summarize(
        .by = c("site", "morphotype"),
        dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE))
      )
  }

  .dilp_summarize_traits <- function(processed_specimen_data) {
    morphotype <- .dilp_morphotype_means(processed_specimen_data)
    margin <- morphotype$margin
    if (length(unique(stats::na.omit(margin))) > 2) {
      morphotype$margin[margin > 0 & margin < 1] <- 0.5
    }
    if (any(margin < 0 | margin > 1, na.rm = TRUE)) {
      warning("Morphotype margin scores outside [0, 1] are present.", call. = FALSE)
    }

    # Average morphotypes equally. The DiLP site equations use percent margin.
    site <- morphotype %>%
      dplyr::group_by(.data$site) %>%
      dplyr::select(-"morphotype") %>%
      dplyr::summarize(
        dplyr::across(dplyr::where(~ !is.character(.)), \(x) mean(x, na.rm = TRUE))
      )
    site$margin <- site$margin * 100
    site$tc_ip <- ifelse(site$margin == 100, 0, site$tc_ip)
    site$perimeter_ratio <- ifelse(site$margin == 100, 1, site$perimeter_ratio)

    list(morphotype = morphotype, site = site)
  }

  function(specimen_data) {
    colnames(specimen_data) <- colnameClean(specimen_data)
    leaf <- dilp_processing(specimen_data)
    traits <- .dilp_summarize_traits(leaf)
    list(processed_leaf_data = leaf, processed_morphotype_data = traits$morphotype,
         processed_site_data = traits$site)
  }
})
