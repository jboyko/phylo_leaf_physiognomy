# Shared fossil taxonomy rules.
#
# Dana Royer (email, April 2026) specified that order, family, or genus names
# inside quotation marks are informal and must not be used for phylogenetic
# placement. The primary `formal_only` scenario censors each quoted rank. The
# `include_informal` sensitivity scenario removes the quote marks and treats the
# reported name provisionally as a placement hypothesis.

taxonomy_ranks <- c("genus", "family", "order")

is_informal_taxon <- function(x) {
  x <- as.character(x)
  grepl('"', x, fixed = TRUE) |
    grepl("\u201c", x, fixed = TRUE) |
    grepl("\u201d", x, fixed = TRUE) |
    grepl("\u2018", x, fixed = TRUE) |
    grepl("\u2019", x, fixed = TRUE)
}

strip_taxon_quotes <- function(x) {
  x <- as.character(x)
  for (mark in c('"', "\u201c", "\u201d", "\u2018", "\u2019")) {
    x <- gsub(mark, "", x, fixed = TRUE)
  }
  trimws(gsub("[[:space:]]+", " ", x))
}

normalise_missing_taxon <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | tolower(x) == "unknown"] <- "unknown"
  x
}

is_known_taxon <- function(x) {
  x <- normalise_missing_taxon(x)
  x != "unknown"
}

taxonomy_for_scenario <- function(data,
                                  scenario = c("formal_only",
                                               "include_informal")) {
  scenario <- match.arg(scenario)
  out <- data

  for (rank in taxonomy_ranks) {
    reported_col <- paste0(rank, "_reported")
    informal_col <- paste0(rank, "_informal")

    reported <- if (reported_col %in% names(out)) {
      out[[reported_col]]
    } else {
      out[[rank]]
    }
    reported <- normalise_missing_taxon(reported)

    informal <- if (informal_col %in% names(out)) {
      as.logical(out[[informal_col]])
    } else {
      is_informal_taxon(reported)
    }
    informal[is.na(informal)] <- FALSE

    cleaned <- normalise_missing_taxon(strip_taxon_quotes(reported))
    if (scenario == "formal_only") {
      cleaned[informal] <- "unknown"
    }

    out[[reported_col]] <- reported
    out[[informal_col]] <- informal
    out[[rank]] <- cleaned
  }

  out$taxonomy_scenario <- scenario
  out
}

fossil_site_ages <- function() {
  c(
    "Fox Hills" = 66.5,
    "Williston Basin I" = 64.75,
    "Williston Basin II" = 63.5,
    "Williston Basin III" = 59.75,
    "Palacio de los Loros PL1" = 61.7,
    "Palacio de los Loros PL2" = 61.7,
    "Cerrejon" = 58.0,
    "Hubble Bubble" = 55.8,
    "Laguna del Hunco" = 51.9,
    "Republic" = 49.4,
    "Bonanza" = 47.3
  )
}

palacio_analysis_site <- function(morphotype) {
  # Subsite assignments preserved from Table S7 of
  # "Fossil specimen scores for phylo-DiLP (June 2012)_adjusted.xlsx".
  pl2_morphotypes <- c("SA050", "SA054", "SA055", "SA057", "SA058", "SA059")
  ifelse(
    trimws(as.character(morphotype)) %in% pl2_morphotypes,
    "Palacio de los Loros PL2",
    "Palacio de los Loros PL1"
  )
}

read_mixed_utf8_csv <- function(path) {
  # The file is UTF-8 except for a few legacy bytes in free-text comments.
  # Preserve valid UTF-8 taxonomy quotes and discard only invalid bytes.
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  lines <- iconv(lines, from = "UTF-8", to = "UTF-8", sub = "")
  read.csv(
    text = lines,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}
