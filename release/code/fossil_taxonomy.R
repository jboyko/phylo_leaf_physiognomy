# Shared fossil taxonomy rules.
#
# Dana Royer (email, April 2026) specified that quoted order, family, or genus
# names are informal. A trailing question mark also marks a rank as uncertain.
# The primary `formal_only` scenario skips marked ranks. The `include_informal`
# sensitivity scenario removes the marks and tries the names provisionally.

taxonomy_ranks <- c("genus", "family", "order")
taxon_question_suffix <- "[[:space:]]*[?]+$"

is_informal_taxon <- function(x) {
  x <- as.character(x)
  grepl('"', x, fixed = TRUE) |
    grepl("\u201c", x, fixed = TRUE) |
    grepl("\u201d", x, fixed = TRUE) |
    grepl("\u2018", x, fixed = TRUE) |
    grepl("\u2019", x, fixed = TRUE) |
    has_trailing_taxon_question(x)
}

strip_taxon_quotes <- function(x) {
  x <- as.character(x)
  for (mark in c('"', "\u201c", "\u201d", "\u2018", "\u2019")) {
    x <- gsub(mark, "", x, fixed = TRUE)
  }
  trimws(gsub("[[:space:]]+", " ", x))
}

has_trailing_taxon_question <- function(x) {
  grepl(taxon_question_suffix, strip_taxon_quotes(x))
}

strip_taxon_markers <- function(x) {
  sub(taxon_question_suffix, "", strip_taxon_quotes(x))
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

    informal <- is_informal_taxon(reported)
    if (informal_col %in% names(out)) {
      informal <- informal | (as.logical(out[[informal_col]]) %in% TRUE)
    }

    cleaned <- normalise_missing_taxon(strip_taxon_markers(reported))
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
  ages <- read.csv("data/fossil_site_ages.csv", stringsAsFactors = FALSE)
  stopifnot(nrow(ages) > 0L, !anyNA(ages$site), all(nzchar(ages$site)),
            !anyDuplicated(ages$site), all(is.finite(ages$age_ma)),
            all(ages$age_ma >= 0))
  setNames(ages$age_ma, ages$site)
}

normalise_fossil_site <- function(site) {
  site <- trimws(as.character(site))
  site[site %in% c("Palacio de los Loros PL1", "Palacio de los Loros PL2")] <-
    "Palacio de los Loros"
  site
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
