# =============================================================================
# KNMI Earthquake API Functions
#
# Dit bestand bevat functies voor het ophalen van aardbevingsdata
# van de KNMI FDSN API.
#
# Author: Jesse Frederik
# =============================================================================

library(httr2)
library(purrr)
library(lubridate)
library(dplyr)
library(tibble)

# KNMI FDSN Event API endpoint
knmi_event_url <- "https://rdsa.knmi.nl/fdsnws/event/1/query"

parse_knmi_features <- function(x) {
  # x is the parsed JSON (list). Features are GeoJSON-like.
  feats <- x$features
  if (is.null(feats) || length(feats) == 0) {
    return(tibble())
  }

  map_dfr(feats, function(f) {
    p <- f$properties
    event_id <- f$id %||% NA_character_

    tibble(
      event_id = event_id,
      time_utc = ymd_hms(p$time, tz = "UTC", quiet = TRUE),
      lat = as.numeric(p$lat),
      lon = as.numeric(p$lon),
      depth_km = as.numeric(p$depth),
      mag = as.numeric(p$mag),
      mag_type = as.character(p$mag_type),
      location = as.character(p$location),
      event_type = as.character(p$event_type),
      # Construct a stable “link” back to the service for that event_id:
      link = paste0(
        "https://rdsa.knmi.nl/fdsnws/event/1/query?eventid=",
        event_id,
        "&format=json"
      )
    )
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

fetch_knmi_events <- function(
  eventtype = "induced or triggered event",
  starttime = NULL,
  endtime = NULL,
  limit = 20000
) {
  out <- list()
  offset <- 1

  repeat {
    req <- request(knmi_event_url) |>
      req_url_query(
        format = "json",
        eventtype = eventtype,
        starttime = starttime,
        endtime = endtime,
        limit = limit,
        offset = offset
      ) |>
      req_retry(max_tries = 5)

    resp <- req_perform(req)
    x <- resp_body_json(resp, simplifyVector = FALSE)

    df <- parse_knmi_features(x)
    n <- nrow(df)

    if (n == 0) {
      break
    }

    out[[length(out) + 1]] <- df

    # next page
    offset <- offset + n
    if (n < limit) break
  }

  bind_rows(out)
}
