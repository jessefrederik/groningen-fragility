# =============================================================================
# Ground Motion Model (PGV) - Bommer et al. (2022)
#
# Dit bestand bevat functies voor het berekenen van Peak Ground Velocity (PGV)
# op basis van het grondbewegingsmodel van Bommer et al. (2022).
#
# @author Jesse Frederik
# @references Bommer, J.J., et al. (2022). Ground-motion model for induced
#             seismicity in the Groningen gas field.
#  @link https://doi.org/10.1007/s10950-022-10120-w
# =============================================================================

#' Get PGV Model Coefficients
#'
#' @description Returns the coefficients for the Bommer et al. (2022) Ground
#'   Motion Model. These coefficients are used to calculate Peak Ground
#'   Velocity (PGV) based on magnitude, distance, and site conditions.
#'
#' @param comp Character. Component definition for PGV calculation.
#'   One of:
#'   \itemize{
#'     \item \code{"GEO"}: Geometric mean of horizontal components
#'     \item \code{"LRG"}: Largest horizontal component (conservative)
#'     \item \code{"MAX"}: Maximum rotated component (most conservative)
#'   }
#'
#' @return A named list containing:
#'   \itemize{
#'     \item \code{c1} through \code{c8}: Model coefficients
#'     \item \code{tau}: Inter-event variability (between-earthquake)
#'     \item \code{phiS2S}: Site-to-site variability
#'     \item \code{phiSS}: Single-station within-event variability
#'   }
#'
#' @examples
#' coef <- get_pgv_coefficients("LRG")
#' coef$c1
#' coef$tau
#'
#' @export
get_pgv_coefficients <- function(comp) {
  coefficients <- list(
    GEO = list(
      c1 = -3.9045,
      c2 = 2.3004,
      c3 = -2.6496,
      c4 = -1.0908,
      c5 = -2.0089,
      c6 = -3.3276,
      c7 = 1.1513,
      c8 = -0.2977,
      tau = 0.2488,
      phiS2S = 0.242,
      phiSS = 0.416
    ),
    LRG = list(
      c1 = -3.3996,
      c2 = 2.3258,
      c3 = -2.8522,
      c4 = -1.0151,
      c5 = -2.1002,
      c6 = -3.4407,
      c7 = 1.1513,
      c8 = -0.3295,
      tau = 0.2448,
      phiS2S = 0.2406,
      phiSS = 0.4569
    ),
    MAX = list(
      c1 = -3.2738,
      c2 = 2.3343,
      c3 = -2.8857,
      c4 = -1.006,
      c5 = -2.1016,
      c6 = -3.394,
      c7 = 1.1513,
      c8 = -0.3354,
      tau = 0.247,
      phiS2S = 0.2442,
      phiSS = 0.453
    )
  )

  if (!comp %in% names(coefficients)) {
    stop("ERROR: argument 'comp' must be one of 'GEO', 'LRG', 'MAX'")
  }

  return(coefficients[[comp]])
}

#' Convert Epicentral to Hypocentral Distance
#'
#' @description Calculates the hypocentral (slant) distance from the epicentral
#'   (surface) distance and earthquake depth using the Pythagorean theorem.
#'
#' @param epi_dist_km Numeric. Epicentral distance in kilometers (horizontal
#'   distance from earthquake epicenter to site).
#' @param depth_km Numeric. Earthquake focal depth in kilometers.
#'
#' @return Numeric. Hypocentral distance in kilometers.
#'
#' @examples
#' # 10 km from epicenter, 3 km deep
#' epi_to_hypo_dist(10, 3)  # Returns ~10.44 km
#'
#' @export
epi_to_hypo_dist <- function(epi_dist_km, depth_km) {
  sqrt(epi_dist_km^2 + depth_km^2)
}

#' Convert Hypocentral to Epicentral Distance
#'
#' @description Calculates the epicentral (surface) distance from the
#'   hypocentral distance and earthquake depth.
#'
#' @param hypo_dist_km Numeric. Hypocentral distance in kilometers.
#' @param depth_km Numeric. Earthquake focal depth in kilometers.
#'
#' @return Numeric. Epicentral distance in kilometers.
#'
#' @examples
#' hypo_to_epi_dist(10.44, 3)  # Returns ~10 km
#'
#' @export
hypo_to_epi_dist <- function(hypo_dist_km, depth_km) {
  sqrt(hypo_dist_km^2 - depth_km^2)
}

#' Calculate Total Sigma for PGV Model
#'
#' @description Computes the total standard deviation (sigma) of the PGV model
#'   by combining inter-event, site-to-site, and single-station variability
#'   components in quadrature.
#'
#' @param comp Character. Component definition: \code{"GEO"}, \code{"LRG"},
#'   or \code{"MAX"}.
#'
#' @return Numeric. Total sigma value (natural log units).
#'
#' @details The total sigma is calculated as:
#'   \deqn{\sigma_{total} = \sqrt{\tau^2 + \phi_{S2S}^2 + \phi_{SS}^2}}
#'
#' @examples
#' pgv_sigma("LRG")  # Returns ~0.56
#'
#' @export
pgv_sigma <- function(comp) {
  coef <- get_pgv_coefficients(comp)
  sqrt(coef$tau^2 + coef$phiS2S^2 + coef$phiSS^2)
}

#' Calculate Peak Ground Velocity (PGV)
#'
#' @description Computes the Peak Ground Velocity using the Bommer et al. (2022)
#'   Ground Motion Model developed specifically for induced seismicity in the
#'   Groningen gas field.
#'
#' @param M Numeric. Local magnitude (ML) of the earthquake.
#' @param Rhyp Numeric. Hypocentral distance in kilometers.
#' @param comp Character. Component definition: \code{"GEO"}, \code{"LRG"},
#'   or \code{"MAX"}.
#' @param VS30 Numeric. Time-averaged shear-wave velocity in the upper 30m
#'   of soil, in m/s. Typical range: 150-300 m/s for Groningen.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item \code{PGV}: Peak Ground Velocity in cm/s
#'     \item \code{lnPGV}: Natural logarithm of PGV
#'     \item \code{tau}: Inter-event variability component
#'     \item \code{phiS2S}: Site-to-site variability component
#'     \item \code{phiSS}: Single-station within-event variability
#'     \item \code{sigma}: Total standard deviation
#'   }
#'
#' @details The model uses a piecewise distance function with hinges at 7 km
#'   and 12 km. The effective distance includes a magnitude-dependent
#'   pseudo-depth term to account for finite rupture effects.
#'
#'   Formula:
#'   \deqn{ln(PGV) = c_1 + c_2 M + c_3 ln(min(R,7)) + c_4 ln(max(min(R,12),7)/7)
#'                 + c_5 ln(max(R,12)/12) + c_8 ln(VS30/200)}
#'
#' @examples
#' # M3.6 earthquake, 5km away, VS30=200 m/s
#' result <- pgv_model(M = 3.6, Rhyp = 5, comp = "LRG", VS30 = 200)
#' result$PGV * 10  # Convert to mm/s
#'
#' @references
#' Bommer, J.J., et al. (2022). Development of Version 7 GMPEs for Response
#' Spectral Accelerations and Significant Durations from Induced Earthquakes
#' in the Groningen Field.
#'
#' @export
pgv_model <- function(M, Rhyp, comp, VS30) {
  coef <- get_pgv_coefficients(comp)

  # Calculate effective distance with magnitude-dependent correction
  hM <- exp(coef$c6 + coef$c7 * M)
  R <- sqrt(Rhyp^2 + hM^2)

  # Calculate ln(PGV) with piecewise distance function
  lnPGV <- coef$c1 +
    coef$c2 * M +
    coef$c3 * log(pmin(R, 7)) +
    coef$c4 * log(pmax(pmin(R, 12), 7) / 7) +
    coef$c5 * log(pmax(R, 12) / 12) +
    coef$c8 * log(VS30 / 200)

  PGV <- exp(lnPGV)
  sigma <- sqrt(coef$tau^2 + coef$phiS2S^2 + coef$phiSS^2)

  return(list(
    PGV = PGV,
    lnPGV = lnPGV,
    tau = coef$tau,
    phiS2S = coef$phiS2S,
    phiSS = coef$phiSS,
    sigma = sigma
  ))
}

#' Calculate PGV for a Data Frame
#'
#' @description Vectorized wrapper to calculate PGV for multiple buildings
#'   or observations stored in a data frame.
#'
#' @param data A data frame containing columns:
#'   \itemize{
#'     \item \code{mag}: Earthquake magnitude
#'     \item \code{distance_km}: Epicentral distance in km
#'     \item \code{depth_km}: Earthquake depth in km
#'     \item \code{VS30}: Shear-wave velocity in m/s
#'   }
#' @param comp Character. Component definition: \code{"GEO"}, \code{"LRG"},
#'   or \code{"MAX"}. Default is \code{"LRG"}.
#'
#' @return The input data frame with additional columns:
#'   \itemize{
#'     \item \code{hypo_dist}: Hypocentral distance in km
#'     \item \code{ln_pgv}: Natural log of PGV
#'     \item \code{pgv}: PGV in cm/s
#'     \item \code{sigma}: Total model uncertainty
#'   }
#'
#' @examples
#' df <- data.frame(
#'   mag = c(3.0, 3.5),
#'   distance_km = c(5, 10),
#'   depth_km = c(3, 3),
#'   VS30 = c(200, 180)
#' )
#' result <- calculate_pgv_for_dataframe(df, comp = "LRG")
#'
#' @export
calculate_pgv_for_dataframe <- function(data, comp = "LRG") {
  # Calculate hypocentral distance
  hypo_dist <- epi_to_hypo_dist(data$distance_km, data$depth_km)

  # Call pgv_model (already vectorized)
  result <- pgv_model(data$mag, hypo_dist, comp, data$VS30)

  # Add results to dataframe
  data$hypo_dist <- hypo_dist
  data$ln_pgv <- result$lnPGV
  data$pgv <- result$PGV
  data$sigma <- result$sigma

  data
}
