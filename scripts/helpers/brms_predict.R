# =============================================================================
# brms_predict.R
#
# Wrapper functions for using brms hurdle-gamma model predictions.
# Maps continuous PGV to ordered factor levels and handles posterior sampling.
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

#' Map continuous PGV values to ordered factor levels
#'
#' @param pgv Numeric vector of PGV values (mm/s)
#' @return Ordered factor with levels matching brms model
map_pgv_to_factor <- function(pgv) {
  pgv_levels <- c(2, 4, 8, 16, 32, 64, 96, 128)

  # Handle single values
  if (length(pgv) == 1) {
    pgv_mapped <- pgv_levels[which.min(abs(pgv_levels - pgv))]
  } else {
    pgv_mapped <- sapply(pgv, function(p) {
      pgv_levels[which.min(abs(pgv_levels - p))]
    })
  }

  factor(pgv_mapped, levels = pgv_levels, ordered = TRUE)
}

#' Map continuous N values to ordered factor levels
#'
#' @param n Numeric vector of N values
#' @return Ordered factor with levels matching brms model
map_n_to_factor <- function(n) {
  n_levels <- c(1, 2, 3, 4, 8)

  # Handle single values
  if (length(n) == 1) {
    n_mapped <- n_levels[which.min(abs(n_levels - n))]
  } else {
    n_mapped <- sapply(n, function(x) {
      n_levels[which.min(abs(n_levels - x))]
    })
  }

  factor(n_mapped, levels = n_levels, ordered = TRUE)
}

#' Map material strength to ordered factor levels
#'
#' @param material Numeric vector of material values
#' @return Ordered factor with levels matching brms model
map_material_to_factor <- function(material) {
  m_levels <- c(0.7, 1.0, 1.3)

  # Handle single values
  if (length(material) == 1) {
    m_mapped <- m_levels[which.min(abs(m_levels - material))]
  } else {
    m_mapped <- sapply(material, function(m) {
      m_levels[which.min(abs(m_levels - m))]
    })
  }

  factor(m_mapped, levels = m_levels, ordered = TRUE)
}

#' Predict damage from brms hurdle-gamma model
#'
#' @param fit_brms brms model object
#' @param pgv PGV values (mm/s)
#' @param n Number of earthquake repetitions (default 1)
#' @param initial_psi Initial damage state
#' @param material Material strength factor
#' @param facade_type "A" or "B"
#' @param soil_profile "A" or "B"
#' @param earthquake_type Earthquake record type (default "ZN")
#' @param ndraws Number of posterior draws to use
#' @param summary If TRUE, return summary statistics; if FALSE, return draws
#' @return Data frame with predictions
predict_damage_brms <- function(
  fit_brms,
  pgv,
  n = 1,
  initial_psi = 0,
  material = 1.0,
  facade_type = "A",
  soil_profile = "A",
  earthquake_type = "ZN",
  ndraws = 100,
  summary = TRUE
) {
  # Get standardization parameters
  std_params <- readRDS(here::here("outputs", "models", "standardization_params.rds"))

  # Prepare new data
  newdata <- data.frame(
    pgv_ord = map_pgv_to_factor(pgv),
    n_ord = map_n_to_factor(n),
    material_ord = map_material_to_factor(material),
    initial_psi_c = initial_psi - std_params$initial_psi_mean,
    FacadeType = factor(facade_type, levels = c("A", "B")),
    SoilProfile = factor(soil_profile, levels = c("A", "B")),
    EarthquakeType = factor(earthquake_type, levels = c("ZN", "ZF", "WN", "WF"))
  )

  # Get posterior predictions
  pred <- posterior_predict(fit_brms, newdata = newdata, ndraws = ndraws)

  if (summary) {
    # Return summary statistics
    data.frame(
      pgv = pgv,
      mean = colMeans(pred),
      median = apply(pred, 2, median),
      sd = apply(pred, 2, sd),
      q05 = apply(pred, 2, quantile, 0.05),
      q25 = apply(pred, 2, quantile, 0.25),
      q75 = apply(pred, 2, quantile, 0.75),
      q95 = apply(pred, 2, quantile, 0.95),
      p_zero = colMeans(pred == 0),
      p_visible = colMeans(pred >= 1)
    )
  } else {
    # Return raw draws (ndraws x n_obs matrix)
    pred
  }
}

#' Simulate damage increment for a single event using brms
#'
#' @param fit_brms brms model object
#' @param pgv PGV value (mm/s) - single value or vector matching n_buildings
#' @param current_psi Current damage state (vector for each building)
#' @param material Material factor (vector)
#' @param facade_type Facade type (vector)
#' @param soil_profile Soil profile (vector)
#' @param ndraws Number of posterior draws per building
#' @return Matrix of damage increments (ndraws x n_buildings)
simulate_damage_increment_brms <- function(
  fit_brms,
  pgv,
  current_psi,
  material = 1.0,
  facade_type = "A",
  soil_profile = "A",
  ndraws = 1
) {
  n_buildings <- length(current_psi)

  # Recycle single values
  if (length(pgv) == 1) pgv <- rep(pgv, n_buildings)
  if (length(material) == 1) material <- rep(material, n_buildings)
  if (length(facade_type) == 1) facade_type <- rep(facade_type, n_buildings)
  if (length(soil_profile) == 1) soil_profile <- rep(soil_profile, n_buildings)

  # Get standardization parameters
  std_params <- readRDS(here::here("outputs", "models", "standardization_params.rds"))

  # Prepare new data
  newdata <- data.frame(
    pgv_ord = map_pgv_to_factor(pgv),
    n_ord = factor(1, levels = c(1, 2, 3, 4, 8), ordered = TRUE),  # Single event
    material_ord = map_material_to_factor(material),
    initial_psi_c = current_psi - std_params$initial_psi_mean,
    FacadeType = factor(facade_type, levels = c("A", "B")),
    SoilProfile = factor(soil_profile, levels = c("A", "B")),
    EarthquakeType = factor("ZN", levels = c("ZN", "ZF", "WN", "WF"))
  )

  # Get posterior predictions - returns matrix (ndraws x n_buildings)
  pred <- posterior_predict(fit_brms, newdata = newdata, ndraws = ndraws)

  # Return as vector if ndraws=1, matrix otherwise
  if (ndraws == 1) {
    as.vector(pred)
  } else {
    pred
  }
}

#' Get P(damage > 0) from brms model
#'
#' @param fit_brms brms model object
#' @param pgv PGV values
#' @param initial_psi Initial damage state
#' @param ndraws Number of draws
#' @return Vector of probabilities
get_p_damage_brms <- function(fit_brms, pgv, initial_psi = 0, ndraws = 500) {
  std_params <- readRDS(here::here("outputs", "models", "standardization_params.rds"))

  newdata <- data.frame(
    pgv_ord = map_pgv_to_factor(pgv),
    n_ord = factor(1, levels = c(1, 2, 3, 4, 8), ordered = TRUE),
    material_ord = factor(1.0, levels = c(0.7, 1.0, 1.3), ordered = TRUE),
    initial_psi_c = initial_psi - std_params$initial_psi_mean,
    FacadeType = factor("A", levels = c("A", "B")),
    SoilProfile = factor("A", levels = c("A", "B")),
    EarthquakeType = factor("ZN", levels = c("ZN", "ZF", "WN", "WF"))
  )

  # Get hurdle (hu) predictions - probability of zero
  pred_hu <- posterior_epred(fit_brms, newdata = newdata, ndraws = ndraws, dpar = "hu")

  # hu is P(Y=0), so P(Y>0) = 1 - hu
  1 - colMeans(pred_hu)
}
