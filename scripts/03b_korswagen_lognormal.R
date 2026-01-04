# =============================================================================
# 03b_korswagen_lognormal.R
#
# Predict cumulative earthquake damage for Groningen buildings using:
# - Korswagen et al. (2022) lognormal fragility curves (Table 10)
# - Bommer et al. (2022) ground motion model with full uncertainty
# - State-dependent damage accumulation
#
# METHODOLOGY:
# 1. For each building and earthquake event:
#    a. Calculate PGV with GMM uncertainty (tau, phiS2S, phiSS)
#    b. Use current Psi to select appropriate fragility parameters
#    c. Sample Delta-Psi from the exceedance probability distribution
#    d. Update: Psi(t+1) = Psi(t) + Delta-Psi
#
# KEY DIFFERENCE FROM BRMS MODEL:
# - Uses published lognormal fragility (Table 10) instead of fitted hurdle-gamma
# - Direct physics-based relationship from Korswagen's Monte Carlo simulations
# - No posterior uncertainty, but includes GMM uncertainty convolution
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(sf)
library(here)
library(matrixStats)

# Source helper functions
source(here("helperfuncties", "fetch_knmi_events.R"))
source(here("helperfuncties", "pgv_model_original.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

PSI_MAX <- 3.0          # Max damage (from FEM data range in Korswagen)
MAG_THRESHOLD <- 2.0    # Minimum magnitude
PGV_THRESHOLD <- 2.0    # mm/s - lower bound for potential damage
L_GMM <- 10             # Number of GMM realizations
SAMPLE_FRAC <- 1.0      # Full dataset (set to 0.1 for testing)

cat("=== Korswagen Lognormal Fragility Analysis ===\n")
cat("Configuration:\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  PGV threshold:", PGV_THRESHOLD, "mm/s\n")
cat("  Max Psi:", PSI_MAX, "\n")
cat("  Sample fraction:", SAMPLE_FRAC, "\n\n")

# =============================================================================
# KORSWAGEN TABLE 10: LOGNORMAL FRAGILITY PARAMETERS
# =============================================================================
#
# From Korswagen et al. (2022) Table 10, pp. 6218
# mu is the lognormal location parameter in ln(PGV) space
# sigma is the lognormal shape parameter (log-standard deviation)
#
# Interpretation: P(Psi_final >= threshold | PGV, Psi_0) = pnorm((ln(PGV) - mu) / sigma)
#
# =============================================================================

korswagen_params <- tribble(
  ~psi_0, ~psi_final, ~mu, ~sigma,
  # Psi_0 = 0 (undamaged)
  0.0, 0.5, 2.604, 0.624,
  0.0, 1.0, 3.370, 0.725,
  0.0, 1.5, 4.323, 0.860,
  0.0, 2.0, 5.518, 1.082,
  0.0, 2.5, 7.416, 1.516,
  0.0, 3.0, 8.869, 1.729,
  # Psi_0 = 0.5 (light imperceptible damage)
  0.5, 1.0, 3.002, 0.931,
  0.5, 1.5, 4.070, 0.983,
  0.5, 2.0, 5.523, 1.259,
  0.5, 2.5, 7.758, 1.780,
  0.5, 3.0, 11.287, 2.606,
  # Psi_0 = 1.0 (visible damage)
  1.0, 1.5, 3.411, 1.034,
  1.0, 2.0, 5.293, 1.444,
  1.0, 2.5, 7.955, 2.091,
  1.0, 3.0, 11.913, 3.034,
  # Psi_0 = 1.5 (clearly visible damage)
  1.5, 2.0, 4.631, 1.860,
  1.5, 2.5, 7.929, 2.585,
  1.5, 3.0, 12.363, 3.610
)

cat("Korswagen Table 10 parameters loaded:\n")
cat("  Rows:", nrow(korswagen_params), "\n")
cat("  Psi_0 levels: 0, 0.5, 1.0, 1.5\n")
cat("  Psi_final range: 0.5 to 3.0\n\n")

# Show median PGV for key thresholds
cat("Median PGV (mm/s) for key transitions:\n")
korswagen_params |>
  filter(psi_0 == 0) |>
  mutate(median_pgv = exp(mu)) |>
  select(psi_final, median_pgv) |>
  print()
cat("\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Get the discrete Psi_0 level for parameter lookup
#'
#' Floors to nearest discrete level (0, 0.5, 1.0, 1.5) and caps at 1.5
get_discrete_psi0 <- function(psi) {
  psi_capped <- pmin(psi, 1.5)
  floor(psi_capped * 2) / 2
}

#' Sample Delta-Psi from Korswagen lognormal fragility
#'
#' Uses the step-function inverse CDF approach:
#' 1. Calculate exceedance probability for each threshold above current psi
#' 2. Sample U ~ Uniform(0,1)
#' 3. Find highest threshold where P(exceed) >= U
#' 4. Delta-Psi = threshold - current_psi
#'
#' @param pgv_mm_s PGV in mm/s
#' @param current_psi Current damage state
#' @param params_table Korswagen parameter table
#' @return Sampled Delta-Psi
sample_delta_psi_korswagen <- function(pgv_mm_s, current_psi, params_table = korswagen_params) {
  # Get discrete Psi_0 for parameter lookup
  psi_0_discrete <- get_discrete_psi0(current_psi)

  # Get available final states (must exceed current psi)
  available <- params_table |>
    filter(psi_0 == psi_0_discrete, psi_final > current_psi) |>
    arrange(psi_final)

  if (nrow(available) == 0) {
    return(0)  # At or beyond maximum damage state
  }

  # Calculate exceedance probabilities for each level
  ln_pgv <- log(pgv_mm_s)
  p_exceed <- pnorm((ln_pgv - available$mu) / available$sigma)

  # Sample: find highest level exceeded
  u <- runif(1)

  # Find all levels where P(exceed) >= u
  exceeded_idx <- which(p_exceed >= u)

  if (length(exceeded_idx) == 0) {
    return(0)  # No damage increment
  }

  # Take the highest exceeded level
  psi_final_sampled <- available$psi_final[max(exceeded_idx)]
  delta_psi <- psi_final_sampled - current_psi

  return(max(0, delta_psi))
}

#' Vectorized version for multiple buildings
#'
#' @param pgv_vec Vector of PGV values (mm/s)
#' @param psi_vec Vector of current damage states
#' @return Vector of sampled Delta-Psi values
sample_delta_psi_vec <- function(pgv_vec, psi_vec, params_table = korswagen_params) {
  n <- length(pgv_vec)
  delta <- numeric(n)
  ln_pgv <- log(pgv_vec)

  # Process by discrete Psi_0 level for efficiency
  psi_0_levels <- c(0, 0.5, 1.0, 1.5)

  for (psi_0 in psi_0_levels) {
    # Find buildings in this Psi_0 bracket
    psi_discrete <- get_discrete_psi0(psi_vec)
    mask <- psi_discrete == psi_0

    if (!any(mask)) next

    # Get parameters for this Psi_0
    params_psi0 <- params_table |> filter(psi_0 == !!psi_0)

    if (nrow(params_psi0) == 0) next

    # For each building in this bracket
    idx_buildings <- which(mask)
    for (j in idx_buildings) {
      # Filter to levels above current psi
      valid_levels <- params_psi0 |> filter(psi_final > psi_vec[j])
      if (nrow(valid_levels) == 0) next

      # Exceedance probabilities
      p_exceed <- pnorm((ln_pgv[j] - valid_levels$mu) / valid_levels$sigma)

      # Sample
      u <- runif(1)
      exceeded_idx <- which(p_exceed >= u)

      if (length(exceeded_idx) > 0) {
        psi_final <- valid_levels$psi_final[max(exceeded_idx)]
        delta[j] <- psi_final - psi_vec[j]
      }
    }
  }

  return(delta)
}

# =============================================================================
# 1. LOAD EARTHQUAKE CATALOGUE
# =============================================================================

cat("Step 1: Loading earthquake catalogue...\n")

eq_cache <- here("outputs", "models", "groningen_earthquakes.rds")
if (file.exists(eq_cache)) {
  groningen_eq <- readRDS(eq_cache)
  cat("Loaded cached catalogue:", nrow(groningen_eq), "earthquakes\n")
} else {
  earthquakes <- fetch_knmi_events(
    eventtype = "induced or triggered event",
    starttime = "1991-01-01",
    endtime = "2024-12-31"
  )
  groningen_eq <- earthquakes |>
    filter(lat >= 53.0, lat <= 53.6, lon >= 6.4, lon <= 7.3, mag >= MAG_THRESHOLD) |>
    arrange(time_utc)
  saveRDS(groningen_eq, eq_cache)
  cat("Fetched and cached:", nrow(groningen_eq), "earthquakes\n")
}

n_events <- nrow(groningen_eq)
cat("Earthquakes (M >=", MAG_THRESHOLD, "):", n_events, "\n\n")

# =============================================================================
# 2. LOAD BUILDINGS
# =============================================================================

cat("Step 2: Loading building data...\n")

buildings <- st_read(here("datafiles", "bag_selectie.gpkg"), quiet = TRUE)
n_buildings_full <- nrow(buildings)
cat("Total buildings:", format(n_buildings_full, big.mark = ","), "\n")

# Load Vs30
vs30 <- read_delim(
  here("datafiles", "vs30.csv"), delim = ";",
  locale = locale(decimal_mark = ","),
  col_names = c("postcode", "vs30", "x1", "x2", "x3"),
  skip = 1, show_col_types = FALSE
) |> select(postcode, vs30) |> mutate(postcode = as.character(postcode))

# Join Vs30 and fill missing
buildings <- buildings |>
  mutate(pc4 = substr(postcode, 1, 4)) |>
  left_join(vs30, by = c("pc4" = "postcode"))

mean_vs30 <- mean(buildings$vs30, na.rm = TRUE)
buildings$vs30[is.na(buildings$vs30)] <- mean_vs30

# Transform to WGS84 and extract centroids
buildings_wgs84 <- st_transform(buildings, 4326)
centroids <- st_centroid(st_geometry(buildings_wgs84))
coords <- st_coordinates(centroids)
buildings$lon <- coords[, 1]
buildings$lat <- coords[, 2]

# Distance helper (Haversine)
calc_dist_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  lat1_rad <- lat1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# Distance to epicenter cluster (Loppersum)
buildings$dist_to_center <- calc_dist_km(buildings$lat, buildings$lon, 53.35, 6.72)

# Stratified sampling for validation
if (SAMPLE_FRAC < 1.0) {
  cat("Stratified sampling:", SAMPLE_FRAC * 100, "% of buildings\n")

  set.seed(123)

  # Stratify by distance zone
  near <- buildings |> filter(dist_to_center < 10)
  mid <- buildings |> filter(dist_to_center >= 10, dist_to_center < 25)
  far <- buildings |> filter(dist_to_center >= 25)

  n_sample <- round(n_buildings_full * SAMPLE_FRAC)
  n_near <- min(round(n_sample * 0.4), nrow(near))
  n_mid <- min(round(n_sample * 0.4), nrow(mid))
  n_far <- min(n_sample - n_near - n_mid, nrow(far))

  buildings <- bind_rows(
    slice_sample(near, n = n_near),
    slice_sample(mid, n = n_mid),
    slice_sample(far, n = n_far)
  )
  cat("  Near field:", n_near, "| Mid field:", n_mid, "| Far field:", n_far, "\n")
}

n_buildings <- nrow(buildings)
cat("Buildings for analysis:", format(n_buildings, big.mark = ","), "\n")

# Add metadata columns
buildings <- buildings |>
  mutate(
    bouwjaar_numeric = as.numeric(bouwjaar),
    age_category = case_when(
      bouwjaar_numeric < 1920 ~ "Pre-1920",
      bouwjaar_numeric < 1950 ~ "1920-1949",
      bouwjaar_numeric < 1970 ~ "1950-1969",
      bouwjaar_numeric < 1991 ~ "1970-1990",
      TRUE ~ "Post-1991"
    ),
    zone = case_when(
      dist_to_center < 10 ~ "Near (<10km)",
      dist_to_center < 25 ~ "Mid (10-25km)",
      TRUE ~ "Far (>25km)"
    )
  )

buildings_df <- st_drop_geometry(buildings)
cat("Buildings prepared\n\n")

# =============================================================================
# 3. GMM UNCERTAINTY COMPONENTS
# =============================================================================

cat("Step 3: Setting up GMM uncertainty...\n")

gmm_coef <- get_pgv_coefficients("LRG")
gmm_tau <- gmm_coef$tau        # Inter-event (0.2448)
gmm_phiS2S <- gmm_coef$phiS2S  # Site-to-site (0.2406)
gmm_phiSS <- gmm_coef$phiSS    # Within-event (0.4569)

cat("  tau (inter-event):", round(gmm_tau, 4), "\n")
cat("  phiS2S (site-to-site):", round(gmm_phiS2S, 4), "\n")
cat("  phiSS (within-event):", round(gmm_phiSS, 4), "\n\n")

# =============================================================================
# 4. PRE-COMPUTE MEDIAN PGV MATRIX
# =============================================================================

cat("Step 4: Computing median PGV matrix...\n")

# Sort events chronologically
event_order <- order(groningen_eq$time_utc)
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

# Pre-compute median lnPGV for all events (N x J matrix)
lnPGV_median_all <- matrix(NA, n_buildings, n_events)
for (j in 1:n_events) {
  epi_dist <- calc_dist_km(buildings_df$lat, buildings_df$lon, eq_lat[j], eq_lon[j])
  hypo_dist <- sqrt(epi_dist^2 + eq_depth[j]^2)
  pgv_result <- pgv_model(M = eq_mag[j], Rhyp = hypo_dist, comp = "LRG", VS30 = buildings_df$vs30)
  lnPGV_median_all[, j] <- pgv_result$lnPGV
}

cat("  PGV matrix size:", n_buildings, "x", n_events, "\n")
cat("  Memory:", round(object.size(lnPGV_median_all) / 1e6, 1), "MB\n\n")

# =============================================================================
# 5. MAIN SIMULATION LOOP WITH GMM UNCERTAINTY
# =============================================================================

cat("Step 5: Running Korswagen lognormal simulation...\n")
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("  Events:", n_events, "\n")
cat("  GMM realizations (L):", L_GMM, "\n\n")

# Storage for final Psi across all GMM realizations
psi_final_all <- matrix(NA, n_buildings, L_GMM)
max_pgv_all <- matrix(NA, n_buildings, L_GMM)
n_damaging_events_all <- matrix(0, n_buildings, L_GMM)

set.seed(42)
global_start <- Sys.time()

# =============================================================================
# OUTER LOOP: L GMM REALIZATIONS
# =============================================================================

for (l in 1:L_GMM) {
  cat(sprintf("GMM realization %d/%d", l, L_GMM))
  l_start <- Sys.time()

  # Generate GMM uncertainty terms for this realization
  # 1. Site terms (phiS2S) - one per building, constant across events
  eta_site <- rnorm(n_buildings, 0, gmm_phiS2S)

  # 2. Event terms (tau) - one per earthquake
  eta_event <- rnorm(n_events, 0, gmm_tau)

  # Initialize damage state - all buildings start at Psi = 0 (virgin walls)
  psi <- rep(0, n_buildings)
  max_pgv <- rep(0, n_buildings)
  n_damaging <- rep(0, n_buildings)

  # INNER LOOP: Events in chronological order
  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]

    # 3. Within-event terms (phiSS) - different for each building-event
    eta_within <- rnorm(n_buildings, 0, gmm_phiSS)

    # Realized PGV for this GMM realization
    lnPGV_realized <- lnPGV_median_all[, j] + eta_event[j] + eta_site + eta_within
    pgv_event <- exp(lnPGV_realized) * 10  # Convert cm/s to mm/s

    # Update max PGV
    max_pgv <- pmax(max_pgv, pgv_event)

    # Find affected buildings (above threshold)
    affected <- which(pgv_event >= PGV_THRESHOLD)

    if (length(affected) > 0) {
      # Sample damage increment for affected buildings
      delta_psi <- sample_delta_psi_vec(
        pgv_vec = pgv_event[affected],
        psi_vec = psi[affected]
      )

      # Update damage state
      psi[affected] <- pmin(psi[affected] + delta_psi, PSI_MAX)

      # Count damaging events (those with delta > 0)
      n_damaging[affected] <- n_damaging[affected] + (delta_psi > 0)
    }
  }

  # Store results for this GMM realization
  psi_final_all[, l] <- psi
  max_pgv_all[, l] <- max_pgv
  n_damaging_events_all[, l] <- n_damaging

  l_elapsed <- as.numeric(difftime(Sys.time(), l_start, units = "secs"))
  cat(sprintf(" - %.1f sec\n", l_elapsed))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\nSimulation complete! Total time: %.1f min\n\n", total_elapsed))

# =============================================================================
# 6. COMPUTE SUMMARY STATISTICS
# =============================================================================

cat("Step 6: Computing summary statistics...\n")

buildings_df <- buildings_df |>
  mutate(
    # Point estimates
    psi_mean = rowMeans(psi_final_all),
    psi_sd = rowSds(psi_final_all),

    # Percentiles (from GMM uncertainty)
    psi_p10 = rowQuantiles(psi_final_all, probs = 0.10),
    psi_p50 = rowQuantiles(psi_final_all, probs = 0.50),
    psi_p90 = rowQuantiles(psi_final_all, probs = 0.90),

    # Exceedance probabilities
    p_visible = rowMeans(psi_final_all >= 1.0),   # P(Psi >= 1)
    p_moderate = rowMeans(psi_final_all >= 2.0),  # P(Psi >= 2)
    p_severe = rowMeans(psi_final_all >= 3.0),    # P(Psi >= 3)

    # Event and PGV statistics
    n_damaging_events = rowMeans(n_damaging_events_all),
    max_pgv = rowMaxs(max_pgv_all),
    max_pgv_mean = rowMeans(max_pgv_all),
    max_pgv_p10 = rowQuantiles(max_pgv_all, probs = 0.10),
    max_pgv_p90 = rowQuantiles(max_pgv_all, probs = 0.90)
  )

# =============================================================================
# 7. RESULTS SUMMARY
# =============================================================================

cat("\n=== RESULTS SUMMARY (Korswagen Lognormal) ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("GMM realizations:", L_GMM, "\n")
cat("Earthquakes:", n_events, "\n\n")

cat("--- DAMAGE STATISTICS ---\n")
cat("Mean Psi:", round(mean(buildings_df$psi_mean), 4), "\n")
cat("  p10:", round(mean(buildings_df$psi_p10), 4),
    "| p50:", round(mean(buildings_df$psi_p50), 4),
    "| p90:", round(mean(buildings_df$psi_p90), 4), "\n")
cat("Mean SD:", round(mean(buildings_df$psi_sd), 4), "\n\n")

cat("--- EXCEEDANCE PROBABILITIES ---\n")
cat("P(Psi >= 1) visible:", round(100 * mean(buildings_df$p_visible), 2), "%\n")
cat("P(Psi >= 2) moderate:", round(100 * mean(buildings_df$p_moderate), 2), "%\n")
cat("P(Psi >= 3) severe:", round(100 * mean(buildings_df$p_severe), 2), "%\n\n")

cat("--- BY DISTANCE ZONE ---\n")
by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    psi_mean = mean(psi_mean),
    psi_p10 = mean(psi_p10),
    psi_p90 = mean(psi_p90),
    p_visible = mean(p_visible),
    max_pgv_mean = mean(max_pgv_mean),
    .groups = "drop"
  )
print(by_zone)

cat("\n--- BY AGE CATEGORY ---\n")
by_age <- buildings_df |>
  group_by(age_category) |>
  summarise(
    n = n(),
    psi_mean = mean(psi_mean),
    p_visible = mean(p_visible),
    .groups = "drop"
  )
print(by_age)

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

cat("\nStep 7: Saving results...\n")

results <- buildings_df |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, zone, age_category,
    # Korswagen lognormal predictions
    psi_mean, psi_sd, psi_p10, psi_p50, psi_p90,
    # Exceedance probabilities
    p_visible, p_moderate, p_severe,
    # Event statistics
    n_damaging_events, max_pgv, max_pgv_mean, max_pgv_p10, max_pgv_p90
  )

output_file <- here("outputs", "models", "korswagen_lognormal_results.rds")
saveRDS(results, output_file)
cat("Saved:", output_file, "\n")

# Summary tables
summary_tables <- list(
  by_zone = by_zone,
  by_age = by_age,
  overall = tibble(
    method = "korswagen_lognormal",
    n_buildings = n_buildings,
    n_events = n_events,
    l_gmm = L_GMM,
    psi_mean = mean(buildings_df$psi_mean),
    psi_p10 = mean(buildings_df$psi_p10),
    psi_p90 = mean(buildings_df$psi_p90),
    psi_sd = mean(buildings_df$psi_sd),
    p_visible = mean(buildings_df$p_visible),
    p_moderate = mean(buildings_df$p_moderate),
    p_severe = mean(buildings_df$p_severe)
  ),
  korswagen_params = korswagen_params
)

summary_file <- here("outputs", "models", "korswagen_lognormal_summaries.rds")
saveRDS(summary_tables, summary_file)
cat("Saved:", summary_file, "\n")

cat("\n=== DONE ===\n")
