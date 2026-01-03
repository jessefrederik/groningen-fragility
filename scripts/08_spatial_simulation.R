# =============================================================================
# 08_spatial_simulation.R
#
# Monte Carlo spatial damage simulation over earthquake catalogue.
# Integrates GMM, vulnerability sampling, and damage accumulation.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(sf)
library(mgcv)
library(truncnorm)

# Source helper functions
source(here::here("helperfuncties", "pgv_model_original.R"))
source(here::here("helperfuncties", "fetch_knmi_events.R"))

# Load saved models and functions
load(here::here("outputs", "models", "accumulation_functions.RData"))
load(here::here("outputs", "models", "building_mapping_functions.RData"))

m1_S <- readRDS(here::here("outputs", "models", "accumulation_stage1_S.rds"))
m2_S <- readRDS(here::here("outputs", "models", "accumulation_stage2_S.rds"))
accum_params <- readRDS(here::here("outputs", "models", "accumulation_params.rds"))

# Load buildings
buildings <- readRDS(here::here("outputs", "models", "buildings_with_vs30.rds"))

fig_dir <- here::here("outputs", "figures")

cat("Loaded", nrow(buildings), "buildings\n")

# -----------------------------------------------------------------------------
# 1. Fetch earthquake catalogue
# -----------------------------------------------------------------------------

cat("\n=== Fetching Earthquake Catalogue ===\n")

# Fetch recent induced seismicity events
# Note: Adjust date range as needed
earthquakes <- tryCatch({
  fetch_knmi_events(
    eventtype = "induced",
    starttime = "2015-01-01",
    endtime = "2024-12-31",
    limit = 500
  )
}, error = function(e) {
  cat("Error fetching KNMI data:", e$message, "\n")
  cat("Using example earthquake catalogue instead.\n")

  # Example catalogue if API fails
  tibble(
    event_id = paste0("EQ", 1:20),
    time_utc = seq(as.POSIXct("2020-01-01"), as.POSIXct("2023-12-31"), length.out = 20),
    lat = runif(20, 53.1, 53.5),
    lon = runif(20, 6.5, 7.0),
    depth_km = runif(20, 2.5, 3.5),
    magnitude = runif(20, 1.5, 3.5),
    location = "Groningen",
    event_type = "induced"
  )
})

cat("Loaded", nrow(earthquakes), "earthquake events\n")
cat("Magnitude range:", range(earthquakes$magnitude), "\n")
cat("Date range:", as.character(range(earthquakes$time_utc)), "\n")

# Filter to significant events (M >= 1.5)
earthquakes <- earthquakes |>
  filter(magnitude >= 1.5) |>
  arrange(time_utc)

cat("After M >= 1.5 filter:", nrow(earthquakes), "events\n")

# -----------------------------------------------------------------------------
# 2. Calculate PGV at each building for each earthquake
# -----------------------------------------------------------------------------

cat("\n=== Calculating PGV at Buildings ===\n")

# Get building coordinates
building_coords <- buildings |>
  st_centroid() |>
  st_coordinates() |>
  as_tibble() |>
  rename(lon = X, lat = Y)

buildings_df <- buildings |>
  st_drop_geometry() |>
  bind_cols(building_coords) |>
  select(identificatie, bouwjaar, oppervlakte, gebruiksdoel, vs30, postcode, lon, lat)

# Function to calculate distance between earthquake and buildings
calc_epi_distance <- function(eq_lat, eq_lon, bldg_lat, bldg_lon) {
  # Haversine formula for distance in km
  R <- 6371  # Earth radius in km

  lat1 <- eq_lat * pi / 180
  lat2 <- bldg_lat * pi / 180
  dlat <- (bldg_lat - eq_lat) * pi / 180
  dlon <- (bldg_lon - eq_lon) * pi / 180

  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))

  R * c
}

# Calculate PGV for a subset of buildings (for demonstration)
# In production, you'd do all buildings or use spatial filtering
set.seed(42)
sample_buildings <- buildings_df |>
  slice_sample(n = min(1000, nrow(buildings_df)))

cat("Simulating for", nrow(sample_buildings), "buildings and", nrow(earthquakes), "events\n")

# Pre-compute PGV matrix: buildings x earthquakes
n_buildings <- nrow(sample_buildings)
n_events <- nrow(earthquakes)

pgv_matrix <- matrix(NA, nrow = n_buildings, ncol = n_events)
colnames(pgv_matrix) <- earthquakes$event_id

cat("Computing PGV matrix...\n")
pb <- txtProgressBar(min = 0, max = n_events, style = 3)

for (j in 1:n_events) {
  eq <- earthquakes[j, ]

  # Distance from earthquake to each building
  distances <- calc_epi_distance(
    eq$lat, eq$lon,
    sample_buildings$lat, sample_buildings$lon
  )

  # Hypocentral distance
  hypo_distances <- epi_to_hypo_dist(distances, eq$depth_km)

  # PGV calculation (using LRG component)
  pgv_result <- pgv_model(
    M = eq$magnitude,
    Rhyp = hypo_distances,
    comp = "LRG",
    VS30 = sample_buildings$vs30
  )

  # Convert from cm/s to mm/s
  pgv_matrix[, j] <- pgv_result$PGV * 10

  setTxtProgressBar(pb, j)
}
close(pb)

cat("PGV matrix computed.\n")
cat("PGV range across all building-events:", range(pgv_matrix), "mm/s\n")

# -----------------------------------------------------------------------------
# 3. Monte Carlo damage simulation
# -----------------------------------------------------------------------------

cat("\n=== Running Monte Carlo Damage Simulation ===\n")

n_sims <- 100  # Monte Carlo simulations
gamma_shape <- accum_params$gamma_shape
alpha <- accum_params$alpha
pgv_ref <- accum_params$pgv_ref

# Generate vulnerability parameters for each building
cat("Generating vulnerability samples...\n")
vuln_params <- generate_vulnerability_batch(sample_buildings, n_sims)

# Initialize damage state matrix: buildings x simulations
psi_final <- matrix(NA, nrow = n_buildings, ncol = n_sims)

cat("Running simulations...\n")
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (sim in 1:n_sims) {
  # Extract vulnerability parameters for this simulation
  material_sim <- vuln_params$material[, sim]
  facade_sim <- vuln_params$facade_type[, sim]
  soil_sim <- vuln_params$soil_profile[, sim]
  psi0_sim <- vuln_params$initial_psi[, sim]

  # Initialize damage state
  psi_current <- psi0_sim

  # Process events chronologically
  for (j in 1:n_events) {
    pgv_event <- pgv_matrix[, j]

    # Calculate cumulative S up to this event
    # For simplicity, we treat each event independently (not truly cumulative)
    # In a full implementation, you'd track S per building
    S_event <- (pgv_event / pgv_ref)^alpha

    # Prepare prediction data
    pred_data <- tibble(
      log_S = log(S_event + 0.01),
      InitialPsi = psi_current,
      Material = material_sim,
      FacadeType = facade_sim,
      SoilProfile = soil_sim,
      EarthquakeType = "ZN"  # Placeholder
    )

    # Predict damage
    p_damage <- predict(m1_S, newdata = pred_data, type = "response")
    mu_damage <- predict(m2_S, newdata = pred_data, type = "response")

    # Simulate damage occurrence and amount
    damage_occurs <- rbinom(n_buildings, 1, p_damage)
    damage_amount <- rgamma(n_buildings, shape = gamma_shape, rate = gamma_shape / mu_damage)
    delta_psi <- damage_occurs * damage_amount

    # Update damage state
    psi_current <- psi_current + delta_psi
  }

  psi_final[, sim] <- psi_current
  setTxtProgressBar(pb, sim)
}
close(pb)

cat("Simulation complete.\n")

# -----------------------------------------------------------------------------
# 4. Summarize results
# -----------------------------------------------------------------------------

cat("\n=== Summarizing Results ===\n")

# Compute summary statistics per building
sample_buildings <- sample_buildings |>
  mutate(
    psi_mean = rowMeans(psi_final),
    psi_median = apply(psi_final, 1, median),
    psi_sd = apply(psi_final, 1, sd),
    psi_q05 = apply(psi_final, 1, quantile, 0.05),
    psi_q95 = apply(psi_final, 1, quantile, 0.95),
    p_visible_damage = rowMeans(psi_final >= 1),
    p_moderate_damage = rowMeans(psi_final >= 2)
  )

cat("\nFinal damage summary (Psi):\n")
cat("Mean:", mean(sample_buildings$psi_mean), "\n")
cat("Median:", median(sample_buildings$psi_median), "\n")
cat("P(Psi >= 1) overall:", mean(sample_buildings$p_visible_damage), "\n")
cat("P(Psi >= 2) overall:", mean(sample_buildings$p_moderate_damage), "\n")

# Summary by building age
cat("\nDamage by construction period:\n")
sample_buildings |>
  mutate(
    period = cut(bouwjaar,
      breaks = c(0, 1920, 1960, 1990, 2025),
      labels = c("<1920", "1920-60", "1960-90", ">1990"))
  ) |>
  group_by(period) |>
  summarise(
    n = n(),
    mean_psi = mean(psi_mean),
    p_visible = mean(p_visible_damage)
  ) |>
  print()

# -----------------------------------------------------------------------------
# 5. Visualize spatial results
# -----------------------------------------------------------------------------

cat("\n=== Generating Spatial Visualizations ===\n")

# Join results back to spatial data
results_sf <- sample_buildings |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Map of expected damage
p_map_damage <- ggplot(results_sf) +
  geom_sf(aes(color = psi_mean), size = 0.5) +
  scale_color_viridis_c(option = "C", limits = c(0, 2)) +
  labs(
    title = "Expected Final Damage (Psi)",
    subtitle = sprintf("After %d earthquakes (M >= 1.5)", n_events),
    color = expression(E[Psi])
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_expected_damage.png"), p_map_damage,
       width = 10, height = 8, dpi = 150)
cat("Saved: spatial_expected_damage.png\n")

# Map of P(visible damage)
p_map_prob <- ggplot(results_sf) +
  geom_sf(aes(color = p_visible_damage), size = 0.5) +
  scale_color_viridis_c(option = "B", limits = c(0, 1)) +
  labs(
    title = "Probability of Visible Damage (Psi >= 1)",
    subtitle = sprintf("After %d earthquakes", n_events),
    color = expression(P(Psi >= 1))
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_prob_visible_damage.png"), p_map_prob,
       width = 10, height = 8, dpi = 150)
cat("Saved: spatial_prob_visible_damage.png\n")

# Histogram of final damage
p_hist <- sample_buildings |>
  ggplot(aes(x = psi_mean)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed",
             linewidth = 1) +
  annotate("text", x = 1.1, y = Inf, label = "Visible damage",
           color = "red", vjust = 2, hjust = 0) +
  labs(
    title = "Distribution of Expected Final Damage",
    x = expression(E[Psi]),
    y = "Count"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_damage_histogram.png"), p_hist,
       width = 8, height = 5, dpi = 150)
cat("Saved: spatial_damage_histogram.png\n")

# Damage by distance from earthquake centroid
eq_centroid <- earthquakes |>
  summarise(
    lat = mean(lat),
    lon = mean(lon)
  )

sample_buildings <- sample_buildings |>
  mutate(
    dist_to_centroid = calc_epi_distance(
      eq_centroid$lat, eq_centroid$lon,
      lat, lon
    )
  )

p_dist <- sample_buildings |>
  ggplot(aes(x = dist_to_centroid, y = psi_mean)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  labs(
    title = "Expected Damage vs Distance from Earthquake Centroid",
    x = "Distance (km)",
    y = expression(E[Psi])
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_damage_vs_distance.png"), p_dist,
       width = 8, height = 5, dpi = 150)
cat("Saved: spatial_damage_vs_distance.png\n")

# -----------------------------------------------------------------------------
# 6. Aggregate to postcode level
# -----------------------------------------------------------------------------

cat("\n=== Aggregating to Postcode Level ===\n")

pc4_summary <- sample_buildings |>
  mutate(pc4 = substr(postcode, 1, 4)) |>
  group_by(pc4) |>
  summarise(
    n_buildings = n(),
    mean_psi = mean(psi_mean),
    sd_psi = sd(psi_mean),
    p_visible = mean(p_visible_damage),
    mean_bouwjaar = mean(bouwjaar, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nTop 10 postcodes by damage:\n")
pc4_summary |>
  arrange(desc(mean_psi)) |>
  head(10) |>
  print()

# Save PC4 summary
saveRDS(pc4_summary, here::here("outputs", "models", "pc4_damage_summary.rds"))

# -----------------------------------------------------------------------------
# 7. Save full results
# -----------------------------------------------------------------------------

# Save building-level results
saveRDS(sample_buildings, here::here("outputs", "models", "building_damage_results.rds"))

# Save raw simulation matrix
saveRDS(psi_final, here::here("outputs", "models", "psi_simulation_matrix.rds"))

cat("\nResults saved to outputs/models/\n")
cat("Spatial simulation complete.\n")
