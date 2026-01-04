# =============================================================================
# 08_spatial_simulation_brms.R
#
# Monte Carlo spatial damage simulation using brms hierarchical model.
# Full posterior uncertainty propagation.
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(sf)
library(brms)
library(truncnorm)

# Source helper functions
source(here::here("scripts", "helpers", "brms_predict.R"))
source(here::here("helperfuncties", "pgv_model_original.R"))
source(here::here("helperfuncties", "fetch_knmi_events.R"))

# Load saved models and functions
load(here::here("outputs", "models", "building_mapping_functions.RData"))

# Load brms model
fit_brms <- readRDS(here::here("outputs", "models", "brms_hurdle_gamma_A.rds"))

# Load buildings
buildings <- readRDS(here::here("outputs", "models", "buildings_with_vs30.rds"))

fig_dir <- here::here("outputs", "figures")

cat("Loaded", nrow(buildings), "buildings\n")
cat("Loaded brms model with", nrow(as.data.frame(fit_brms)), "posterior draws\n")

# -----------------------------------------------------------------------------
# 1. Fetch earthquake catalogue
# -----------------------------------------------------------------------------

cat("\n=== Fetching Earthquake Catalogue ===\n")

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

# Filter to significant events
earthquakes <- earthquakes |>
  filter(magnitude >= 1.5) |>
  arrange(time_utc)

cat("After M >= 1.5 filter:", nrow(earthquakes), "events\n")

# -----------------------------------------------------------------------------
# 2. Prepare building sample
# -----------------------------------------------------------------------------

cat("\n=== Preparing Building Sample ===\n")

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

# Distance calculation
calc_epi_distance <- function(eq_lat, eq_lon, bldg_lat, bldg_lon) {
  R <- 6371
  lat1 <- eq_lat * pi / 180
  lat2 <- bldg_lat * pi / 180
  dlat <- (bldg_lat - eq_lat) * pi / 180
  dlon <- (bldg_lon - eq_lon) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  R * c
}

# Sample buildings for simulation
set.seed(42)
sample_buildings <- buildings_df |>
  slice_sample(n = min(500, nrow(buildings_df)))  # Smaller sample for brms speed

cat("Simulating for", nrow(sample_buildings), "buildings and", nrow(earthquakes), "events\n")

# -----------------------------------------------------------------------------
# 3. Calculate PGV matrix
# -----------------------------------------------------------------------------

cat("\n=== Calculating PGV at Buildings ===\n")

n_buildings <- nrow(sample_buildings)
n_events <- nrow(earthquakes)

pgv_matrix <- matrix(NA, nrow = n_buildings, ncol = n_events)
colnames(pgv_matrix) <- earthquakes$event_id

pb <- txtProgressBar(min = 0, max = n_events, style = 3)

for (j in 1:n_events) {
  eq <- earthquakes[j, ]

  distances <- calc_epi_distance(
    eq$lat, eq$lon,
    sample_buildings$lat, sample_buildings$lon
  )

  hypo_distances <- epi_to_hypo_dist(distances, eq$depth_km)

  pgv_result <- pgv_model(
    M = eq$magnitude,
    Rhyp = hypo_distances,
    comp = "LRG",
    VS30 = sample_buildings$vs30
  )

  pgv_matrix[, j] <- pgv_result$PGV * 10  # cm/s to mm/s

  setTxtProgressBar(pb, j)
}
close(pb)

cat("PGV matrix computed. Range:", range(pgv_matrix), "mm/s\n")

# -----------------------------------------------------------------------------
# 4. Monte Carlo damage simulation with brms
# -----------------------------------------------------------------------------

cat("\n=== Running Monte Carlo Damage Simulation (brms) ===\n")

n_sims <- 50  # Fewer sims since brms is slower but has posterior uncertainty

# Generate vulnerability parameters for each building
cat("Generating vulnerability samples...\n")
vuln_params <- generate_vulnerability_batch(sample_buildings, n_sims)

# Initialize damage state matrices
psi_trajectories <- array(NA, dim = c(n_buildings, n_events + 1, n_sims))

cat("Running simulations with brms posterior sampling...\n")
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (sim in 1:n_sims) {
  # Extract vulnerability parameters for this simulation
  material_sim <- vuln_params$material[, sim]
  facade_sim <- vuln_params$facade_type[, sim]
  soil_sim <- vuln_params$soil_profile[, sim]
  psi0_sim <- vuln_params$initial_psi[, sim]

  # Initialize damage state
  psi_current <- psi0_sim
  psi_trajectories[, 1, sim] <- psi_current

  # Process events chronologically
  for (j in 1:n_events) {
    pgv_event <- pgv_matrix[, j]

    # Sample damage increment from brms posterior
    delta_psi <- simulate_damage_increment_brms(
      fit_brms = fit_brms,
      pgv = pgv_event,
      current_psi = psi_current,
      material = material_sim,
      facade_type = facade_sim,
      soil_profile = soil_sim,
      ndraws = 1  # Single draw per event per sim
    )

    # Update damage state
    psi_current <- psi_current + delta_psi
    psi_trajectories[, j + 1, sim] <- psi_current
  }

  setTxtProgressBar(pb, sim)
}
close(pb)

cat("Simulation complete.\n")

# Final damage state
psi_final <- psi_trajectories[, n_events + 1, ]

# -----------------------------------------------------------------------------
# 5. Summarize results
# -----------------------------------------------------------------------------

cat("\n=== Summarizing Results ===\n")

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

cat("\nFinal damage summary (Psi) using brms:\n")
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
    sd_psi = mean(psi_sd),
    p_visible = mean(p_visible_damage),
    .groups = "drop"
  ) |>
  print()

# -----------------------------------------------------------------------------
# 6. Visualize spatial results
# -----------------------------------------------------------------------------

cat("\n=== Generating Spatial Visualizations ===\n")

results_sf <- sample_buildings |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Map of expected damage
p_map_damage <- ggplot(results_sf) +
  geom_sf(aes(color = psi_mean), size = 1) +
  scale_color_viridis_c(option = "C", limits = c(0, 3)) +
  labs(
    title = "Expected Final Damage (Psi) - brms Model",
    subtitle = sprintf("After %d earthquakes (M >= 1.5), %d simulations", n_events, n_sims),
    color = expression(E[Psi])
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_expected_damage_brms.png"), p_map_damage,
       width = 10, height = 8, dpi = 150)
cat("Saved: spatial_expected_damage_brms.png\n")

# Map of damage uncertainty
p_map_uncertainty <- ggplot(results_sf) +
  geom_sf(aes(color = psi_sd), size = 1) +
  scale_color_viridis_c(option = "B") +
  labs(
    title = "Damage Uncertainty (SD) - brms Model",
    subtitle = "Posterior predictive standard deviation",
    color = expression(SD[Psi])
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_damage_uncertainty_brms.png"), p_map_uncertainty,
       width = 10, height = 8, dpi = 150)
cat("Saved: spatial_damage_uncertainty_brms.png\n")

# Map of P(visible damage)
p_map_prob <- ggplot(results_sf) +
  geom_sf(aes(color = p_visible_damage), size = 1) +
  scale_color_viridis_c(option = "D", limits = c(0, 1)) +
  labs(
    title = "Probability of Visible Damage (Psi >= 1) - brms",
    subtitle = sprintf("After %d earthquakes", n_events),
    color = expression(P(Psi >= 1))
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_prob_visible_damage_brms.png"), p_map_prob,
       width = 10, height = 8, dpi = 150)
cat("Saved: spatial_prob_visible_damage_brms.png\n")

# Distribution of final damage with uncertainty bands
p_hist <- sample_buildings |>
  ggplot(aes(x = psi_mean)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1.1, y = Inf, label = "Visible damage",
           color = "red", vjust = 2, hjust = 0) +
  labs(
    title = "Distribution of Expected Final Damage - brms",
    x = expression(E[Psi]),
    y = "Count"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_damage_histogram_brms.png"), p_hist,
       width = 8, height = 5, dpi = 150)
cat("Saved: spatial_damage_histogram_brms.png\n")

# -----------------------------------------------------------------------------
# 7. Aggregate to postcode level
# -----------------------------------------------------------------------------

cat("\n=== Aggregating to Postcode Level ===\n")

pc4_summary <- sample_buildings |>
  mutate(pc4 = substr(postcode, 1, 4)) |>
  group_by(pc4) |>
  summarise(
    n_buildings = n(),
    mean_psi = mean(psi_mean),
    sd_psi = mean(psi_sd),
    p_visible = mean(p_visible_damage),
    mean_bouwjaar = mean(bouwjaar, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nTop 10 postcodes by damage:\n")
pc4_summary |>
  arrange(desc(mean_psi)) |>
  head(10) |>
  print()

saveRDS(pc4_summary, here::here("outputs", "models", "pc4_damage_summary_brms.rds"))

# -----------------------------------------------------------------------------
# 8. Save full results
# -----------------------------------------------------------------------------

saveRDS(sample_buildings, here::here("outputs", "models", "building_damage_results_brms.rds"))
saveRDS(psi_trajectories, here::here("outputs", "models", "psi_trajectories_brms.rds"))

cat("\nResults saved to outputs/models/\n")
cat("Spatial simulation with brms complete.\n")
