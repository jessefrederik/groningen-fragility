# =============================================================================
# 11_spatial_damage_prediction.R
#
# Predict cumulative earthquake damage for Groningen buildings using:
# - KNMI induced earthquake catalogue (1991-2024)
# - Bommer et al. (2022) ground motion model
# - Fitted brms hurdle-gamma fragility model
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(sf)
library(brms)
library(here)

# Source helper functions
source(here("helperfuncties", "fetch_knmi_events.R"))
source(here("helperfuncties", "pgv_model_original.R"))

fig_dir <- here("outputs", "figures")

cat("=== Groningen Spatial Damage Prediction ===\n\n")

# -----------------------------------------------------------------------------
# 1. Fetch KNMI Earthquake Catalogue
# -----------------------------------------------------------------------------

cat("Step 1: Fetching KNMI earthquake catalogue...\n")

earthquakes <- tryCatch({
  fetch_knmi_events(
    eventtype = "induced or triggered event",  # Full QuakeML event type
    starttime = "1991-01-01",
    endtime = "2024-12-31"
  )
}, error = function(e) {
  cat("Warning: Could not fetch KNMI data:", conditionMessage(e), "\n")
  cat("Using synthetic data for testing.\n")
  NULL
})

if (is.null(earthquakes) || nrow(earthquakes) == 0) {
  cat("Creating synthetic earthquake catalogue for testing...\n")
  # Synthetic data based on known major events
  earthquakes <- tibble(
    event_id = paste0("synth_", 1:50),
    time_utc = seq(as.POSIXct("1991-01-01"), as.POSIXct("2024-01-01"), length.out = 50),
    lat = runif(50, 53.2, 53.45),
    lon = runif(50, 6.6, 6.9),
    depth_km = runif(50, 2, 4),
    mag = c(3.6, 3.4, 3.2, 3.0, rep(runif(46, 1.5, 2.8))),  # Include major events
    location = "Groningen"
  )
  # Add known major events
  earthquakes$lat[1] <- 53.35; earthquakes$lon[1] <- 6.67  # Huizinge 2012
  earthquakes$lat[2] <- 53.36; earthquakes$lon[2] <- 6.78  # Zeerijp 2018
}

# Filter Groningen region + magnitude threshold
groningen_eq <- earthquakes |>
  filter(
    lat >= 53.0, lat <= 53.6,
    lon >= 6.4, lon <= 7.3,
    mag >= 2.0
  ) |>
  arrange(time_utc)

cat("Found", nrow(groningen_eq), "induced earthquakes (M >= 2.0)\n")
cat("Date range:", as.character(min(groningen_eq$time_utc)), "to",
    as.character(max(groningen_eq$time_utc)), "\n")
cat("Magnitude range:", round(min(groningen_eq$mag), 1), "to",
    round(max(groningen_eq$mag), 1), "\n\n")

# -----------------------------------------------------------------------------
# 2. Load Buildings and Vs30
# -----------------------------------------------------------------------------

cat("Step 2: Loading building data...\n")

buildings <- st_read(here("datafiles", "bag_selectie.gpkg"), quiet = TRUE)
cat("Loaded", nrow(buildings), "buildings from BAG\n")

# Load Vs30 by postcode
vs30 <- read_delim(
  here("datafiles", "vs30.csv"),
  delim = ";",
  locale = locale(decimal_mark = ","),
  col_names = c("postcode", "vs30", "x1", "x2", "x3"),
  skip = 1,
  show_col_types = FALSE
) |>
  select(postcode, vs30) |>
  mutate(postcode = as.character(postcode))

cat("Loaded Vs30 for", nrow(vs30), "postcodes\n")

# Join Vs30 to buildings
buildings <- buildings |>
  mutate(pc4 = substr(postcode, 1, 4)) |>
  left_join(vs30, by = c("pc4" = "postcode"))

# Fill missing Vs30 with regional mean
mean_vs30 <- mean(buildings$vs30, na.rm = TRUE)
n_missing <- sum(is.na(buildings$vs30))
buildings <- buildings |>
  mutate(vs30 = coalesce(vs30, mean_vs30))

cat("Filled", n_missing, "missing Vs30 values with regional mean:",
    round(mean_vs30, 0), "m/s\n")

# Transform to WGS84 and extract centroids
cat("Transforming coordinates to WGS84...\n")
buildings_wgs84 <- st_transform(buildings, 4326)
centroids <- st_centroid(st_geometry(buildings_wgs84))
coords <- st_coordinates(centroids)
buildings$lon <- coords[, 1]
buildings$lat <- coords[, 2]

cat("Coordinate range: lat", round(min(buildings$lat), 2), "-",
    round(max(buildings$lat), 2), ", lon", round(min(buildings$lon), 2), "-",
    round(max(buildings$lon), 2), "\n\n")

# -----------------------------------------------------------------------------
# 3. Stratified Sampling (enriched for high-damage areas)
# -----------------------------------------------------------------------------

cat("Step 3: Creating stratified building sample...\n")

# Haversine distance function
calc_dist_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  lat1_rad <- lat1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# Epicenter cluster centroid (Loppersum area)
epicenter_lat <- 53.35
epicenter_lon <- 6.72

buildings$dist_to_center <- calc_dist_km(
  buildings$lat, buildings$lon,
  epicenter_lat, epicenter_lon
)

# Stratified sample: oversample near-field
set.seed(42)

near_field <- buildings |> filter(dist_to_center < 10)
mid_field <- buildings |> filter(dist_to_center >= 10, dist_to_center < 25)
far_field <- buildings |> filter(dist_to_center >= 25)

cat("Near-field (<10km):", nrow(near_field), "buildings\n")
cat("Mid-field (10-25km):", nrow(mid_field), "buildings\n")
cat("Far-field (>25km):", nrow(far_field), "buildings\n")

sample_buildings <- bind_rows(
  slice_sample(near_field, n = min(2000, nrow(near_field))),
  slice_sample(mid_field, n = min(2000, nrow(mid_field))),
  slice_sample(far_field, n = min(1000, nrow(far_field)))
) |>
  st_drop_geometry()  # Drop geometry for faster processing

cat("Sampled", nrow(sample_buildings), "buildings for analysis\n\n")

# -----------------------------------------------------------------------------
# 3b. Age-Based Pre-Damage Assignment
# -----------------------------------------------------------------------------

cat("Step 3b: Assigning initial damage based on building age...\n")

# Pre-damage heuristic based on construction year
# Older buildings: higher initial Ψ₀ (prior earthquake exposure + age-related wear)
# Reference: Gas extraction began 1963, major seismicity from ~1991
assign_initial_damage <- function(bouwjaar) {
  case_when(
    bouwjaar < 1920 ~ 0.8,   # Very old: significant prior wear
    bouwjaar < 1950 ~ 0.5,   # Pre-war: moderate prior damage
    bouwjaar < 1970 ~ 0.3,   # Post-war: some exposure
    bouwjaar < 1991 ~ 0.15,  # Built before major seismicity
    TRUE ~ 0                  # Recent construction: virgin
  )
}

sample_buildings <- sample_buildings |>
  mutate(
    bouwjaar_numeric = as.numeric(bouwjaar),
    initial_psi = assign_initial_damage(bouwjaar_numeric),
    age_category = case_when(
      bouwjaar_numeric < 1920 ~ "Pre-1920",
      bouwjaar_numeric < 1950 ~ "1920-1949",
      bouwjaar_numeric < 1970 ~ "1950-1969",
      bouwjaar_numeric < 1991 ~ "1970-1990",
      TRUE ~ "Post-1991"
    )
  )

# Summary of initial damage distribution
age_summary <- sample_buildings |>
  group_by(age_category) |>
  summarise(
    n = n(),
    mean_initial_psi = mean(initial_psi),
    .groups = "drop"
  )

cat("\nInitial damage by building age:\n")
print(age_summary)
cat("\n")

# -----------------------------------------------------------------------------
# 4. Calculate PGV Matrix
# -----------------------------------------------------------------------------

cat("Step 4: Calculating PGV at each building for each event...\n")

n_buildings <- nrow(sample_buildings)
n_events <- nrow(groningen_eq)

pgv_matrix <- matrix(NA, n_buildings, n_events)

pb <- txtProgressBar(min = 0, max = n_events, style = 3)
for (j in seq_len(n_events)) {
  eq <- groningen_eq[j, ]

  # Epicentral distance
  epi_dist <- calc_dist_km(
    sample_buildings$lat, sample_buildings$lon,
    eq$lat, eq$lon
  )

  # Hypocentral distance
  hypo_dist <- sqrt(epi_dist^2 + eq$depth_km^2)

  # PGV from Bommer et al. (2022)
  result <- pgv_model(
    M = eq$mag,
    Rhyp = hypo_dist,
    comp = "LRG",
    VS30 = sample_buildings$vs30
  )

  pgv_matrix[, j] <- result$PGV * 10  # Convert cm/s to mm/s

  setTxtProgressBar(pb, j)
}
close(pb)

cat("PGV matrix:", n_buildings, "buildings ×", n_events, "events\n")
cat("PGV range:", round(min(pgv_matrix), 2), "to",
    round(max(pgv_matrix), 2), "mm/s\n\n")

# -----------------------------------------------------------------------------
# 5. Load brms Model and Extract Posteriors for Fast Prediction
# -----------------------------------------------------------------------------

cat("Step 5: Loading brms model and extracting posteriors...\n")

fit_brms <- readRDS(here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here("outputs", "models", "standardization_params.rds"))

# Extract posterior draws
post <- as_draws_df(fit_brms)
n_post <- nrow(post)

# Fixed effects
b_intercept <- post$b_Intercept
b_initial_psi <- post$b_initial_psi_c
bsp_pgv <- post$bsp_mopgv_ord
bsp_n <- post$bsp_mon_ord
bsp_material <- post$bsp_momaterial_ord
shape <- post$shape

hu_intercept <- post$b_hu_Intercept
hu_initial_psi <- post$b_hu_initial_psi_c
hu_bsp_pgv <- post$bsp_hu_mopgv_ord
hu_bsp_n <- post$bsp_hu_mon_ord
hu_bsp_material <- post$bsp_hu_momaterial_ord

# Random effect SDs
sd_eq <- post$`sd_EarthquakeType__Intercept`
sd_fs <- post$`sd_FacadeType:SoilProfile__Intercept`
hu_sd_eq <- post$`sd_EarthquakeType__hu_Intercept`
hu_sd_fs <- post$`sd_FacadeType:SoilProfile__hu_Intercept`

# Monotonic simplex weights (scaled by K)
pgv_K <- 7
pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) {
  pgv_simplex[, i] <- post[[paste0("simo_mopgv_ord1[", i, "]")]]
}
pgv_cumsum <- pgv_K * t(apply(pgv_simplex, 1, cumsum))

hu_pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) {
  hu_pgv_simplex[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
}
hu_pgv_cumsum <- pgv_K * t(apply(hu_pgv_simplex, 1, cumsum))

# Material simplex (K=2 for 3 levels)
mat_K <- 2
mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) {
  mat_simplex[, i] <- post[[paste0("simo_momaterial_ord1[", i, "]")]]
}
mat_cumsum <- mat_K * t(apply(mat_simplex, 1, cumsum))

hu_mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) {
  hu_mat_simplex[, i] <- post[[paste0("simo_hu_momaterial_ord1[", i, "]")]]
}
hu_mat_cumsum <- mat_K * t(apply(hu_mat_simplex, 1, cumsum))

cat("Extracted", n_post, "posterior draws\n\n")

# -----------------------------------------------------------------------------
# 6. Fast Prediction Function
# -----------------------------------------------------------------------------

# Map PGV to level index (1-8)
pgv_to_idx <- function(pgv) {
  pgv_levels <- c(2, 4, 8, 16, 32, 64, 96, 128)
  sapply(pgv, function(p) {
    idx <- which.min(abs(pgv_levels - p))
    return(idx)
  })
}

# Fast damage prediction (single building, single draw)
fast_predict_single <- function(pgv, initial_psi, material_idx = 2) {
  # Sample single posterior index

idx <- sample(n_post, 1)

  # PGV level
  pgv_idx <- pgv_to_idx(pgv)
  if (pgv_idx == 1) {
    pgv_effect <- 0
    hu_pgv_effect <- 0
  } else {
    pgv_effect <- pgv_cumsum[idx, pgv_idx - 1]
    hu_pgv_effect <- hu_pgv_cumsum[idx, pgv_idx - 1]
  }

  # Material effect
  if (material_idx == 1) {
    mat_effect <- 0
    hu_mat_effect <- 0
  } else {
    mat_effect <- mat_cumsum[idx, material_idx - 1]
    hu_mat_effect <- hu_mat_cumsum[idx, material_idx - 1]
  }

  # Random effects (marginal)
  r_eq <- rnorm(1, 0, sd_eq[idx])
  r_fs <- rnorm(1, 0, sd_fs[idx])
  hu_r_eq <- rnorm(1, 0, hu_sd_eq[idx])
  hu_r_fs <- rnorm(1, 0, hu_sd_fs[idx])

  # Center initial_psi
  psi_c <- initial_psi - std_params$initial_psi_mean

  # Linear predictors
  eta_mu <- b_intercept[idx] +
            bsp_pgv[idx] * pgv_effect +
            bsp_material[idx] * mat_effect +
            b_initial_psi[idx] * psi_c +
            r_eq + r_fs

  eta_hu <- hu_intercept[idx] +
            hu_bsp_pgv[idx] * hu_pgv_effect +
            hu_bsp_material[idx] * hu_mat_effect +
            hu_initial_psi[idx] * psi_c +
            hu_r_eq + hu_r_fs

  # Transform
  mu <- exp(eta_mu)
  p_zero <- plogis(eta_hu)

  # Sample from hurdle-gamma
  if (runif(1) < p_zero) {
    return(0)
  } else {
    return(rgamma(1, shape = shape[idx], rate = shape[idx] / mu))
  }
}

# -----------------------------------------------------------------------------
# 7. Cumulative Damage Simulation (Both Scenarios)
# -----------------------------------------------------------------------------

cat("Step 6: Running cumulative damage simulation (both scenarios)...\n")
cat("Processing", n_buildings, "buildings through", n_events, "events...\n")

# Sort events chronologically
event_order <- order(groningen_eq$time_utc)

# Damage threshold for counting events
pgv_threshold <- 2  # mm/s - below this, assume no damage

# Function to run simulation with given initial state
run_damage_simulation <- function(initial_psi_vec, scenario_name) {
  cat("\n  Running", scenario_name, "scenario...\n")

  # Initialize damage state
  psi <- initial_psi_vec
  n_damaging_events <- rep(0, length(psi))

  pb <- txtProgressBar(min = 0, max = n_events, style = 3)
  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]
    pgv_event <- pgv_matrix[, j]

    # Only process buildings with PGV above threshold
    affected <- which(pgv_event > pgv_threshold)

    if (length(affected) > 0) {
      for (i in affected) {
        delta <- fast_predict_single(
          pgv = pgv_event[i],
          initial_psi = psi[i],
          material_idx = 2  # Average material
        )
        psi[i] <- psi[i] + delta
        if (delta > 0) {
          n_damaging_events[i] <- n_damaging_events[i] + 1
        }
      }
    }

    setTxtProgressBar(pb, j_idx)
  }
  close(pb)

  list(psi = psi, n_damaging_events = n_damaging_events)
}

# Scenario A: Virgin walls (Ψ₀ = 0)
set.seed(42)  # For reproducibility
results_virgin <- run_damage_simulation(
  initial_psi_vec = rep(0, n_buildings),
  scenario_name = "Virgin Walls"
)

# Scenario B: Age-based pre-damage
set.seed(42)  # Same seed for fair comparison
results_predamage <- run_damage_simulation(
  initial_psi_vec = sample_buildings$initial_psi,
  scenario_name = "Age-Based Pre-Damage"
)

# Store results
sample_buildings$psi_virgin <- results_virgin$psi
sample_buildings$psi_predamage <- results_predamage$psi
sample_buildings$n_damaging_events <- results_virgin$n_damaging_events
sample_buildings$delta_psi_virgin <- results_virgin$psi  # Same as psi_virgin since start from 0
sample_buildings$delta_psi_predamage <- results_predamage$psi - sample_buildings$initial_psi

cat("\n\n")

# -----------------------------------------------------------------------------
# 8. Summary Statistics
# -----------------------------------------------------------------------------

cat("Step 7: Computing summary statistics...\n")

# Per-building statistics
sample_buildings$max_pgv <- apply(pgv_matrix, 1, max)
sample_buildings$sum_pgv <- apply(pgv_matrix, 1, sum)

# Damage thresholds for both scenarios
sample_buildings$visible_virgin <- as.integer(sample_buildings$psi_virgin >= 1)
sample_buildings$visible_predamage <- as.integer(sample_buildings$psi_predamage >= 1)
sample_buildings$moderate_virgin <- as.integer(sample_buildings$psi_virgin >= 2)
sample_buildings$moderate_predamage <- as.integer(sample_buildings$psi_predamage >= 2)

# Overall summary
cat("\n=== Results Summary ===\n")
cat("Buildings analyzed:", nrow(sample_buildings), "\n")
cat("Earthquakes processed:", n_events, "\n\n")

cat("--- SCENARIO A: Virgin Walls (Ψ₀ = 0) ---\n")
cat("  Mean cumulative Ψ:", round(mean(sample_buildings$psi_virgin), 3), "\n")
cat("  Median cumulative Ψ:", round(median(sample_buildings$psi_virgin), 3), "\n")
cat("  Max cumulative Ψ:", round(max(sample_buildings$psi_virgin), 3), "\n")
cat("  % with visible damage (Ψ≥1):",
    round(100 * mean(sample_buildings$visible_virgin), 1), "%\n")
cat("  % with moderate damage (Ψ≥2):",
    round(100 * mean(sample_buildings$moderate_virgin), 1), "%\n\n")

cat("--- SCENARIO B: Age-Based Pre-Damage ---\n")
cat("  Mean initial Ψ₀:", round(mean(sample_buildings$initial_psi), 3), "\n")
cat("  Mean final Ψ:", round(mean(sample_buildings$psi_predamage), 3), "\n")
cat("  Mean ΔΨ (earthquake-only):", round(mean(sample_buildings$delta_psi_predamage), 3), "\n")
cat("  % with visible damage (Ψ≥1):",
    round(100 * mean(sample_buildings$visible_predamage), 1), "%\n")
cat("  % with moderate damage (Ψ≥2):",
    round(100 * mean(sample_buildings$moderate_predamage), 1), "%\n\n")

cat("--- COMPARISON ---\n")
cat("  Additional visible damage from pre-damage:",
    round(100 * (mean(sample_buildings$visible_predamage) - mean(sample_buildings$visible_virgin)), 1),
    "percentage points\n")
cat("  ΔΨ reduction due to saturation:",
    round(100 * (1 - mean(sample_buildings$delta_psi_predamage) / mean(sample_buildings$delta_psi_virgin)), 1),
    "% less earthquake damage\n\n")

cat("PGV statistics:\n")
cat("  Max PGV experienced:", round(max(sample_buildings$max_pgv), 1), "mm/s\n")
cat("  Mean max PGV:", round(mean(sample_buildings$max_pgv), 1), "mm/s\n\n")

# By distance zone - Virgin walls
by_zone_virgin <- sample_buildings |>
  mutate(zone = case_when(
    dist_to_center < 10 ~ "Near-field (<10km)",
    dist_to_center < 25 ~ "Mid-field (10-25km)",
    TRUE ~ "Far-field (>25km)"
  )) |>
  group_by(zone) |>
  summarise(
    n = n(),
    mean_psi = mean(psi_virgin),
    p_visible = mean(visible_virgin),
    p_moderate = mean(moderate_virgin),
    mean_max_pgv = mean(max_pgv),
    .groups = "drop"
  )

cat("Results by distance zone (Virgin Walls):\n")
print(by_zone_virgin)
cat("\n")

# By distance zone - Pre-damage
by_zone_predamage <- sample_buildings |>
  mutate(zone = case_when(
    dist_to_center < 10 ~ "Near-field (<10km)",
    dist_to_center < 25 ~ "Mid-field (10-25km)",
    TRUE ~ "Far-field (>25km)"
  )) |>
  group_by(zone) |>
  summarise(
    n = n(),
    mean_psi = mean(psi_predamage),
    mean_initial = mean(initial_psi),
    p_visible = mean(visible_predamage),
    p_moderate = mean(moderate_predamage),
    .groups = "drop"
  )

cat("Results by distance zone (Pre-Damage):\n")
print(by_zone_predamage)

# By building age
by_age <- sample_buildings |>
  group_by(age_category) |>
  summarise(
    n = n(),
    mean_initial = mean(initial_psi),
    mean_final = mean(psi_predamage),
    p_visible_virgin = mean(visible_virgin),
    p_visible_predamage = mean(visible_predamage),
    increase = mean(visible_predamage) - mean(visible_virgin),
    .groups = "drop"
  )

cat("\nResults by building age:\n")
print(by_age)

# -----------------------------------------------------------------------------
# 9. Visualizations
# -----------------------------------------------------------------------------

cat("\n\nStep 8: Generating visualizations...\n")

# Combined histogram: Virgin vs Pre-damage
plot_data_hist <- sample_buildings |>
  select(psi_virgin, psi_predamage) |>
  pivot_longer(everything(), names_to = "scenario", values_to = "psi") |>
  mutate(scenario = case_when(
    scenario == "psi_virgin" ~ "Virgin Walls",
    scenario == "psi_predamage" ~ "Age-Based Pre-Damage"
  ))

p_hist <- ggplot(plot_data_hist, aes(x = psi, fill = scenario)) +
  geom_histogram(binwidth = 0.2, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1.1, y = Inf, label = "Visible damage",
           color = "red", hjust = 0, vjust = 2) +
  scale_fill_manual(values = c("Virgin Walls" = "steelblue",
                               "Age-Based Pre-Damage" = "darkorange")) +
  labs(
    title = "Distribution of Cumulative Damage: Virgin vs Pre-Damage",
    subtitle = paste(n_events, "earthquakes (M >= 2.0), 1991-2024"),
    x = expression(Cumulative~Psi),
    y = "Number of buildings",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "spatial_damage_histogram.png"), p_hist,
       width = 10, height = 6, dpi = 150)
cat("Saved: spatial_damage_histogram.png\n")

# Damage vs distance - both scenarios
plot_data_dist <- sample_buildings |>
  select(dist_to_center, psi_virgin, psi_predamage) |>
  pivot_longer(c(psi_virgin, psi_predamage), names_to = "scenario", values_to = "psi") |>
  mutate(scenario = case_when(
    scenario == "psi_virgin" ~ "Virgin Walls",
    scenario == "psi_predamage" ~ "Pre-Damage"
  ))

p_dist <- ggplot(plot_data_dist, aes(x = dist_to_center, y = psi, color = scenario)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("Virgin Walls" = "steelblue",
                                "Pre-Damage" = "darkorange")) +
  labs(
    title = "Cumulative Damage vs Distance: Virgin vs Pre-Damage",
    x = "Distance to Loppersum (km)",
    y = expression(Cumulative~Psi),
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "spatial_damage_vs_distance.png"), p_dist,
       width = 10, height = 6, dpi = 150)
cat("Saved: spatial_damage_vs_distance.png\n")

# Damage vs max PGV - virgin walls
p_pgv <- sample_buildings |>
  ggplot(aes(x = max_pgv, y = psi_virgin)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "loess", color = "red") +
  geom_vline(xintercept = 15.3, linetype = "dashed", color = "blue") +
  annotate("text", x = 16, y = max(sample_buildings$psi_virgin) * 0.9,
           label = "10% threshold\n(15.3 mm/s)", color = "blue", hjust = 0, size = 3) +
  labs(
    title = "Cumulative Damage vs Maximum PGV (Virgin Walls)",
    x = "Maximum PGV (mm/s)",
    y = expression(Cumulative~Psi)
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "spatial_damage_vs_pgv.png"), p_pgv,
       width = 8, height = 5, dpi = 150)
cat("Saved: spatial_damage_vs_pgv.png\n")

# P(visible damage) by scenario and zone
p_compare_zones <- sample_buildings |>
  mutate(zone = case_when(
    dist_to_center < 10 ~ "Near (<10km)",
    dist_to_center < 25 ~ "Mid (10-25km)",
    TRUE ~ "Far (>25km)"
  )) |>
  group_by(zone) |>
  summarise(
    `Virgin Walls` = mean(visible_virgin),
    `Pre-Damage` = mean(visible_predamage),
    .groups = "drop"
  ) |>
  pivot_longer(-zone, names_to = "scenario", values_to = "p_visible") |>
  ggplot(aes(x = zone, y = p_visible, fill = scenario)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = scales::percent(p_visible, accuracy = 0.1)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = c("Virgin Walls" = "steelblue",
                               "Pre-Damage" = "darkorange")) +
  labs(
    title = "P(Visible Damage) by Zone: Virgin vs Pre-Damage",
    subtitle = "Effect of age-based initial damage assumption",
    x = "Distance Zone",
    y = "P(Ψ ≥ 1)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "spatial_fragility_by_zone.png"), p_compare_zones,
       width = 10, height = 6, dpi = 150)
cat("Saved: spatial_fragility_by_zone.png\n")

# Damage by building age
p_age <- sample_buildings |>
  group_by(age_category) |>
  summarise(
    `Virgin Walls` = mean(visible_virgin),
    `Pre-Damage` = mean(visible_predamage),
    `Initial (Ψ₀)` = mean(initial_psi > 0),  # Fraction with pre-damage
    .groups = "drop"
  ) |>
  mutate(age_category = factor(age_category,
    levels = c("Pre-1920", "1920-1949", "1950-1969", "1970-1990", "Post-1991"))) |>
  pivot_longer(-age_category, names_to = "metric", values_to = "value") |>
  ggplot(aes(x = age_category, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Virgin Walls" = "steelblue",
                               "Pre-Damage" = "darkorange",
                               "Initial (Ψ₀)" = "gray50")) +
  labs(
    title = "P(Visible Damage) by Building Age",
    subtitle = "Older buildings assumed to have higher initial damage",
    x = "Construction Period",
    y = "Proportion",
    fill = "Metric"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "spatial_damage_by_age.png"), p_age,
       width = 10, height = 6, dpi = 150)
cat("Saved: spatial_damage_by_age.png\n")

# -----------------------------------------------------------------------------
# 10. Save Results
# -----------------------------------------------------------------------------

cat("\nStep 9: Saving results...\n")

results <- sample_buildings |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, age_category,
    initial_psi,
    psi_virgin, psi_predamage,
    delta_psi_virgin, delta_psi_predamage,
    max_pgv, sum_pgv, n_damaging_events,
    visible_virgin, visible_predamage,
    moderate_virgin, moderate_predamage
  )

saveRDS(results, here("outputs", "models", "spatial_damage_results.rds"))
cat("Saved: outputs/models/spatial_damage_results.rds\n")

# Save earthquake catalogue used
saveRDS(groningen_eq, here("outputs", "models", "groningen_earthquakes.rds"))
cat("Saved: outputs/models/groningen_earthquakes.rds\n")

# Save summary tables
summary_tables <- list(
  by_zone_virgin = by_zone_virgin,
  by_zone_predamage = by_zone_predamage,
  by_age = by_age,
  age_initial_damage = age_summary
)
saveRDS(summary_tables, here("outputs", "models", "spatial_damage_summaries.rds"))
cat("Saved: outputs/models/spatial_damage_summaries.rds\n")

cat("\n=== Spatial Damage Prediction Complete ===\n")
cat("Scenarios: Virgin Walls vs Age-Based Pre-Damage\n")
cat("Magnitude threshold: M >= 2.0\n")
