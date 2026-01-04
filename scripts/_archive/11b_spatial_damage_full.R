# =============================================================================
# 11b_spatial_damage_full.R
#
# Full-scale spatial damage prediction for ALL Groningen buildings
# Optimized for ~543k buildings with batched processing
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

# =============================================================================
# CONFIGURATION
# =============================================================================

PSI_MAX <- 4.0          # Physical damage cap (FEM max ~3.6)
MAG_THRESHOLD <- 2.0    # Minimum magnitude
PGV_THRESHOLD <- 8.0    # mm/s - model well-calibrated for PGV≥8 (overpredicts 1.8x at PGV<8)
BATCH_SIZE <- 10000     # Buildings per batch for memory management

cat("=== Full-Scale Groningen Spatial Damage Prediction ===\n")
cat("Configuration:\n")
cat("  PSI_MAX:", PSI_MAX, "\n")
cat("  MAG_THRESHOLD:", MAG_THRESHOLD, "\n")
cat("  BATCH_SIZE:", BATCH_SIZE, "\n\n")

# =============================================================================
# 1. LOAD EARTHQUAKE CATALOGUE
# =============================================================================

cat("Step 1: Loading earthquake catalogue...\n")

# Try to load cached catalogue first
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
}

n_events <- nrow(groningen_eq)
cat("Earthquakes (M >=", MAG_THRESHOLD, "):", n_events, "\n")
cat("Date range:", as.character(min(groningen_eq$time_utc)), "to",
    as.character(max(groningen_eq$time_utc)), "\n\n")

# =============================================================================
# 2. LOAD BUILDINGS
# =============================================================================

cat("Step 2: Loading building data...\n")

buildings <- st_read(here("datafiles", "bag_selectie.gpkg"), quiet = TRUE)
n_buildings <- nrow(buildings)
cat("Total buildings:", format(n_buildings, big.mark = ","), "\n")

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
cat("Extracting centroids...\n")
buildings_wgs84 <- st_transform(buildings, 4326)
centroids <- st_centroid(st_geometry(buildings_wgs84))
coords <- st_coordinates(centroids)
buildings$lon <- coords[, 1]
buildings$lat <- coords[, 2]

# Distance helper
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

# Age-based initial damage
assign_initial_damage <- function(bouwjaar) {
  case_when(
    bouwjaar < 1920 ~ 0.8,
    bouwjaar < 1950 ~ 0.5,
    bouwjaar < 1970 ~ 0.3,
    bouwjaar < 1991 ~ 0.15,
    TRUE ~ 0
  )
}

buildings <- buildings |>
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

# Drop geometry for processing
buildings_df <- st_drop_geometry(buildings)
cat("Buildings prepared for processing\n\n")

# =============================================================================
# 2b. PRE-GENERATE GMM SITE TERMS (φS2S)
# =============================================================================

cat("Step 2b: Generating site-specific GMM uncertainty terms...\n")

# Get GMM uncertainty components from Bommer et al. (2022)
gmm_coef <- get_pgv_coefficients("LRG")
gmm_tau <- gmm_coef$tau        # Inter-event (τ = 0.2448)
gmm_phiS2S <- gmm_coef$phiS2S  # Site-to-site (φS2S = 0.2406)
gmm_phiSS <- gmm_coef$phiSS    # Within-event single-station (φSS = 0.4569)

cat("  GMM uncertainty components (ln units):\n")
cat("    τ (inter-event):", round(gmm_tau, 4), "\n")
cat("    φS2S (site-to-site):", round(gmm_phiS2S, 4), "\n")
cat("    φSS (within-event):", round(gmm_phiSS, 4), "\n")
cat("    σ_total:", round(sqrt(gmm_tau^2 + gmm_phiS2S^2 + gmm_phiSS^2), 4), "\n")

# Pre-generate site terms - these persist across all earthquakes for each building
set.seed(123)  # Reproducible site terms
buildings_df$eta_site <- rnorm(n_buildings, 0, gmm_phiS2S)
cat("  Generated", format(n_buildings, big.mark = ","), "site-specific η_S2S terms\n\n")

# =============================================================================
# 3. LOAD BRMS MODEL
# =============================================================================

cat("Step 3: Loading brms model and extracting posteriors...\n")

fit_brms <- readRDS(here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here("outputs", "models", "standardization_params.rds"))

post <- as_draws_df(fit_brms)
n_post <- nrow(post)

# Extract all parameters once
b_intercept <- post$b_Intercept
b_initial_psi <- post$b_initial_psi_c
bsp_pgv <- post$bsp_mopgv_ord
bsp_material <- post$bsp_momaterial_ord
shape <- post$shape

hu_intercept <- post$b_hu_Intercept
hu_initial_psi <- post$b_hu_initial_psi_c
hu_bsp_pgv <- post$bsp_hu_mopgv_ord
hu_bsp_material <- post$bsp_hu_momaterial_ord

sd_eq <- post$`sd_EarthquakeType__Intercept`
sd_fs <- post$`sd_FacadeType:SoilProfile__Intercept`
hu_sd_eq <- post$`sd_EarthquakeType__hu_Intercept`
hu_sd_fs <- post$`sd_FacadeType:SoilProfile__hu_Intercept`

# Monotonic simplex (PGV: K=7)
pgv_K <- 7
pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) pgv_simplex[, i] <- post[[paste0("simo_mopgv_ord1[", i, "]")]]
pgv_cumsum <- pgv_K * t(apply(pgv_simplex, 1, cumsum))

hu_pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) hu_pgv_simplex[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
hu_pgv_cumsum <- pgv_K * t(apply(hu_pgv_simplex, 1, cumsum))

# Material simplex (K=2)
mat_K <- 2
mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) mat_simplex[, i] <- post[[paste0("simo_momaterial_ord1[", i, "]")]]
mat_cumsum <- mat_K * t(apply(mat_simplex, 1, cumsum))

hu_mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) hu_mat_simplex[, i] <- post[[paste0("simo_hu_momaterial_ord1[", i, "]")]]
hu_mat_cumsum <- mat_K * t(apply(hu_mat_simplex, 1, cumsum))

cat("Extracted", n_post, "posterior draws\n\n")

# =============================================================================
# 4. FAST PREDICTION FUNCTIONS
# =============================================================================

pgv_levels <- c(2, 4, 8, 16, 32, 64, 96, 128)

# Vectorized PGV to index
pgv_to_idx <- function(pgv) {
  idx <- findInterval(pgv, c(-Inf, 3, 6, 12, 24, 48, 80, 112, Inf))
  pmax(1, pmin(idx, 8))
}

# Vectorized damage prediction for a batch
predict_damage_batch <- function(pgv_vec, psi_vec, n_samples = 1) {
  n <- length(pgv_vec)
  delta <- numeric(n)

  # Sample posterior indices
  idx <- sample(n_post, n, replace = TRUE)

  # PGV levels
  pgv_idx <- pgv_to_idx(pgv_vec)

  # Center initial psi
  psi_c <- psi_vec - std_params$initial_psi_mean

  for (i in seq_len(n)) {
    if (pgv_vec[i] < PGV_THRESHOLD) {
      delta[i] <- 0
      next
    }

    j <- idx[i]

    # PGV effect
    if (pgv_idx[i] == 1) {
      pgv_eff <- 0; hu_pgv_eff <- 0
    } else {
      pgv_eff <- pgv_cumsum[j, pgv_idx[i] - 1]
      hu_pgv_eff <- hu_pgv_cumsum[j, pgv_idx[i] - 1]
    }

    # Material effect (average = level 2)
    mat_eff <- mat_cumsum[j, 1]
    hu_mat_eff <- hu_mat_cumsum[j, 1]

    # Random effects
    r_eq <- rnorm(1, 0, sd_eq[j])
    r_fs <- rnorm(1, 0, sd_fs[j])
    hu_r_eq <- rnorm(1, 0, hu_sd_eq[j])
    hu_r_fs <- rnorm(1, 0, hu_sd_fs[j])

    # Linear predictors
    eta_mu <- b_intercept[j] + bsp_pgv[j] * pgv_eff + bsp_material[j] * mat_eff +
              b_initial_psi[j] * psi_c[i] + r_eq + r_fs
    eta_hu <- hu_intercept[j] + hu_bsp_pgv[j] * hu_pgv_eff + hu_bsp_material[j] * hu_mat_eff +
              hu_initial_psi[j] * psi_c[i] + hu_r_eq + hu_r_fs

    # Sample from hurdle-gamma
    if (runif(1) < plogis(eta_hu)) {
      delta[i] <- 0
    } else {
      delta[i] <- rgamma(1, shape = shape[j], rate = shape[j] / exp(eta_mu))
    }
  }

  return(delta)
}

# =============================================================================
# 5. MAIN SIMULATION LOOP
# =============================================================================

cat("Step 4: Running cumulative damage simulation...\n")
cat("Processing", format(n_buildings, big.mark = ","), "buildings ×", n_events, "events\n")
cat("Estimated batches:", ceiling(n_buildings / BATCH_SIZE), "\n\n")

# Initialize results
buildings_df$psi_virgin <- 0
buildings_df$psi_predamage <- buildings_df$initial_psi
buildings_df$max_pgv <- 0
buildings_df$n_events_felt <- 0

# Sort events chronologically
event_order <- order(groningen_eq$time_utc)

# Pre-compute earthquake data
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

set.seed(42)
start_time <- Sys.time()

for (j_idx in seq_along(event_order)) {
  j <- event_order[j_idx]

  # Calculate MEDIAN PGV for all buildings (vectorized)
  epi_dist <- calc_dist_km(buildings_df$lat, buildings_df$lon, eq_lat[j], eq_lon[j])
  hypo_dist <- sqrt(epi_dist^2 + eq_depth[j]^2)

  pgv_result <- pgv_model(M = eq_mag[j], Rhyp = hypo_dist, comp = "LRG", VS30 = buildings_df$vs30)
  lnPGV_median <- pgv_result$lnPGV  # Natural log of median PGV (cm/s)

  # === PROPER UNCERTAINTY CONVOLUTION ===
  # 1. Event term (τ) - same for all buildings this earthquake
  eta_event <- rnorm(1, 0, gmm_tau)

  # 2. Site term (φS2S) - pre-generated, persists across earthquakes
  # Already in buildings_df$eta_site

  # 3. Within-event term (φSS) - different for each building-event pair
  eta_within <- rnorm(n_buildings, 0, gmm_phiSS)

  # Realized lnPGV = median + event + site + within-event
  lnPGV_realized <- lnPGV_median + eta_event + buildings_df$eta_site + eta_within

  # Convert to mm/s (original is cm/s)
  pgv_event <- exp(lnPGV_realized) * 10  # cm/s -> mm/s

  # Update max PGV (using realized values)
  buildings_df$max_pgv <- pmax(buildings_df$max_pgv, pgv_event)

  # Find affected buildings
  affected <- which(pgv_event >= PGV_THRESHOLD)
  n_affected <- length(affected)

  if (n_affected > 0) {
    buildings_df$n_events_felt[affected] <- buildings_df$n_events_felt[affected] + 1

    # Process in batches
    n_batches <- ceiling(n_affected / BATCH_SIZE)

    for (b in seq_len(n_batches)) {
      batch_start <- (b - 1) * BATCH_SIZE + 1
      batch_end <- min(b * BATCH_SIZE, n_affected)
      batch_idx <- affected[batch_start:batch_end]

      # Virgin walls prediction
      delta_virgin <- predict_damage_batch(
        pgv_vec = pgv_event[batch_idx],
        psi_vec = buildings_df$psi_virgin[batch_idx]
      )
      buildings_df$psi_virgin[batch_idx] <- pmin(
        buildings_df$psi_virgin[batch_idx] + delta_virgin, PSI_MAX
      )

      # Pre-damage prediction
      delta_predamage <- predict_damage_batch(
        pgv_vec = pgv_event[batch_idx],
        psi_vec = buildings_df$psi_predamage[batch_idx]
      )
      buildings_df$psi_predamage[batch_idx] <- pmin(
        buildings_df$psi_predamage[batch_idx] + delta_predamage, PSI_MAX
      )
    }
  }

  # Progress
  if (j_idx %% 10 == 0 || j_idx == n_events) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    eta <- elapsed / j_idx * (n_events - j_idx)
    cat(sprintf("\r  Event %d/%d (%.1f min elapsed, %.1f min remaining)    ",
                j_idx, n_events, elapsed, eta))
  }
}

cat("\n\nSimulation complete!\n\n")

# =============================================================================
# 6. COMPUTE SUMMARY STATISTICS
# =============================================================================

cat("Step 5: Computing statistics...\n")

buildings_df <- buildings_df |>
  mutate(
    delta_psi_virgin = psi_virgin,
    delta_psi_predamage = psi_predamage - initial_psi,
    visible_virgin = as.integer(psi_virgin >= 1),
    visible_predamage = as.integer(psi_predamage >= 1),
    moderate_virgin = as.integer(psi_virgin >= 2),
    moderate_predamage = as.integer(psi_predamage >= 2),
    severe_virgin = as.integer(psi_virgin >= 3),
    severe_predamage = as.integer(psi_predamage >= 3),
    zone = case_when(
      dist_to_center < 10 ~ "Near (<10km)",
      dist_to_center < 25 ~ "Mid (10-25km)",
      TRUE ~ "Far (>25km)"
    )
  )

# Summary output
cat("\n=== RESULTS SUMMARY ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("Earthquakes:", n_events, "(M >=", MAG_THRESHOLD, ")\n")
cat("Damage cap (Ψ_max):", PSI_MAX, "\n\n")

cat("--- VIRGIN WALLS ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_virgin), 3), "\n")
cat("% visible (Ψ≥1):", round(100 * mean(buildings_df$visible_virgin), 2), "%\n")
cat("% moderate (Ψ≥2):", round(100 * mean(buildings_df$moderate_virgin), 2), "%\n")
cat("% severe (Ψ≥3):", round(100 * mean(buildings_df$severe_virgin), 2), "%\n")
cat("% at cap (Ψ=4):", round(100 * mean(buildings_df$psi_virgin >= PSI_MAX), 2), "%\n\n")

cat("--- PRE-DAMAGE ---\n")
cat("Mean initial Ψ₀:", round(mean(buildings_df$initial_psi), 3), "\n")
cat("Mean final Ψ:", round(mean(buildings_df$psi_predamage), 3), "\n")
cat("% visible (Ψ≥1):", round(100 * mean(buildings_df$visible_predamage), 2), "%\n")
cat("% moderate (Ψ≥2):", round(100 * mean(buildings_df$moderate_predamage), 2), "%\n")
cat("% severe (Ψ≥3):", round(100 * mean(buildings_df$severe_predamage), 2), "%\n\n")

# By zone
by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    mean_psi_virgin = mean(psi_virgin),
    mean_psi_predamage = mean(psi_predamage),
    p_visible_virgin = mean(visible_virgin),
    p_visible_predamage = mean(visible_predamage),
    mean_max_pgv = mean(max_pgv),
    .groups = "drop"
  )

cat("By distance zone:\n")
print(by_zone)

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

cat("\nStep 6: Saving results...\n")

results <- buildings_df |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, zone, age_category,
    initial_psi, psi_virgin, psi_predamage,
    delta_psi_virgin, delta_psi_predamage,
    max_pgv, n_events_felt,
    visible_virgin, visible_predamage,
    moderate_virgin, moderate_predamage,
    severe_virgin, severe_predamage
  )

saveRDS(results, here("outputs", "models", "spatial_damage_full.rds"))
cat("Saved: outputs/models/spatial_damage_full.rds\n")

# Save summary
summary_tables <- list(
  by_zone = by_zone,
  by_age = buildings_df |>
    group_by(age_category) |>
    summarise(
      n = n(),
      mean_initial = mean(initial_psi),
      mean_final_virgin = mean(psi_virgin),
      mean_final_predamage = mean(psi_predamage),
      p_visible_virgin = mean(visible_virgin),
      p_visible_predamage = mean(visible_predamage),
      .groups = "drop"
    ),
  by_postcode = buildings_df |>
    group_by(pc4) |>
    summarise(
      n = n(),
      mean_psi_virgin = mean(psi_virgin),
      mean_psi_predamage = mean(psi_predamage),
      p_visible_virgin = mean(visible_virgin),
      p_visible_predamage = mean(visible_predamage),
      mean_max_pgv = mean(max_pgv),
      mean_lat = mean(lat),
      mean_lon = mean(lon),
      .groups = "drop"
    )
)

saveRDS(summary_tables, here("outputs", "models", "spatial_damage_full_summaries.rds"))
cat("Saved: outputs/models/spatial_damage_full_summaries.rds\n")

elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n=== Complete! Total time:", round(elapsed_total, 1), "minutes ===\n")
