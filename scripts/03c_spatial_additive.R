# =============================================================================
# 03c_spatial_additive.R
#
# Full uncertainty propagation for Groningen building damage predictions
# Vectorized across buildings AND posterior draws
#
# UNCERTAINTY SOURCES (fully integrated):
# 1. GMM uncertainty: τ (inter-event), φS2S (site-to-site), φSS (within-event)
#    → L independent GMM realizations
# 2. Fragility posterior uncertainty: M draws from brms posterior
#    → Different parameter sets for hurdle-gamma model
# 3. Stochastic sampling: hurdle + gamma draws
#    → Inherent randomness in hurdle-gamma predictions
#
# TOTAL SAMPLES: L × M = e.g., 10 × 50 = 500 per building
#
# RANDOM EFFECTS NOTE:
# The brms model has group-level random effects for EarthquakeType (4 levels)
# and FacadeType:SoilProfile (4 levels). Since we don't have group membership
# for new buildings, we use MARGINAL prediction: sample fresh random effects
# from their estimated distributions for each prediction, effectively
# integrating over unknown group membership.
#
# MONOTONIC EFFECTS NOTE:
# brms mo() effects use cumulative simplex × K where K = #levels - 1:
#   - PGV (8 levels) → K = 7
#   - Material (3 levels) → K = 2
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(sf)
library(brms)
library(here)
library(matrixStats)

# Source helper functions
source(here("helperfuncties", "fetch_knmi_events.R"))
source(here("helperfuncties", "pgv_model_original.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

PSI_MAX <- 4.0          # Physical damage cap
MAG_THRESHOLD <- 2.0    # Minimum magnitude
PGV_THRESHOLD <- 8.0    # mm/s - model threshold
L_GMM <- 10             # Number of GMM realizations (outer loop)
M_DRAWS <- 50           # Number of posterior draws (inner loop)
SAMPLE_FRAC <- 1.0      # Full dataset (set to 0.1 for validation)

cat("=== Spatial Damage Prediction with Full Uncertainty ===\n")
cat("Configuration:\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  Posterior draws (M):", M_DRAWS, "\n")
cat("  Total samples per building:", L_GMM * M_DRAWS, "\n")
cat("  Sample fraction:", SAMPLE_FRAC, "\n")
cat("  PGV threshold:", PGV_THRESHOLD, "mm/s\n\n")

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
}

n_events <- nrow(groningen_eq)
cat("Earthquakes (M >=", MAG_THRESHOLD, "):", n_events, "\n\n")

# =============================================================================
# 2. LOAD BUILDINGS (with optional stratified sampling)
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

cat("  τ (inter-event):", round(gmm_tau, 4), "\n")
cat("  φS2S (site-to-site):", round(gmm_phiS2S, 4), "\n")
cat("  φSS (within-event):", round(gmm_phiSS, 4), "\n")

# Pre-generate site terms (φS2S)
set.seed(123)
eta_site <- rnorm(n_buildings, 0, gmm_phiS2S)
cat("  Generated site terms for", n_buildings, "buildings\n\n")

# =============================================================================
# 4. LOAD BRMS MODEL AND EXTRACT POSTERIORS
# =============================================================================

cat("Step 4: Loading brms model and sampling posteriors...\n")

fit_brms <- readRDS(here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here("outputs", "models", "standardization_params.rds"))

post <- as_draws_df(fit_brms)
n_post <- nrow(post)

# Sample M posterior draws (fixed for entire simulation)
set.seed(42)
posterior_idx <- sample(n_post, M_DRAWS, replace = FALSE)
cat("  Sampled", M_DRAWS, "posterior draws from", n_post, "total\n")

# Extract parameters for selected draws only (length M each)
b_intercept <- post$b_Intercept[posterior_idx]
b_initial_psi <- post$b_initial_psi_c[posterior_idx]
bsp_pgv <- post$bsp_mopgv_ord[posterior_idx]
bsp_material <- post$bsp_momaterial_ord[posterior_idx]
shape <- post$shape[posterior_idx]

hu_intercept <- post$b_hu_Intercept[posterior_idx]
hu_initial_psi <- post$b_hu_initial_psi_c[posterior_idx]
hu_bsp_pgv <- post$bsp_hu_mopgv_ord[posterior_idx]
hu_bsp_material <- post$bsp_hu_momaterial_ord[posterior_idx]

sd_eq <- post$`sd_EarthquakeType__Intercept`[posterior_idx]
sd_fs <- post$`sd_FacadeType:SoilProfile__Intercept`[posterior_idx]
hu_sd_eq <- post$`sd_EarthquakeType__hu_Intercept`[posterior_idx]
hu_sd_fs <- post$`sd_FacadeType:SoilProfile__hu_Intercept`[posterior_idx]

# Monotonic simplex for PGV (K=7)
pgv_K <- 7
pgv_simplex_full <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) pgv_simplex_full[, i] <- post[[paste0("simo_mopgv_ord1[", i, "]")]]
pgv_cumsum_full <- pgv_K * t(apply(pgv_simplex_full, 1, cumsum))
pgv_cumsum <- pgv_cumsum_full[posterior_idx, ]  # M × 7

hu_pgv_simplex_full <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) hu_pgv_simplex_full[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
hu_pgv_cumsum_full <- pgv_K * t(apply(hu_pgv_simplex_full, 1, cumsum))
hu_pgv_cumsum <- hu_pgv_cumsum_full[posterior_idx, ]  # M × 7

# Material simplex (K=2, we use level 2 = average)
mat_K <- 2
mat_simplex_full <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) mat_simplex_full[, i] <- post[[paste0("simo_momaterial_ord1[", i, "]")]]
mat_cumsum_full <- mat_K * t(apply(mat_simplex_full, 1, cumsum))
mat_eff <- mat_cumsum_full[posterior_idx, 1]  # M-vector (material level 2)

hu_mat_simplex_full <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) hu_mat_simplex_full[, i] <- post[[paste0("simo_hu_momaterial_ord1[", i, "]")]]
hu_mat_cumsum_full <- mat_K * t(apply(hu_mat_simplex_full, 1, cumsum))
hu_mat_eff <- hu_mat_cumsum_full[posterior_idx, 1]  # M-vector

cat("  Posterior parameters extracted\n\n")

# =============================================================================
# 5. VECTORIZED PREDICTION FUNCTION
# =============================================================================

pgv_to_idx <- function(pgv) {
  idx <- findInterval(pgv, c(-Inf, 3, 6, 12, 24, 48, 80, 112, Inf))
  pmax(1, pmin(idx, 8))
}

# Fully vectorized damage prediction: returns N × M matrix
predict_damage_matrix <- function(pgv_vec, psi_matrix) {
  N <- length(pgv_vec)
  M <- M_DRAWS

  # PGV levels for all buildings (length N)
  pgv_idx <- pgv_to_idx(pgv_vec)

  # Center psi (N × M matrix)
  psi_c <- psi_matrix - std_params$initial_psi_mean

  # PGV effect: for each building, get cumsum value for its pgv level
  # pgv_cumsum is M × 7, we need N × M
  # For buildings with pgv_idx == 1, effect is 0
  pgv_eff_mat <- matrix(0, N, M)
  hu_pgv_eff_mat <- matrix(0, N, M)
  for (k in 2:8) {
    mask <- pgv_idx == k
    if (any(mask)) {
      # pgv_cumsum[, k-1] is M-vector; replicate for each matching building
      pgv_eff_mat[mask, ] <- matrix(pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
      hu_pgv_eff_mat[mask, ] <- matrix(hu_pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
    }
  }

  # Material effect (same for all buildings, M-vector → N × M)
  mat_eff_mat <- matrix(mat_eff, N, M, byrow = TRUE)
  hu_mat_eff_mat <- matrix(hu_mat_eff, N, M, byrow = TRUE)

  # Random effects: sample N × M matrices
  r_eq <- matrix(rnorm(N * M, 0, rep(sd_eq, each = N)), N, M)
  r_fs <- matrix(rnorm(N * M, 0, rep(sd_fs, each = N)), N, M)
  hu_r_eq <- matrix(rnorm(N * M, 0, rep(hu_sd_eq, each = N)), N, M)
  hu_r_fs <- matrix(rnorm(N * M, 0, rep(hu_sd_fs, each = N)), N, M)

  # Intercepts and slopes: M-vectors → N × M matrices
  b_int_mat <- matrix(b_intercept, N, M, byrow = TRUE)
  bsp_pgv_mat <- matrix(bsp_pgv, N, M, byrow = TRUE)
  bsp_mat_mat <- matrix(bsp_material, N, M, byrow = TRUE)
  b_psi_mat <- matrix(b_initial_psi, N, M, byrow = TRUE)

  hu_int_mat <- matrix(hu_intercept, N, M, byrow = TRUE)
  hu_bsp_pgv_mat <- matrix(hu_bsp_pgv, N, M, byrow = TRUE)
  hu_bsp_mat_mat <- matrix(hu_bsp_material, N, M, byrow = TRUE)
  hu_b_psi_mat <- matrix(hu_initial_psi, N, M, byrow = TRUE)

  # Linear predictors (N × M)
  eta_mu <- b_int_mat + bsp_pgv_mat * pgv_eff_mat + bsp_mat_mat * mat_eff_mat +
            b_psi_mat * psi_c + r_eq + r_fs

  eta_hu <- hu_int_mat + hu_bsp_pgv_mat * hu_pgv_eff_mat + hu_bsp_mat_mat * hu_mat_eff_mat +
            hu_b_psi_mat * psi_c + hu_r_eq + hu_r_fs

  # Hurdle probability (N × M)
  hurdle_prob <- 1 / (1 + exp(-eta_hu))  # plogis
  is_zero <- matrix(runif(N * M), N, M) < hurdle_prob

  # Gamma mean and shape (N × M)
  gamma_mean <- exp(eta_mu)
  shape_mat <- matrix(shape, N, M, byrow = TRUE)

  # Sample from gamma (N × M)
  delta <- matrix(rgamma(N * M, shape = shape_mat, rate = shape_mat / gamma_mean), N, M)

  # Apply hurdle zeros
  delta[is_zero] <- 0

  # Apply PGV threshold (buildings below threshold get 0)
  below_threshold <- pgv_vec < PGV_THRESHOLD
  delta[below_threshold, ] <- 0

  return(delta)
}

# =============================================================================
# 6. MAIN SIMULATION LOOP
# =============================================================================

cat("Step 5: Running simulation with full uncertainty propagation...\n")
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("  Events:", n_events, "\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  Posterior draws (M):", M_DRAWS, "\n")
cat("  Total samples per building:", L_GMM * M_DRAWS, "\n")
cat("  Memory for results:",
    round(n_buildings * L_GMM * M_DRAWS * 8 * 2 / 1e6, 1), "MB\n\n")

# Sort events chronologically
event_order <- order(groningen_eq$time_utc)
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

# Pre-compute median lnPGV for all events (N × J matrix)
cat("  Pre-computing median lnPGV for all events...\n")
lnPGV_median_all <- matrix(NA, n_buildings, n_events)
for (j in 1:n_events) {
  epi_dist <- calc_dist_km(buildings_df$lat, buildings_df$lon, eq_lat[j], eq_lon[j])
  hypo_dist <- sqrt(epi_dist^2 + eq_depth[j]^2)
  pgv_result <- pgv_model(M = eq_mag[j], Rhyp = hypo_dist, comp = "LRG", VS30 = buildings_df$vs30)
  lnPGV_median_all[, j] <- pgv_result$lnPGV
}
cat("  Done.\n\n")

# Initialize storage for all L × M samples
# We'll collect results as a list and combine at the end
all_psi_virgin <- vector("list", L_GMM)
all_psi_predamage <- vector("list", L_GMM)
max_pgv_all <- matrix(0, n_buildings, L_GMM)

set.seed(42)
global_start <- Sys.time()

# =============================================================================
# OUTER LOOP: L GMM REALIZATIONS
# =============================================================================
for (l in 1:L_GMM) {
  cat(sprintf("GMM realization %d/%d\n", l, L_GMM))
  l_start <- Sys.time()

  # Generate GMM uncertainty terms for this realization
  # 1. Site terms (φS2S) - persist across all events for this GMM realization
  eta_site_l <- rnorm(n_buildings, 0, gmm_phiS2S)

  # 2. Event terms (τ) - one per earthquake
  eta_event_l <- rnorm(n_events, 0, gmm_tau)

  # Initialize N × M matrices for this GMM realization
  psi_virgin <- matrix(0, nrow = n_buildings, ncol = M_DRAWS)
  psi_predamage <- matrix(buildings_df$initial_psi, nrow = n_buildings, ncol = M_DRAWS)
  max_pgv_l <- rep(0, n_buildings)
  n_events_felt_l <- rep(0, n_buildings)

  # INNER LOOP: Events (with M posterior draws vectorized)
  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]

    # 3. Within-event terms (φSS) - different for each building-event
    eta_within <- rnorm(n_buildings, 0, gmm_phiSS)

    # Realized PGV for this GMM realization
    lnPGV_realized <- lnPGV_median_all[, j] + eta_event_l[j] + eta_site_l + eta_within
    pgv_event <- exp(lnPGV_realized) * 10  # mm/s

    # Update max PGV
    max_pgv_l <- pmax(max_pgv_l, pgv_event)

    # Find affected buildings
    affected <- which(pgv_event >= PGV_THRESHOLD)
    n_affected <- length(affected)

    if (n_affected > 0) {
      n_events_felt_l[affected] <- n_events_felt_l[affected] + 1

      # Predict damage for affected buildings (N_affected × M matrix)
      delta_virgin <- predict_damage_matrix(
        pgv_vec = pgv_event[affected],
        psi_matrix = psi_virgin[affected, , drop = FALSE]
      )
      psi_virgin[affected, ] <- pmin(psi_virgin[affected, ] + delta_virgin, PSI_MAX)

      delta_predamage <- predict_damage_matrix(
        pgv_vec = pgv_event[affected],
        psi_matrix = psi_predamage[affected, , drop = FALSE]
      )
      psi_predamage[affected, ] <- pmin(psi_predamage[affected, ] + delta_predamage, PSI_MAX)
    }
  }

  # Store results for this GMM realization
  all_psi_virgin[[l]] <- psi_virgin
  all_psi_predamage[[l]] <- psi_predamage
  max_pgv_all[, l] <- max_pgv_l

  l_elapsed <- as.numeric(difftime(Sys.time(), l_start, units = "mins"))
  cat(sprintf("  Completed in %.1f min\n", l_elapsed))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\nAll simulations complete! Total time: %.1f min\n\n", total_elapsed))

# =============================================================================
# 7. COMBINE L × M SAMPLES AND COMPUTE PERCENTILES
# =============================================================================

cat("Step 6: Combining", L_GMM, "×", M_DRAWS, "=", L_GMM * M_DRAWS, "samples and computing percentiles...\n")

# Combine all L GMM realizations into N × (L*M) matrices
# Each GMM realization l gives us an N × M matrix
psi_virgin_all <- do.call(cbind, all_psi_virgin)   # N × (L*M)
psi_predamage_all <- do.call(cbind, all_psi_predamage)  # N × (L*M)

cat("  Combined matrix size:", nrow(psi_virgin_all), "×", ncol(psi_virgin_all), "\n")

# Compute percentiles over ALL samples (GMM + posterior + stochastic)
buildings_df <- buildings_df |>
  mutate(
    # Point estimates (mean across ALL samples)
    psi_virgin_mean = rowMeans(psi_virgin_all),
    psi_predamage_mean = rowMeans(psi_predamage_all),

    # Full predictive percentiles (integrating GMM + posterior uncertainty)
    psi_virgin_p10 = rowQuantiles(psi_virgin_all, probs = 0.10),
    psi_virgin_p50 = rowQuantiles(psi_virgin_all, probs = 0.50),
    psi_virgin_p90 = rowQuantiles(psi_virgin_all, probs = 0.90),
    psi_virgin_sd = rowSds(psi_virgin_all),

    psi_predamage_p10 = rowQuantiles(psi_predamage_all, probs = 0.10),
    psi_predamage_p50 = rowQuantiles(psi_predamage_all, probs = 0.50),
    psi_predamage_p90 = rowQuantiles(psi_predamage_all, probs = 0.90),
    psi_predamage_sd = rowSds(psi_predamage_all),

    # Delta Ψ (for virgin, delta = psi; for predamage, delta = psi - initial)
    delta_psi_virgin_mean = psi_virgin_mean,
    delta_psi_predamage_mean = psi_predamage_mean - initial_psi,

    # Exceedance probabilities (fraction of ALL L×M samples)
    p_visible_virgin = rowMeans(psi_virgin_all >= 1.0),
    p_visible_predamage = rowMeans(psi_predamage_all >= 1.0),
    p_moderate_virgin = rowMeans(psi_virgin_all >= 2.0),
    p_moderate_predamage = rowMeans(psi_predamage_all >= 2.0),
    p_severe_virgin = rowMeans(psi_virgin_all >= 3.0),
    p_severe_predamage = rowMeans(psi_predamage_all >= 3.0),

    # PGV (max across all GMM realizations)
    max_pgv = rowMaxs(max_pgv_all),
    max_pgv_mean = rowMeans(max_pgv_all),
    max_pgv_p10 = rowQuantiles(max_pgv_all, probs = 0.10),
    max_pgv_p90 = rowQuantiles(max_pgv_all, probs = 0.90)
  )

# =============================================================================
# 8. RESULTS SUMMARY
# =============================================================================

cat("\n=== RESULTS SUMMARY ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("GMM realizations:", L_GMM, "\n")
cat("Posterior draws:", M_DRAWS, "\n")
cat("Total samples per building:", L_GMM * M_DRAWS, "\n")
cat("Earthquakes:", n_events, "\n\n")

cat("--- VIRGIN WALLS ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_virgin_mean), 3), "\n")
cat("  p10:", round(mean(buildings_df$psi_virgin_p10), 3),
    "| p50:", round(mean(buildings_df$psi_virgin_p50), 3),
    "| p90:", round(mean(buildings_df$psi_virgin_p90), 3), "\n")
cat("% visible (Ψ≥1):", round(100 * mean(buildings_df$p_visible_virgin), 2), "%\n")
cat("% moderate (Ψ≥2):", round(100 * mean(buildings_df$p_moderate_virgin), 2), "%\n\n")

cat("--- PRE-DAMAGE ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_predamage_mean), 3), "\n")
cat("  p10:", round(mean(buildings_df$psi_predamage_p10), 3),
    "| p50:", round(mean(buildings_df$psi_predamage_p50), 3),
    "| p90:", round(mean(buildings_df$psi_predamage_p90), 3), "\n")
cat("% visible (Ψ≥1):", round(100 * mean(buildings_df$p_visible_predamage), 2), "%\n\n")

cat("--- UNCERTAINTY MAGNITUDE ---\n")
cat("Mean predictive SD (virgin):", round(mean(buildings_df$psi_virgin_sd), 3), "\n")
cat("Mean p90/p50 ratio (virgin):",
    round(mean(buildings_df$psi_virgin_p90 / pmax(buildings_df$psi_virgin_p50, 0.01)), 2), "\n")
cat("Mean p90-p10 range (virgin):",
    round(mean(buildings_df$psi_virgin_p90 - buildings_df$psi_virgin_p10), 3), "\n\n")

cat("NOTE: Credible intervals integrate:\n")
cat("  - GMM uncertainty (τ, φS2S, φSS)\n")
cat("  - Fragility posterior uncertainty (M draws)\n")
cat("  - Stochastic hurdle-gamma sampling\n\n")

# By zone
by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    psi_virgin_mean = mean(psi_virgin_mean),
    psi_virgin_p10 = mean(psi_virgin_p10),
    psi_virgin_p90 = mean(psi_virgin_p90),
    psi_virgin_sd = mean(psi_virgin_sd),
    p_visible_virgin = mean(p_visible_virgin),
    .groups = "drop"
  )

cat("By distance zone:\n")
print(by_zone)

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("\nStep 7: Saving results...\n")

results <- buildings_df |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, zone, age_category,
    initial_psi,
    # Point estimates
    psi_virgin_mean, psi_predamage_mean,
    delta_psi_virgin_mean, delta_psi_predamage_mean,
    # Full predictive percentiles (GMM + posterior + stochastic)
    psi_virgin_p10, psi_virgin_p50, psi_virgin_p90, psi_virgin_sd,
    psi_predamage_p10, psi_predamage_p50, psi_predamage_p90, psi_predamage_sd,
    # Exceedance probabilities
    p_visible_virgin, p_visible_predamage,
    p_moderate_virgin, p_moderate_predamage,
    p_severe_virgin, p_severe_predamage,
    # PGV with uncertainty
    max_pgv, max_pgv_mean, max_pgv_p10, max_pgv_p90
  )

# Choose filename based on sample fraction
if (SAMPLE_FRAC < 1.0) {
  output_file <- here("outputs", "models", "spatial_damage_uncertainty_sample.rds")
} else {
  output_file <- here("outputs", "models", "spatial_damage_uncertainty_full.rds")
}

saveRDS(results, output_file)
cat("Saved:", output_file, "\n")

# Summary tables
summary_tables <- list(
  by_zone = by_zone,
  by_age = buildings_df |>
    group_by(age_category) |>
    summarise(
      n = n(),
      psi_virgin_mean = mean(psi_virgin_mean),
      psi_virgin_p10 = mean(psi_virgin_p10),
      psi_virgin_p90 = mean(psi_virgin_p90),
      psi_virgin_sd = mean(psi_virgin_sd),
      p_visible_virgin = mean(p_visible_virgin),
      .groups = "drop"
    ),
  overall = tibble(
    n_buildings = n_buildings,
    n_events = n_events,
    l_gmm = L_GMM,
    m_draws = M_DRAWS,
    total_samples = L_GMM * M_DRAWS,
    psi_virgin_mean = mean(buildings_df$psi_virgin_mean),
    psi_virgin_p10 = mean(buildings_df$psi_virgin_p10),
    psi_virgin_p90 = mean(buildings_df$psi_virgin_p90),
    psi_virgin_sd = mean(buildings_df$psi_virgin_sd),
    p_visible_virgin = mean(buildings_df$p_visible_virgin)
  )
)

summary_file <- gsub("\\.rds$", "_summaries.rds", output_file)
saveRDS(summary_tables, summary_file)
cat("Saved:", summary_file, "\n")

cat("\n=== DONE ===\n")
