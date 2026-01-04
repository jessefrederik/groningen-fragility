# =============================================================================
# 11d_spatial_damage_ratchet.R
#
# RATCHET MODEL for Groningen building damage predictions
#
# KEY INSIGHT: FEM data shows N=8 repetitions give essentially the same damage
# as N=1 (ratio ≈ 1.00-1.09). This implies damage only increases when loading
# EXCEEDS previous maximum - the "ratchet" mechanism.
#
# CRITICAL CHANGE FROM 11c:
# - OLD (additive): ΔΨ = predict(PGV, current_Ψ); Ψ += ΔΨ for each event
# - NEW (ratchet):  Only add damage when PGV > max_PGV_so_far
#                   ΔΨ = damage(new_max) - damage(old_max) at baseline Ψ
#
# IMPLEMENTATION NOTES (per ChatGPT review):
# 1. Use conditional mean E[ΔΨ] = (1 - p_zero) × μ, not stochastic sampling
#    → Avoids negative deltas from differencing stochastic draws
# 2. The "exceeds max" mask is same for all M posterior draws within GMM-l
#    → PGV depends only on GMM, not posterior
# 3. Compute damage at baseline initial_psi (0 for virgin, building's initial
#    for predamage), not at current accumulated Ψ
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
PGV_THRESHOLD <- 2.0    # mm/s - lower threshold for ratchet model
L_GMM <- 10             # Number of GMM realizations (outer loop)
M_DRAWS <- 50           # Number of posterior draws (inner loop)
SAMPLE_FRAC <- 1.0      # Full dataset (set to 0.1 for validation)

cat("=== Spatial Damage Prediction with RATCHET MODEL ===\n")
cat("Configuration:\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  Posterior draws (M):", M_DRAWS, "\n")
cat("  Total samples per building:", L_GMM * M_DRAWS, "\n")
cat("  Sample fraction:", SAMPLE_FRAC, "\n")
cat("  PGV threshold:", PGV_THRESHOLD, "mm/s\n")
cat("  Model: RATCHET (damage only when exceeding previous max)\n\n")

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

buildings$dist_to_center <- calc_dist_km(buildings$lat, buildings$lon, 53.35, 6.72)

# Stratified sampling for validation
if (SAMPLE_FRAC < 1.0) {
  cat("Stratified sampling:", SAMPLE_FRAC * 100, "% of buildings\n")
  set.seed(123)
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
gmm_tau <- gmm_coef$tau
gmm_phiS2S <- gmm_coef$phiS2S
gmm_phiSS <- gmm_coef$phiSS

cat("  τ (inter-event):", round(gmm_tau, 4), "\n")
cat("  φS2S (site-to-site):", round(gmm_phiS2S, 4), "\n")
cat("  φSS (within-event):", round(gmm_phiSS, 4), "\n\n")

# =============================================================================
# 4. LOAD BRMS MODEL AND EXTRACT POSTERIORS
# =============================================================================

cat("Step 4: Loading brms model and sampling posteriors...\n")

fit_brms <- readRDS(here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here("outputs", "models", "standardization_params.rds"))

post <- as_draws_df(fit_brms)
n_post <- nrow(post)

set.seed(42)
posterior_idx <- sample(n_post, M_DRAWS, replace = FALSE)
cat("  Sampled", M_DRAWS, "posterior draws from", n_post, "total\n")

# Extract parameters for selected draws
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
pgv_cumsum <- pgv_cumsum_full[posterior_idx, ]

hu_pgv_simplex_full <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) hu_pgv_simplex_full[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
hu_pgv_cumsum_full <- pgv_K * t(apply(hu_pgv_simplex_full, 1, cumsum))
hu_pgv_cumsum <- hu_pgv_cumsum_full[posterior_idx, ]

mat_K <- 2
mat_simplex_full <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) mat_simplex_full[, i] <- post[[paste0("simo_momaterial_ord1[", i, "]")]]
mat_cumsum_full <- mat_K * t(apply(mat_simplex_full, 1, cumsum))
mat_eff <- mat_cumsum_full[posterior_idx, 1]

hu_mat_simplex_full <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) hu_mat_simplex_full[, i] <- post[[paste0("simo_hu_momaterial_ord1[", i, "]")]]
hu_mat_cumsum_full <- mat_K * t(apply(hu_mat_simplex_full, 1, cumsum))
hu_mat_eff <- hu_mat_cumsum_full[posterior_idx, 1]

cat("  Posterior parameters extracted\n\n")

# =============================================================================
# 5. PREDICTION FUNCTIONS
# =============================================================================

pgv_to_idx <- function(pgv) {
  idx <- findInterval(pgv, c(-Inf, 3, 6, 12, 24, 48, 80, 112, Inf))
  pmax(1, pmin(idx, 8))
}

# =============================================================================
# CONDITIONAL MEAN PREDICTION (for ratchet differencing)
#
# Returns E[ΔΨ | PGV, initial_psi] = (1 - p_zero) × μ
# This is deterministic given posterior draw, avoiding stochastic noise
# that would cause problems when differencing (old vs new).
#
# pgv_vec: length N vector of PGV values (mm/s)
# psi_baseline: length N vector of BASELINE initial_psi (not current accumulated)
# Returns: N × M matrix of conditional mean damage
# =============================================================================
predict_damage_mean_matrix <- function(pgv_vec, psi_baseline) {
  N <- length(pgv_vec)
  M <- M_DRAWS

  pgv_idx <- pgv_to_idx(pgv_vec)

  # Center psi at baseline (N-vector, replicated to N × M)
  psi_c <- psi_baseline - std_params$initial_psi_mean
  psi_c_mat <- matrix(psi_c, N, M)

  # PGV effect matrices
  pgv_eff_mat <- matrix(0, N, M)
  hu_pgv_eff_mat <- matrix(0, N, M)
  for (k in 2:8) {
    mask <- pgv_idx == k
    if (any(mask)) {
      pgv_eff_mat[mask, ] <- matrix(pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
      hu_pgv_eff_mat[mask, ] <- matrix(hu_pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
    }
  }

  # Material effect (same for all buildings)
  mat_eff_mat <- matrix(mat_eff, N, M, byrow = TRUE)
  hu_mat_eff_mat <- matrix(hu_mat_eff, N, M, byrow = TRUE)

  # Random effects: sample from marginal distributions
  # (For mean prediction, we could use 0, but sampling gives predictive distribution)
  r_eq <- matrix(rnorm(N * M, 0, rep(sd_eq, each = N)), N, M)
  r_fs <- matrix(rnorm(N * M, 0, rep(sd_fs, each = N)), N, M)
  hu_r_eq <- matrix(rnorm(N * M, 0, rep(hu_sd_eq, each = N)), N, M)
  hu_r_fs <- matrix(rnorm(N * M, 0, rep(hu_sd_fs, each = N)), N, M)

  # Fixed effects matrices
  b_int_mat <- matrix(b_intercept, N, M, byrow = TRUE)
  bsp_pgv_mat <- matrix(bsp_pgv, N, M, byrow = TRUE)
  bsp_mat_mat <- matrix(bsp_material, N, M, byrow = TRUE)
  b_psi_mat <- matrix(b_initial_psi, N, M, byrow = TRUE)

  hu_int_mat <- matrix(hu_intercept, N, M, byrow = TRUE)
  hu_bsp_pgv_mat <- matrix(hu_bsp_pgv, N, M, byrow = TRUE)
  hu_bsp_mat_mat <- matrix(hu_bsp_material, N, M, byrow = TRUE)
  hu_b_psi_mat <- matrix(hu_initial_psi, N, M, byrow = TRUE)

  # Linear predictors
  eta_mu <- b_int_mat + bsp_pgv_mat * pgv_eff_mat + bsp_mat_mat * mat_eff_mat +
            b_psi_mat * psi_c_mat + r_eq + r_fs

  eta_hu <- hu_int_mat + hu_bsp_pgv_mat * hu_pgv_eff_mat + hu_bsp_mat_mat * hu_mat_eff_mat +
            hu_b_psi_mat * psi_c_mat + hu_r_eq + hu_r_fs

  # Hurdle probability (probability of ZERO damage)
  p_zero <- 1 / (1 + exp(-eta_hu))  # plogis(eta_hu)

  # Gamma mean
  mu <- exp(eta_mu)

  # Conditional mean: E[ΔΨ] = (1 - p_zero) × μ
  delta_mean <- (1 - p_zero) * mu

  # Apply PGV threshold
  below_threshold <- pgv_vec < PGV_THRESHOLD
  delta_mean[below_threshold, ] <- 0

  return(delta_mean)
}

# =============================================================================
# 6. MAIN SIMULATION LOOP - RATCHET MODEL
# =============================================================================

cat("Step 5: Running RATCHET simulation...\n")
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("  Events:", n_events, "\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  Posterior draws (M):", M_DRAWS, "\n\n")

# Sort events chronologically
event_order <- order(groningen_eq$time_utc)
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

# Pre-compute median lnPGV for all events
cat("  Pre-computing median lnPGV for all events...\n")
lnPGV_median_all <- matrix(NA, n_buildings, n_events)
for (j in 1:n_events) {
  epi_dist <- calc_dist_km(buildings_df$lat, buildings_df$lon, eq_lat[j], eq_lon[j])
  hypo_dist <- sqrt(epi_dist^2 + eq_depth[j]^2)
  pgv_result <- pgv_model(M = eq_mag[j], Rhyp = hypo_dist, comp = "LRG", VS30 = buildings_df$vs30)
  lnPGV_median_all[, j] <- pgv_result$lnPGV
}
cat("  Done.\n\n")

# Storage
all_psi_virgin <- vector("list", L_GMM)
all_psi_predamage <- vector("list", L_GMM)
max_pgv_all <- matrix(0, n_buildings, L_GMM)
n_contributing_all <- matrix(0, n_buildings, L_GMM)

set.seed(42)
global_start <- Sys.time()

# =============================================================================
# OUTER LOOP: L GMM REALIZATIONS
# =============================================================================
for (l in 1:L_GMM) {
  cat(sprintf("GMM realization %d/%d\n", l, L_GMM))
  l_start <- Sys.time()

  # GMM uncertainty terms for this realization
  eta_site_l <- rnorm(n_buildings, 0, gmm_phiS2S)
  eta_event_l <- rnorm(n_events, 0, gmm_tau)

  # Initialize: N × M matrices
  # For ratchet, we track damage-at-max-pgv, not accumulated damage
  psi_virgin <- matrix(0, nrow = n_buildings, ncol = M_DRAWS)
  psi_predamage <- matrix(buildings_df$initial_psi, nrow = n_buildings, ncol = M_DRAWS)

  # Track max PGV per building (same for all M draws within GMM-l)
  max_pgv_l <- rep(0, n_buildings)
  n_contributing_l <- rep(0, n_buildings)

  # Baseline psi (for computing damage curves)
  psi0_virgin <- rep(0, n_buildings)
  psi0_predamage <- buildings_df$initial_psi

  # INNER LOOP: Events (chronological)
  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]

    # Within-event term
    eta_within <- rnorm(n_buildings, 0, gmm_phiSS)

    # Realized PGV for this GMM realization (same for all M draws!)
    lnPGV_realized <- lnPGV_median_all[, j] + eta_event_l[j] + eta_site_l + eta_within
    pgv_event <- exp(lnPGV_realized) * 10  # mm/s

    # =======================================================================
    # RATCHET LOGIC: Only process buildings where this event exceeds their max
    # =======================================================================
    exceeds <- pgv_event > max_pgv_l
    n_exceeds <- sum(exceeds)

    if (n_exceeds > 0) {
      # These buildings set a new max PGV
      idx_exceeds <- which(exceeds)

      # Old and new PGV for exceeding buildings
      pgv_old <- max_pgv_l[idx_exceeds]
      pgv_new <- pgv_event[idx_exceeds]

      # Compute damage at OLD max (baseline psi)
      # For buildings with pgv_old = 0 (first event), this gives 0
      inc_old_virgin <- predict_damage_mean_matrix(pgv_old, psi0_virgin[idx_exceeds])
      inc_old_predamage <- predict_damage_mean_matrix(pgv_old, psi0_predamage[idx_exceeds])

      # Compute damage at NEW max (baseline psi)
      inc_new_virgin <- predict_damage_mean_matrix(pgv_new, psi0_virgin[idx_exceeds])
      inc_new_predamage <- predict_damage_mean_matrix(pgv_new, psi0_predamage[idx_exceeds])

      # Incremental damage = new - old (clamped to >= 0)
      delta_virgin <- pmax(inc_new_virgin - inc_old_virgin, 0)
      delta_predamage <- pmax(inc_new_predamage - inc_old_predamage, 0)

      # Update psi
      psi_virgin[idx_exceeds, ] <- pmin(psi_virgin[idx_exceeds, ] + delta_virgin, PSI_MAX)
      psi_predamage[idx_exceeds, ] <- pmin(psi_predamage[idx_exceeds, ] + delta_predamage, PSI_MAX)

      # Update max PGV and contribution count
      max_pgv_l[idx_exceeds] <- pgv_new
      n_contributing_l[idx_exceeds] <- n_contributing_l[idx_exceeds] + 1
    }
  }

  # Store results
  all_psi_virgin[[l]] <- psi_virgin
  all_psi_predamage[[l]] <- psi_predamage
  max_pgv_all[, l] <- max_pgv_l
  n_contributing_all[, l] <- n_contributing_l

  l_elapsed <- as.numeric(difftime(Sys.time(), l_start, units = "mins"))
  cat(sprintf("  Completed in %.1f min (avg %.1f contributing events/building)\n",
              l_elapsed, mean(n_contributing_l)))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\nAll simulations complete! Total time: %.1f min\n\n", total_elapsed))

# =============================================================================
# 7. COMBINE SAMPLES AND COMPUTE STATISTICS
# =============================================================================

cat("Step 6: Computing statistics...\n")

psi_virgin_all <- do.call(cbind, all_psi_virgin)
psi_predamage_all <- do.call(cbind, all_psi_predamage)

buildings_df <- buildings_df |>
  mutate(
    psi_virgin_mean = rowMeans(psi_virgin_all),
    psi_predamage_mean = rowMeans(psi_predamage_all),

    psi_virgin_p10 = rowQuantiles(psi_virgin_all, probs = 0.10),
    psi_virgin_p50 = rowQuantiles(psi_virgin_all, probs = 0.50),
    psi_virgin_p90 = rowQuantiles(psi_virgin_all, probs = 0.90),
    psi_virgin_sd = rowSds(psi_virgin_all),

    psi_predamage_p10 = rowQuantiles(psi_predamage_all, probs = 0.10),
    psi_predamage_p50 = rowQuantiles(psi_predamage_all, probs = 0.50),
    psi_predamage_p90 = rowQuantiles(psi_predamage_all, probs = 0.90),
    psi_predamage_sd = rowSds(psi_predamage_all),

    delta_psi_virgin_mean = psi_virgin_mean,
    delta_psi_predamage_mean = psi_predamage_mean - initial_psi,

    p_visible_virgin = rowMeans(psi_virgin_all >= 1.0),
    p_visible_predamage = rowMeans(psi_predamage_all >= 1.0),
    p_moderate_virgin = rowMeans(psi_virgin_all >= 2.0),
    p_moderate_predamage = rowMeans(psi_predamage_all >= 2.0),
    p_severe_virgin = rowMeans(psi_virgin_all >= 3.0),
    p_severe_predamage = rowMeans(psi_predamage_all >= 3.0),

    max_pgv = rowMaxs(max_pgv_all),
    max_pgv_mean = rowMeans(max_pgv_all),
    max_pgv_p10 = rowQuantiles(max_pgv_all, probs = 0.10),
    max_pgv_p90 = rowQuantiles(max_pgv_all, probs = 0.90),

    n_contributing_mean = rowMeans(n_contributing_all),
    n_contributing_max = rowMaxs(n_contributing_all)
  )

# =============================================================================
# 8. RESULTS SUMMARY
# =============================================================================

cat("\n=== RATCHET MODEL RESULTS ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("Earthquakes:", n_events, "\n")
cat("Total samples per building:", L_GMM * M_DRAWS, "\n\n")

cat("--- VIRGIN WALLS ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_virgin_mean), 4), "\n")
cat("  p10:", round(mean(buildings_df$psi_virgin_p10), 4),
    "| p50:", round(mean(buildings_df$psi_virgin_p50), 4),
    "| p90:", round(mean(buildings_df$psi_virgin_p90), 4), "\n")
cat("% visible (Ψ≥1):", round(100 * mean(buildings_df$p_visible_virgin), 2), "%\n")
cat("% moderate (Ψ≥2):", round(100 * mean(buildings_df$p_moderate_virgin), 2), "%\n\n")

cat("--- RATCHET STATISTICS ---\n")
cat("Mean contributing events:", round(mean(buildings_df$n_contributing_mean), 2), "\n")
cat("Max contributing events:", max(buildings_df$n_contributing_max), "\n")
cat("Mean max PGV:", round(mean(buildings_df$max_pgv_mean), 2), "mm/s\n\n")

# By zone
by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    psi_virgin_mean = mean(psi_virgin_mean),
    psi_virgin_p90 = mean(psi_virgin_p90),
    p_visible = mean(p_visible_virgin),
    n_contrib = mean(n_contributing_mean),
    max_pgv = mean(max_pgv_mean),
    .groups = "drop"
  )

cat("By distance zone:\n")
print(by_zone)

# =============================================================================
# 9. VALIDATION CHECKS (per ChatGPT recommendations)
# =============================================================================

cat("\n=== VALIDATION CHECKS ===\n")

# Check 1: Contribution count should be small
cat("\n1. CONTRIBUTION COUNT (should be << total events):\n")
cat("   Total events:", n_events, "\n")
cat("   Mean contributing:", round(mean(buildings_df$n_contributing_mean), 2), "\n")
cat("   Max contributing:", max(buildings_df$n_contributing_max), "\n")
if (mean(buildings_df$n_contributing_mean) < n_events / 5) {
  cat("   ✓ PASS: Most events don't contribute (ratchet working)\n")
} else {
  cat("   ⚠ WARN: Many events contributing (check ratchet logic)\n")
}

# Check 2: Final Ψ should be close to single-event at max PGV
cat("\n2. CATALOGUE SANITY (Ψ ≈ damage at max PGV):\n")
# For buildings with high max_pgv, check if damage is reasonable
high_pgv_buildings <- buildings_df |> filter(max_pgv_mean > 50)
if (nrow(high_pgv_buildings) > 0) {
  cat("   Buildings with max PGV > 50 mm/s:", nrow(high_pgv_buildings), "\n")
  cat("   Mean Ψ for these:", round(mean(high_pgv_buildings$psi_virgin_mean), 2), "\n")
  cat("   (Should be ~1-2 based on FEM at PGV=64)\n")
}

# Check 3: Compare to additive (if available)
additive_file <- here("outputs", "models", "spatial_damage_uncertainty_full.rds")
if (file.exists(additive_file)) {
  additive <- readRDS(additive_file)
  cat("\n3. COMPARISON TO ADDITIVE MODEL:\n")
  cat("   Additive mean Ψ:", round(mean(additive$psi_virgin_mean), 4), "\n")
  cat("   Ratchet mean Ψ:", round(mean(buildings_df$psi_virgin_mean), 4), "\n")
  cat("   Reduction factor:", round(mean(additive$psi_virgin_mean) / mean(buildings_df$psi_virgin_mean), 1), "x\n")
}

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("\nStep 7: Saving results...\n")

results <- buildings_df |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, zone, age_category, initial_psi,
    psi_virgin_mean, psi_predamage_mean,
    delta_psi_virgin_mean, delta_psi_predamage_mean,
    psi_virgin_p10, psi_virgin_p50, psi_virgin_p90, psi_virgin_sd,
    psi_predamage_p10, psi_predamage_p50, psi_predamage_p90, psi_predamage_sd,
    p_visible_virgin, p_visible_predamage,
    p_moderate_virgin, p_moderate_predamage,
    p_severe_virgin, p_severe_predamage,
    max_pgv, max_pgv_mean, max_pgv_p10, max_pgv_p90,
    n_contributing_mean, n_contributing_max
  )

output_file <- here("outputs", "models", "spatial_damage_ratchet.rds")
saveRDS(results, output_file)
cat("Saved:", output_file, "\n")

summary_tables <- list(
  by_zone = by_zone,
  by_age = buildings_df |>
    group_by(age_category) |>
    summarise(
      n = n(),
      psi_virgin_mean = mean(psi_virgin_mean),
      p_visible = mean(p_visible_virgin),
      n_contrib = mean(n_contributing_mean),
      .groups = "drop"
    ),
  overall = tibble(
    model = "ratchet",
    n_buildings = n_buildings,
    n_events = n_events,
    psi_virgin_mean = mean(buildings_df$psi_virgin_mean),
    psi_virgin_p90 = mean(buildings_df$psi_virgin_p90),
    p_visible = mean(buildings_df$p_visible_virgin),
    n_contributing_mean = mean(buildings_df$n_contributing_mean)
  )
)

saveRDS(summary_tables, here("outputs", "models", "spatial_damage_ratchet_summaries.rds"))

cat("\n=== DONE ===\n")
