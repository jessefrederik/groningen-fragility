# =============================================================================
# 03d_spatial_daf.R
#
# KORSWAGEN DAF (Damage Accumulation Function) implementation
# Uses brms hurdle-gamma model as surrogate with history-aware differencing
#
# KEY FEATURES:
# 1. Two-regime update rule (Korswagen Eq. 6.7):
#    - Record-breaking: ΔΨ = Surrogate(Ψ, V_new, N=2) - Surrogate(Ψ, V_max, N=1)
#    - Similar events:  ΔΨ = Surrogate(Ψ, V, N_eq+1) - Surrogate(Ψ, V, N_eq)
#
# 2. N_eq computation:
#    - Count prior events "similar" to current (within ρ=0.25 tolerance)
#
# 3. Uses brms N-coefficients (previously ignored in 11c/11d)
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(sf)
library(brms)
library(here)
library(matrixStats)

source(here("helperfuncties", "fetch_knmi_events.R"))
source(here("helperfuncties", "pgv_model_original.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

PSI_MAX <- 4.0
MAG_THRESHOLD <- 2.0
PGV_THRESHOLD <- 2.0    # mm/s
L_GMM <- 10
M_DRAWS <- 50
SAMPLE_FRAC <- 1.0
RHO <- 0.25             # Korswagen similarity tolerance

cat("=== Spatial Damage Prediction with KORSWAGEN DAF ===\n")
cat("Configuration:\n")
cat("  GMM realizations (L):", L_GMM, "\n")
cat("  Posterior draws (M):", M_DRAWS, "\n")
cat("  PGV threshold:", PGV_THRESHOLD, "mm/s\n")
cat("  Similarity tolerance (ρ):", RHO, "\n")
cat("  Model: DAF (two-regime with N_eq)\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("Step 1: Loading data...\n")

eq_cache <- here("outputs", "models", "groningen_earthquakes.rds")
groningen_eq <- readRDS(eq_cache)
n_events <- nrow(groningen_eq)
cat("  Earthquakes:", n_events, "\n")

buildings <- st_read(here("datafiles", "bag_selectie.gpkg"), quiet = TRUE)
vs30 <- read_delim(
  here("datafiles", "vs30.csv"), delim = ";",
  locale = locale(decimal_mark = ","),
  col_names = c("postcode", "vs30", "x1", "x2", "x3"),
  skip = 1, show_col_types = FALSE
) |> select(postcode, vs30) |> mutate(postcode = as.character(postcode))

buildings <- buildings |>
  mutate(pc4 = substr(postcode, 1, 4)) |>
  left_join(vs30, by = c("pc4" = "postcode"))
buildings$vs30[is.na(buildings$vs30)] <- mean(buildings$vs30, na.rm = TRUE)

buildings_wgs84 <- st_transform(buildings, 4326)
coords <- st_coordinates(st_centroid(st_geometry(buildings_wgs84)))
buildings$lon <- coords[, 1]
buildings$lat <- coords[, 2]

calc_dist_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

buildings$dist_to_center <- calc_dist_km(buildings$lat, buildings$lon, 53.35, 6.72)

if (SAMPLE_FRAC < 1.0) {
  set.seed(123)
  n_sample <- round(nrow(buildings) * SAMPLE_FRAC)
  buildings <- slice_sample(buildings, n = n_sample)
}

n_buildings <- nrow(buildings)
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")

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
cat("  Data loaded\n\n")

# =============================================================================
# 2. GMM SETUP
# =============================================================================

cat("Step 2: Setting up GMM...\n")
gmm_coef <- get_pgv_coefficients("LRG")
gmm_tau <- gmm_coef$tau
gmm_phiS2S <- gmm_coef$phiS2S
gmm_phiSS <- gmm_coef$phiSS
cat("  τ:", round(gmm_tau, 4), "| φS2S:", round(gmm_phiS2S, 4), "| φSS:", round(gmm_phiSS, 4), "\n\n")

# =============================================================================
# 3. LOAD BRMS MODEL AND EXTRACT POSTERIORS (INCLUDING N-EFFECT)
# =============================================================================

cat("Step 3: Loading brms model with N-effect...\n")

fit_brms <- readRDS(here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here("outputs", "models", "standardization_params.rds"))

post <- as_draws_df(fit_brms)
n_post <- nrow(post)

set.seed(42)
posterior_idx <- sample(n_post, M_DRAWS, replace = FALSE)

# Standard parameters
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

# PGV monotonic (K=7)
pgv_K <- 7
pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) pgv_simplex[, i] <- post[[paste0("simo_mopgv_ord1[", i, "]")]]
pgv_cumsum <- pgv_K * t(apply(pgv_simplex, 1, cumsum))[posterior_idx, ]

hu_pgv_simplex <- matrix(NA, n_post, pgv_K)
for (i in 1:pgv_K) hu_pgv_simplex[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
hu_pgv_cumsum <- pgv_K * t(apply(hu_pgv_simplex, 1, cumsum))[posterior_idx, ]

# Material monotonic (K=2)
mat_K <- 2
mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) mat_simplex[, i] <- post[[paste0("simo_momaterial_ord1[", i, "]")]]
mat_cumsum <- mat_K * t(apply(mat_simplex, 1, cumsum))
mat_eff <- mat_cumsum[posterior_idx, 1]

hu_mat_simplex <- matrix(NA, n_post, mat_K)
for (i in 1:mat_K) hu_mat_simplex[, i] <- post[[paste0("simo_hu_momaterial_ord1[", i, "]")]]
hu_mat_cumsum <- mat_K * t(apply(hu_mat_simplex, 1, cumsum))
hu_mat_eff <- hu_mat_cumsum[posterior_idx, 1]

# ===== NEW: N-EFFECT COEFFICIENTS =====
bsp_n <- post$bsp_mon_ord[posterior_idx]           # N coefficient (severity)
hu_bsp_n <- post$bsp_hu_mon_ord[posterior_idx]     # N coefficient (hurdle)

# N monotonic simplex (K=4 for 5 levels: 1,2,3,4,8)
n_K <- 4
n_simplex <- matrix(NA, n_post, n_K)
for (i in 1:n_K) n_simplex[, i] <- post[[paste0("simo_mon_ord1[", i, "]")]]
n_cumsum <- n_K * t(apply(n_simplex, 1, cumsum))[posterior_idx, ]

hu_n_simplex <- matrix(NA, n_post, n_K)
for (i in 1:n_K) hu_n_simplex[, i] <- post[[paste0("simo_hu_mon_ord1[", i, "]")]]
hu_n_cumsum <- n_K * t(apply(hu_n_simplex, 1, cumsum))[posterior_idx, ]

cat("  Extracted", M_DRAWS, "posterior draws including N-effect\n")
cat("  Mean N coefficient (severity):", round(mean(bsp_n), 4), "\n")
cat("  N=8 vs N=1 damage ratio:", round(exp(mean(bsp_n) * 4), 3), "\n\n")

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

pgv_to_idx <- function(pgv) {
  idx <- findInterval(pgv, c(-Inf, 3, 6, 12, 24, 48, 80, 112, Inf))
  pmax(1, pmin(idx, 8))
}

# Map N_eq to ordered factor index (1=N1, 2=N2, 3=N3, 4=N4, 5=N8)
n_eq_to_idx <- function(n_eq) {
  idx <- case_when(
    n_eq <= 1 ~ 1L,
    n_eq <= 2 ~ 2L,
    n_eq <= 3 ~ 3L,
    n_eq <= 5 ~ 4L,
    TRUE ~ 5L
  )
  return(idx)
}

# =============================================================================
# PREDICT DAMAGE WITH N-EFFECT
#
# Returns E[ΔΨ | PGV, Ψ_current, N] = (1 - p_zero) × μ
# =============================================================================
predict_damage_mean_with_n <- function(pgv_vec, psi_vec, n_eq_vec) {
  N <- length(pgv_vec)
  M <- M_DRAWS

  pgv_idx <- pgv_to_idx(pgv_vec)
  n_idx <- n_eq_to_idx(n_eq_vec)

  # Center psi
  psi_c <- psi_vec - std_params$initial_psi_mean
  psi_c_mat <- matrix(psi_c, N, M)

  # PGV effect
  pgv_eff_mat <- matrix(0, N, M)
  hu_pgv_eff_mat <- matrix(0, N, M)
  for (k in 2:8) {
    mask <- pgv_idx == k
    if (any(mask)) {
      pgv_eff_mat[mask, ] <- matrix(pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
      hu_pgv_eff_mat[mask, ] <- matrix(hu_pgv_cumsum[, k-1], sum(mask), M, byrow = TRUE)
    }
  }

  # N effect (NEW!)
  n_eff_mat <- matrix(0, N, M)
  hu_n_eff_mat <- matrix(0, N, M)
  for (k in 2:5) {
    mask <- n_idx == k
    if (any(mask)) {
      n_eff_mat[mask, ] <- matrix(n_cumsum[, k-1], sum(mask), M, byrow = TRUE)
      hu_n_eff_mat[mask, ] <- matrix(hu_n_cumsum[, k-1], sum(mask), M, byrow = TRUE)
    }
  }

  # Material effect
  mat_eff_mat <- matrix(mat_eff, N, M, byrow = TRUE)
  hu_mat_eff_mat <- matrix(hu_mat_eff, N, M, byrow = TRUE)

  # Random effects (marginal)
  r_eq <- matrix(rnorm(N * M, 0, rep(sd_eq, each = N)), N, M)
  r_fs <- matrix(rnorm(N * M, 0, rep(sd_fs, each = N)), N, M)
  hu_r_eq <- matrix(rnorm(N * M, 0, rep(hu_sd_eq, each = N)), N, M)
  hu_r_fs <- matrix(rnorm(N * M, 0, rep(hu_sd_fs, each = N)), N, M)

  # Fixed effects
  b_int_mat <- matrix(b_intercept, N, M, byrow = TRUE)
  bsp_pgv_mat <- matrix(bsp_pgv, N, M, byrow = TRUE)
  bsp_n_mat <- matrix(bsp_n, N, M, byrow = TRUE)
  bsp_mat_mat <- matrix(bsp_material, N, M, byrow = TRUE)
  b_psi_mat <- matrix(b_initial_psi, N, M, byrow = TRUE)

  hu_int_mat <- matrix(hu_intercept, N, M, byrow = TRUE)
  hu_bsp_pgv_mat <- matrix(hu_bsp_pgv, N, M, byrow = TRUE)
  hu_bsp_n_mat <- matrix(hu_bsp_n, N, M, byrow = TRUE)
  hu_bsp_mat_mat <- matrix(hu_bsp_material, N, M, byrow = TRUE)
  hu_b_psi_mat <- matrix(hu_initial_psi, N, M, byrow = TRUE)

  # Linear predictors (NOW INCLUDES N!)
  eta_mu <- b_int_mat +
            bsp_pgv_mat * pgv_eff_mat +
            bsp_n_mat * n_eff_mat +        # <-- NEW
            bsp_mat_mat * mat_eff_mat +
            b_psi_mat * psi_c_mat +
            r_eq + r_fs

  eta_hu <- hu_int_mat +
            hu_bsp_pgv_mat * hu_pgv_eff_mat +
            hu_bsp_n_mat * hu_n_eff_mat +  # <-- NEW
            hu_bsp_mat_mat * hu_mat_eff_mat +
            hu_b_psi_mat * psi_c_mat +
            hu_r_eq + hu_r_fs

  # Conditional mean
  p_zero <- 1 / (1 + exp(-eta_hu))
  mu <- exp(eta_mu)
  delta_mean <- (1 - p_zero) * mu

  # Threshold
  below <- pgv_vec < PGV_THRESHOLD
  delta_mean[below, ] <- 0

  return(delta_mean)
}

# =============================================================================
# 5. MAIN SIMULATION
# =============================================================================

cat("Step 4: Running DAF simulation...\n")
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("  Events:", n_events, "\n\n")

event_order <- order(groningen_eq$time_utc)
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

# Pre-compute median lnPGV
cat("  Pre-computing median lnPGV...\n")
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
n_record_all <- matrix(0, n_buildings, L_GMM)
n_similar_all <- matrix(0, n_buildings, L_GMM)

set.seed(42)
global_start <- Sys.time()

for (l in 1:L_GMM) {
  cat(sprintf("GMM realization %d/%d\n", l, L_GMM))
  l_start <- Sys.time()

  eta_site_l <- rnorm(n_buildings, 0, gmm_phiS2S)
  eta_event_l <- rnorm(n_events, 0, gmm_tau)

  psi_virgin <- matrix(0, n_buildings, M_DRAWS)
  psi_predamage <- matrix(buildings_df$initial_psi, n_buildings, M_DRAWS)

  max_pgv_l <- rep(0, n_buildings)
  pgv_history <- matrix(0, n_buildings, n_events)  # Track PGV history
  n_record_l <- rep(0, n_buildings)
  n_similar_l <- rep(0, n_buildings)

  psi0_virgin <- rep(0, n_buildings)
  psi0_predamage <- buildings_df$initial_psi

  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]

    eta_within <- rnorm(n_buildings, 0, gmm_phiSS)
    lnPGV_realized <- lnPGV_median_all[, j] + eta_event_l[j] + eta_site_l + eta_within
    pgv_event <- exp(lnPGV_realized) * 10

    # Compute N_eq: count prior events "similar" to current
    if (j_idx > 1) {
      prior_pgvs <- pgv_history[, 1:(j_idx-1), drop = FALSE]
      similar <- (prior_pgvs > (1 - RHO) * pgv_event) &
                 (prior_pgvs < (1 + RHO) * pgv_event) &
                 (prior_pgvs > 0)
      n_eq <- rowSums(similar)
    } else {
      n_eq <- rep(0, n_buildings)
    }

    # Determine regime
    is_record <- pgv_event > max_pgv_l & pgv_event >= PGV_THRESHOLD
    is_similar <- n_eq > 0 & !is_record & pgv_event >= PGV_THRESHOLD

    # =========================================================================
    # REGIME 1: Record-breaking events
    # ΔΨ = Surrogate(Ψ, V_new, N=2) - Surrogate(Ψ, V_max, N=1)
    # =========================================================================
    if (any(is_record)) {
      idx <- which(is_record)
      pgv_new <- pgv_event[idx]
      pgv_old <- max_pgv_l[idx]

      # Virgin walls
      inc_new_v <- predict_damage_mean_with_n(pgv_new, psi0_virgin[idx], rep(2, length(idx)))
      inc_old_v <- predict_damage_mean_with_n(pgv_old, psi0_virgin[idx], rep(1, length(idx)))
      delta_v <- pmax(inc_new_v - inc_old_v, 0)
      psi_virgin[idx, ] <- pmin(psi_virgin[idx, ] + delta_v, PSI_MAX)

      # Pre-damage
      inc_new_p <- predict_damage_mean_with_n(pgv_new, psi0_predamage[idx], rep(2, length(idx)))
      inc_old_p <- predict_damage_mean_with_n(pgv_old, psi0_predamage[idx], rep(1, length(idx)))
      delta_p <- pmax(inc_new_p - inc_old_p, 0)
      psi_predamage[idx, ] <- pmin(psi_predamage[idx, ] + delta_p, PSI_MAX)

      max_pgv_l[idx] <- pgv_new
      n_record_l[idx] <- n_record_l[idx] + 1
    }

    # =========================================================================
    # REGIME 2: Similar events (small marginal N contribution)
    # ΔΨ = Surrogate(Ψ, V, N_eq+1) - Surrogate(Ψ, V, N_eq)
    # =========================================================================
    if (any(is_similar)) {
      idx <- which(is_similar)
      pgv_sim <- pgv_event[idx]
      n_eq_sim <- n_eq[idx]

      # Virgin walls
      inc_new_v <- predict_damage_mean_with_n(pgv_sim, psi0_virgin[idx], n_eq_sim + 1)
      inc_old_v <- predict_damage_mean_with_n(pgv_sim, psi0_virgin[idx], n_eq_sim)
      delta_v <- pmax(inc_new_v - inc_old_v, 0)
      psi_virgin[idx, ] <- pmin(psi_virgin[idx, ] + delta_v, PSI_MAX)

      # Pre-damage
      inc_new_p <- predict_damage_mean_with_n(pgv_sim, psi0_predamage[idx], n_eq_sim + 1)
      inc_old_p <- predict_damage_mean_with_n(pgv_sim, psi0_predamage[idx], n_eq_sim)
      delta_p <- pmax(inc_new_p - inc_old_p, 0)
      psi_predamage[idx, ] <- pmin(psi_predamage[idx, ] + delta_p, PSI_MAX)

      n_similar_l[idx] <- n_similar_l[idx] + 1
    }

    # Store PGV in history
    pgv_history[, j_idx] <- pgv_event
  }

  all_psi_virgin[[l]] <- psi_virgin
  all_psi_predamage[[l]] <- psi_predamage
  max_pgv_all[, l] <- max_pgv_l
  n_record_all[, l] <- n_record_l
  n_similar_all[, l] <- n_similar_l

  l_elapsed <- as.numeric(difftime(Sys.time(), l_start, units = "mins"))
  cat(sprintf("  %.1f min | record: %.1f | similar: %.1f events/building\n",
              l_elapsed, mean(n_record_l), mean(n_similar_l)))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\nComplete! Total: %.1f min\n\n", total_elapsed))

# =============================================================================
# 6. COMPUTE STATISTICS
# =============================================================================

cat("Step 5: Computing statistics...\n")

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
    n_record_mean = rowMeans(n_record_all),
    n_similar_mean = rowMeans(n_similar_all)
  )

# =============================================================================
# 7. RESULTS
# =============================================================================

cat("\n=== DAF MODEL RESULTS ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("ρ (similarity):", RHO, "\n\n")

cat("--- VIRGIN WALLS ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_virgin_mean), 4), "\n")
cat("P(visible):", round(100 * mean(buildings_df$p_visible_virgin), 2), "%\n")
cat("P(moderate):", round(100 * mean(buildings_df$p_moderate_virgin), 2), "%\n\n")

cat("--- EVENT STATISTICS ---\n")
cat("Mean record-breaking events:", round(mean(buildings_df$n_record_mean), 2), "\n")
cat("Mean similar events:", round(mean(buildings_df$n_similar_mean), 2), "\n\n")

by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    psi = mean(psi_virgin_mean),
    p_vis = mean(p_visible_virgin),
    n_rec = mean(n_record_mean),
    n_sim = mean(n_similar_mean),
    .groups = "drop"
  )

cat("By zone:\n")
print(by_zone)

# =============================================================================
# 8. COMPARISON
# =============================================================================

cat("\n=== COMPARISON ===\n")

ratchet_file <- here("outputs", "models", "spatial_damage_ratchet.rds")
if (file.exists(ratchet_file)) {
  ratchet <- readRDS(ratchet_file)
  cat("DAF vs Hard Ratchet:\n")
  cat("  DAF mean Ψ:", round(mean(buildings_df$psi_virgin_mean), 4), "\n")
  cat("  Ratchet mean Ψ:", round(mean(ratchet$psi_virgin_mean), 4), "\n")
  cat("  Ratio (DAF/Ratchet):", round(mean(buildings_df$psi_virgin_mean) / mean(ratchet$psi_virgin_mean), 3), "\n")
}

# =============================================================================
# 9. SAVE
# =============================================================================

cat("\nStep 6: Saving results...\n")

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
    max_pgv, max_pgv_mean,
    n_record_mean, n_similar_mean
  )

saveRDS(results, here("outputs", "models", "spatial_damage_daf.rds"))
cat("Saved: outputs/models/spatial_damage_daf.rds\n")

cat("\n=== DONE ===\n")
