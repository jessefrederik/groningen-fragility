# =============================================================================
# 16_korswagen_daf.R
#
# Korswagen lognormal fragility with DAF (Damage Accumulation Function) logic:
# - Record-breaking events (V > V_max): Full damage from lognormal fragility
# - Similar events (V ~ prior): Attenuated - marginal damage approaches zero
#
# Key insight from Korswagen's Monte Carlo: If a building experiences the same
# PGV multiple times, most damage occurs on the FIRST event. Subsequent similar
# events cause negligible additional damage (the structure has already "yielded"
# at that intensity).
#
# Implementation:
# 1. Track max PGV experienced (V_max) and similar event count (N_similar)
# 2. Record-breaking: Sample Δψ from Korswagen fragility (full damage)
# 3. Similar events: Δψ ≈ 0 (already exercised at this intensity)
#
# @author Generated with Claude Code
# @date 2026-01-05
# =============================================================================

library(tidyverse)
library(sf)
library(here)
library(matrixStats)

source(here("helperfuncties", "fetch_knmi_events.R"))
source(here("helperfuncties", "pgv_model_original.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

PSI_MAX <- 3.0
MAG_THRESHOLD <- 2.0
PGV_THRESHOLD <- 2.0    # mm/s
L_GMM <- 10             # GMM realizations
RHO <- 0.25             # Similarity tolerance (Korswagen Eq. 6.7)
SAMPLE_FRAC <- 1.0      # Full dataset

cat("=== Korswagen Lognormal with DAF Logic ===\n")
cat("Configuration:\n")
cat("  GMM realizations:", L_GMM, "\n")
cat("  Similarity tolerance (ρ):", RHO, "\n")
cat("  PGV threshold:", PGV_THRESHOLD, "mm/s\n\n")

# =============================================================================
# KORSWAGEN TABLE 10: LOGNORMAL FRAGILITY PARAMETERS
# =============================================================================

korswagen_params <- tribble(
  ~psi_0, ~psi_final, ~mu, ~sigma,
  # Psi_0 = 0
  0.0, 0.5, 2.604, 0.624,
  0.0, 1.0, 3.370, 0.725,
  0.0, 1.5, 4.323, 0.860,
  0.0, 2.0, 5.518, 1.082,
  0.0, 2.5, 7.416, 1.516,
  0.0, 3.0, 8.869, 1.729,
  # Psi_0 = 0.5
  0.5, 1.0, 3.002, 0.931,
  0.5, 1.5, 4.070, 0.983,
  0.5, 2.0, 5.523, 1.259,
  0.5, 2.5, 7.758, 1.780,
  0.5, 3.0, 11.287, 2.606,
  # Psi_0 = 1.0
  1.0, 1.5, 3.411, 1.034,
  1.0, 2.0, 5.293, 1.444,
  1.0, 2.5, 7.955, 2.091,
  1.0, 3.0, 11.913, 3.034,
  # Psi_0 = 1.5
  1.5, 2.0, 4.631, 1.860,
  1.5, 2.5, 7.929, 2.585,
  1.5, 3.0, 12.363, 3.610
)

cat("Korswagen Table 10 loaded\n")
cat("Median PGV for P(Psi>=0.5|Psi_0=0) =", round(exp(2.604), 1), "mm/s\n")
cat("Median PGV for P(Psi>=1.0|Psi_0=0) =", round(exp(3.370), 1), "mm/s\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

get_discrete_psi0 <- function(psi) {
  psi_capped <- pmin(psi, 1.5)
  floor(psi_capped * 2) / 2
}

#' Sample Delta-Psi from Korswagen lognormal (for RECORD-BREAKING events)
sample_delta_psi_korswagen <- function(pgv_mm_s, current_psi, params = korswagen_params) {
  psi_0_discrete <- get_discrete_psi0(current_psi)

  available <- params |>
    filter(psi_0 == psi_0_discrete, psi_final > current_psi) |>
    arrange(psi_final)

  if (nrow(available) == 0) return(0)

  ln_pgv <- log(pgv_mm_s)
  p_exceed <- pnorm((ln_pgv - available$mu) / available$sigma)

  u <- runif(1)
  exceeded_idx <- which(p_exceed >= u)

  if (length(exceeded_idx) == 0) return(0)

  psi_final <- available$psi_final[max(exceeded_idx)]
  max(0, psi_final - current_psi)
}

#' Vectorized sampling for record-breaking events only
sample_delta_psi_record_breaking <- function(pgv_vec, psi_vec, params = korswagen_params) {
  n <- length(pgv_vec)
  delta <- numeric(n)

  for (i in seq_len(n)) {
    delta[i] <- sample_delta_psi_korswagen(pgv_vec[i], psi_vec[i], params)
  }

  delta
}

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
  buildings <- slice_sample(buildings, n = round(nrow(buildings) * SAMPLE_FRAC))
}

n_buildings <- nrow(buildings)
cat("  Buildings:", format(n_buildings, big.mark = ","), "\n")

buildings <- buildings |>
  mutate(
    bouwjaar_numeric = as.numeric(bouwjaar),
    zone = case_when(
      dist_to_center < 10 ~ "Near (<10km)",
      dist_to_center < 25 ~ "Mid (10-25km)",
      TRUE ~ "Far (>25km)"
    )
  )

buildings_df <- st_drop_geometry(buildings)
cat("  Data ready\n\n")

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
# 3. PRE-COMPUTE MEDIAN PGV
# =============================================================================

cat("Step 3: Pre-computing median PGV...\n")

event_order <- order(groningen_eq$time_utc)
eq_lat <- groningen_eq$lat
eq_lon <- groningen_eq$lon
eq_depth <- groningen_eq$depth_km
eq_mag <- groningen_eq$mag

lnPGV_median_all <- matrix(NA, n_buildings, n_events)
for (j in 1:n_events) {
  epi_dist <- calc_dist_km(buildings_df$lat, buildings_df$lon, eq_lat[j], eq_lon[j])
  hypo_dist <- sqrt(epi_dist^2 + eq_depth[j]^2)
  pgv_result <- pgv_model(M = eq_mag[j], Rhyp = hypo_dist, comp = "LRG", VS30 = buildings_df$vs30)
  lnPGV_median_all[, j] <- pgv_result$lnPGV
}
cat("  Done\n\n")

# =============================================================================
# 4. MAIN SIMULATION WITH DAF LOGIC
# =============================================================================

cat("Step 4: Running DAF simulation...\n")
cat("  Key: Only RECORD-BREAKING events cause damage\n")
cat("  Similar events (within ρ =", RHO, ") are attenuated to zero\n\n")

psi_final_all <- matrix(NA, n_buildings, L_GMM)
max_pgv_all <- matrix(NA, n_buildings, L_GMM)
n_record_all <- matrix(0, n_buildings, L_GMM)
n_similar_all <- matrix(0, n_buildings, L_GMM)

set.seed(42)
global_start <- Sys.time()

for (l in 1:L_GMM) {
  cat(sprintf("GMM realization %d/%d", l, L_GMM))
  l_start <- Sys.time()

  # GMM uncertainty
  eta_site <- rnorm(n_buildings, 0, gmm_phiS2S)
  eta_event <- rnorm(n_events, 0, gmm_tau)

  # Initialize
  psi <- rep(0, n_buildings)      # All buildings start undamaged
  max_pgv <- rep(0, n_buildings)  # Track max PGV seen
  n_record <- rep(0, n_buildings)
  n_similar <- rep(0, n_buildings)

  # Track PGV history for similarity detection
  pgv_history <- matrix(0, n_buildings, n_events)

  for (j_idx in seq_along(event_order)) {
    j <- event_order[j_idx]

    # Realized PGV with uncertainty
    eta_within <- rnorm(n_buildings, 0, gmm_phiSS)
    lnPGV_realized <- lnPGV_median_all[, j] + eta_event[j] + eta_site + eta_within
    pgv_event <- exp(lnPGV_realized) * 10  # cm/s to mm/s

    # Store in history
    pgv_history[, j_idx] <- pgv_event

    # Classify events
    above_threshold <- pgv_event >= PGV_THRESHOLD
    is_record <- above_threshold & (pgv_event > max_pgv * (1 + RHO))

    # Count similar events (within ρ tolerance of current)
    if (j_idx > 1) {
      prior_pgvs <- pgv_history[, 1:(j_idx-1), drop = FALSE]
      similar <- (prior_pgvs > (1 - RHO) * pgv_event) &
                 (prior_pgvs < (1 + RHO) * pgv_event) &
                 (prior_pgvs > 0)
      n_prior_similar <- rowSums(similar)
    } else {
      n_prior_similar <- rep(0, n_buildings)
    }

    is_similar <- above_threshold & !is_record & (n_prior_similar > 0)

    # =========================================================================
    # REGIME 1: RECORD-BREAKING EVENTS
    # Full damage from Korswagen fragility
    # =========================================================================
    if (any(is_record)) {
      idx <- which(is_record)

      delta <- sample_delta_psi_record_breaking(
        pgv_vec = pgv_event[idx],
        psi_vec = psi[idx]
      )

      psi[idx] <- pmin(psi[idx] + delta, PSI_MAX)
      max_pgv[idx] <- pgv_event[idx]
      n_record[idx] <- n_record[idx] + 1
    }

    # =========================================================================
    # REGIME 2: SIMILAR EVENTS
    # NO damage - structure already "exercised" at this intensity
    # (This is the key DAF insight: marginal damage → 0 for repeated similar PGV)
    # =========================================================================
    if (any(is_similar)) {
      idx <- which(is_similar)
      n_similar[idx] <- n_similar[idx] + 1
      # Delta = 0, so psi unchanged
    }

    # Note: Events below threshold or between record and similar → no damage
  }

  psi_final_all[, l] <- psi
  max_pgv_all[, l] <- max_pgv
  n_record_all[, l] <- n_record
  n_similar_all[, l] <- n_similar

  l_elapsed <- as.numeric(difftime(Sys.time(), l_start, units = "secs"))
  cat(sprintf(" - %.1fs | record: %.2f | similar: %.2f\n",
              l_elapsed, mean(n_record), mean(n_similar)))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat(sprintf("\nComplete! Total: %.1f min\n\n", total_elapsed))

# =============================================================================
# 5. COMPUTE STATISTICS
# =============================================================================

cat("Step 5: Computing statistics...\n")

buildings_df <- buildings_df |>
  mutate(
    psi_mean = rowMeans(psi_final_all),
    psi_sd = rowSds(psi_final_all),
    psi_p10 = rowQuantiles(psi_final_all, probs = 0.10),
    psi_p50 = rowQuantiles(psi_final_all, probs = 0.50),
    psi_p90 = rowQuantiles(psi_final_all, probs = 0.90),

    p_visible = rowMeans(psi_final_all >= 1.0),
    p_moderate = rowMeans(psi_final_all >= 2.0),
    p_severe = rowMeans(psi_final_all >= 3.0),

    max_pgv = rowMaxs(max_pgv_all),
    max_pgv_mean = rowMeans(max_pgv_all),
    n_record = rowMeans(n_record_all),
    n_similar = rowMeans(n_similar_all)
  )

# =============================================================================
# 6. RESULTS SUMMARY
# =============================================================================

cat("\n=== KORSWAGEN DAF RESULTS ===\n")
cat("Buildings:", format(n_buildings, big.mark = ","), "\n")
cat("Earthquakes:", n_events, "\n")
cat("Similarity tolerance (ρ):", RHO, "\n\n")

cat("--- DAMAGE STATISTICS ---\n")
cat("Mean Ψ:", round(mean(buildings_df$psi_mean), 4), "\n")
cat("  p10:", round(mean(buildings_df$psi_p10), 4),
    "| p50:", round(mean(buildings_df$psi_p50), 4),
    "| p90:", round(mean(buildings_df$psi_p90), 4), "\n\n")

cat("--- EXCEEDANCE PROBABILITIES ---\n")
cat("P(Ψ >= 0.5):", round(100 * mean(buildings_df$psi_mean >= 0.5), 2), "%\n")
cat("P(Ψ >= 1) visible:", round(100 * mean(buildings_df$p_visible), 2), "%\n")
cat("P(Ψ >= 2) moderate:", round(100 * mean(buildings_df$p_moderate), 2), "%\n\n")

cat("--- EVENT BREAKDOWN ---\n")
cat("Mean record-breaking events:", round(mean(buildings_df$n_record), 2), "\n")
cat("Mean similar events (attenuated):", round(mean(buildings_df$n_similar), 2), "\n\n")

cat("--- BY DISTANCE ZONE ---\n")
by_zone <- buildings_df |>
  group_by(zone) |>
  summarise(
    n = n(),
    psi_mean = mean(psi_mean),
    p_visible = mean(p_visible),
    n_record = mean(n_record),
    n_similar = mean(n_similar),
    max_pgv = mean(max_pgv_mean),
    .groups = "drop"
  )
print(by_zone)

# =============================================================================
# 7. COMPARISON WITH ADDITIVE MODEL
# =============================================================================

cat("\n=== COMPARISON ===\n")

# Load additive results if available
additive_file <- here("outputs", "models", "korswagen_lognormal_results.rds")
if (file.exists(additive_file)) {
  additive <- readRDS(additive_file)

  # Match buildings
  common_ids <- intersect(buildings_df$identificatie, additive$identificatie)
  if (length(common_ids) > 1000) {
    daf_sub <- buildings_df |> filter(identificatie %in% common_ids)
    add_sub <- additive |> filter(identificatie %in% common_ids)

    # Sort to align
    daf_sub <- daf_sub |> arrange(identificatie)
    add_sub <- add_sub |> arrange(identificatie)

    cat("DAF vs Additive (same buildings, n =", length(common_ids), "):\n")
    cat("  DAF mean Ψ:", round(mean(daf_sub$psi_mean), 4), "\n")
    cat("  Additive mean Ψ:", round(mean(add_sub$psi_mean), 4), "\n")
    cat("  Ratio (DAF/Additive):", round(mean(daf_sub$psi_mean) / mean(add_sub$psi_mean), 3), "\n")
    cat("  DAF P(visible):", round(100 * mean(daf_sub$p_visible), 2), "%\n")
    cat("  Additive P(visible):", round(100 * mean(add_sub$p_visible), 2), "%\n")
  }
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

cat("\nStep 6: Saving results...\n")

results <- buildings_df |>
  select(
    identificatie, bouwjaar, bouwjaar_numeric, postcode, pc4, vs30,
    lat, lon, dist_to_center, zone,
    psi_mean, psi_sd, psi_p10, psi_p50, psi_p90,
    p_visible, p_moderate, p_severe,
    max_pgv, max_pgv_mean, n_record, n_similar
  )

output_file <- here("outputs", "models", "korswagen_daf_results.rds")
saveRDS(results, output_file)
cat("Saved:", output_file, "\n")

# Summary
summary_list <- list(
  by_zone = by_zone,
  config = list(
    rho = RHO,
    l_gmm = L_GMM,
    n_buildings = n_buildings,
    n_events = n_events
  ),
  overall = tibble(
    method = "korswagen_daf",
    psi_mean = mean(buildings_df$psi_mean),
    p_visible = mean(buildings_df$p_visible),
    p_moderate = mean(buildings_df$p_moderate),
    n_record_mean = mean(buildings_df$n_record),
    n_similar_mean = mean(buildings_df$n_similar)
  )
)

saveRDS(summary_list, here("outputs", "models", "korswagen_daf_summaries.rds"))

cat("\n=== DONE ===\n")
cat("Key insight: Only", round(mean(buildings_df$n_record), 2),
    "record-breaking events/building cause damage\n")
cat("The", round(mean(buildings_df$n_similar), 2),
    "similar events/building are attenuated to zero\n")
