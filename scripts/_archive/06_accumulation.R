# =============================================================================
# 06_accumulation.R
#
# Cumulative damage intensity measure for arbitrary earthquake sequences.
# Bridges from FEM's "N identical events" to real earthquake catalogues.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(mgcv)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load data and models
# -----------------------------------------------------------------------------

fem <- readRDS(here::here("outputs", "models", "fem_prepared.rds"))
m1_init <- readRDS(here::here("outputs", "models", "baseline_stage1_init.rds"))
m2_severity <- readRDS(here::here("outputs", "models", "baseline_stage2_severity.rds"))

fig_dir <- here::here("outputs", "figures")

# -----------------------------------------------------------------------------
# 2. Cumulative Intensity Measure Theory
# -----------------------------------------------------------------------------
#
# The FEM data uses N repetitions of identical ground motion.
# For arbitrary sequences, we define a cumulative intensity:
#
#   S = sum_t (PGV_t / PGV_ref)^alpha
#
# With N identical events at PGV:
#   S = N * (PGV / PGV_ref)^alpha
#
# Inverting: if we know alpha, we can map any sequence to equivalent S.
#
# We estimate alpha by finding the exponent that best collapses the
# (PGV, N) -> DeltaPsi relationship onto a single S -> DeltaPsi curve.
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 3. Estimate cumulative intensity exponent (alpha)
# -----------------------------------------------------------------------------

cat("\n=== Estimating Cumulative Intensity Exponent ===\n")

# Reference PGV (arbitrary choice, affects intercept not shape)
PGV_ref <- 10  # mm/s

# Function to compute cumulative intensity
compute_S <- function(pgv, n, alpha, pgv_ref = PGV_ref) {
  n * (pgv / pgv_ref)^alpha
}

# Grid search for optimal alpha
alpha_grid <- seq(0.5, 3.0, by = 0.1)
fit_results <- list()

for (alpha in alpha_grid) {
  fem_alpha <- fem |>
    mutate(
      S = compute_S(PGV, N, alpha),
      log_S = log(S + 0.01)  # Small offset to handle S near 0
    )

  # Fit model: DeltaPsi ~ S
  # Using subset with damage for cleaner signal
  fem_pos <- fem_alpha |> filter(has_damage == 1)

  m_alpha <- tryCatch({
    gam(delta_psi_pos ~ s(log_S, k = 5) + InitialPsi + Material +
                        FacadeType + SoilProfile + EarthquakeType,
        family = Gamma(link = "log"),
        data = fem_pos,
        method = "REML")
  }, error = function(e) NULL)

  if (!is.null(m_alpha)) {
    fit_results[[as.character(alpha)]] <- list(
      alpha = alpha,
      deviance = deviance(m_alpha),
      aic = AIC(m_alpha),
      gcv = m_alpha$gcv.ubre
    )
  }
}

# Find best alpha
results_df <- bind_rows(fit_results)
best_alpha <- results_df$alpha[which.min(results_df$aic)]

cat("Best alpha (by AIC):", best_alpha, "\n")
cat("AIC at best alpha:", min(results_df$aic), "\n")

# Plot alpha search
p_alpha <- results_df |>
  ggplot(aes(x = alpha, y = aic)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = best_alpha, color = "red", linetype = "dashed") +
  labs(
    title = "Cumulative Intensity Exponent Selection",
    subtitle = sprintf("Optimal alpha = %.2f", best_alpha),
    x = expression(alpha),
    y = "AIC"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "accumulation_alpha_search.png"), p_alpha,
       width = 8, height = 5, dpi = 150)
cat("Saved: accumulation_alpha_search.png\n")

# -----------------------------------------------------------------------------
# 4. Fit model with cumulative intensity S
# -----------------------------------------------------------------------------

cat("\n=== Fitting Cumulative Intensity Model ===\n")

# Add S to data with optimal alpha
fem_S <- fem |>
  mutate(
    S = compute_S(PGV, N, best_alpha),
    log_S = log(S + 0.01)
  )

# Stage 1: P(damage > 0) ~ S
m1_S <- gam(
  has_damage ~ s(log_S, k = 6, bs = "tp") +
               s(InitialPsi, k = 5, bs = "tp") +
               Material +
               FacadeType +
               SoilProfile +
               EarthquakeType,
  family = binomial(link = "logit"),
  data = fem_S,
  method = "REML"
)

cat("\nStage 1 (Cumulative S) Summary:\n")
print(summary(m1_S))

# Stage 2: E[damage | damage > 0] ~ S
fem_S_pos <- fem_S |> filter(has_damage == 1)

m2_S <- gam(
  delta_psi_pos ~ s(log_S, k = 6, bs = "tp") +
                  s(InitialPsi, k = 5, bs = "tp") +
                  Material +
                  FacadeType +
                  SoilProfile +
                  EarthquakeType,
  family = Gamma(link = "log"),
  data = fem_S_pos,
  method = "REML"
)

cat("\nStage 2 (Cumulative S) Summary:\n")
print(summary(m2_S))

# -----------------------------------------------------------------------------
# 5. Visualize S -> DeltaPsi relationship
# -----------------------------------------------------------------------------

cat("\n=== Visualizing Cumulative Intensity Relationship ===\n")

# Plot raw data colored by N
p_S_raw <- fem_S |>
  ggplot(aes(x = S, y = delta_psi_pos, color = factor(N))) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +
  scale_color_viridis_d(option = "D") +
  labs(
    title = sprintf("Damage vs Cumulative Intensity (alpha = %.2f)", best_alpha),
    subtitle = "Points from different N should collapse onto single curve",
    x = expression(S == N %.% (PGV/PGV[ref])^alpha),
    y = expression(Delta * Psi),
    color = "N"
  ) +
  theme_minimal()

# Fitted curve
S_grid <- tibble(
  log_S = seq(min(fem_S$log_S), max(fem_S$log_S), length.out = 100),
  InitialPsi = 0,
  Material = 1.0,
  FacadeType = "A",
  SoilProfile = "A",
  EarthquakeType = "ZN"
) |>
  mutate(S = exp(log_S))

S_grid$pred_p <- predict(m1_S, newdata = S_grid, type = "response")
S_grid$pred_mu <- predict(m2_S, newdata = S_grid, type = "response")
S_grid$pred_E <- S_grid$pred_p * S_grid$pred_mu

p_S_fit <- S_grid |>
  ggplot(aes(x = S, y = pred_E)) +
  geom_line(color = "red", linewidth = 1.5) +
  geom_point(data = fem_S, aes(y = delta_psi_pos), alpha = 0.1) +
  scale_x_log10() +
  labs(
    title = "Fitted Damage vs Cumulative Intensity",
    subtitle = "Red line: E[DeltaPsi] from hurdle model",
    x = expression(S == N %.% (PGV/PGV[ref])^alpha),
    y = expression(E[Delta * Psi])
  ) +
  theme_minimal()

combined_S <- p_S_raw / p_S_fit
ggsave(file.path(fig_dir, "accumulation_S_relationship.png"), combined_S,
       width = 10, height = 10, dpi = 150)
cat("Saved: accumulation_S_relationship.png\n")

# -----------------------------------------------------------------------------
# 6. Compare S-based model to (PGV, N) model
# -----------------------------------------------------------------------------

cat("\n=== Model Comparison: S vs (PGV, N) ===\n")

# Get predictions from both models on same data
fem_compare <- fem_S |>
  filter(has_damage == 1) |>
  mutate(
    # S-based prediction
    pred_S = predict(m2_S, newdata = cur_data(), type = "response"),
    # Original (PGV, N) prediction
    pred_PN = predict(m2_severity, newdata = cur_data(), type = "response"),
    # Residuals
    resid_S = delta_psi_pos - pred_S,
    resid_PN = delta_psi_pos - pred_PN
  )

cat("S-based model RMSE:", sqrt(mean(fem_compare$resid_S^2)), "\n")
cat("(PGV,N)-based model RMSE:", sqrt(mean(fem_compare$resid_PN^2)), "\n")
cat("S-based model R2:", 1 - sum(fem_compare$resid_S^2)/sum((fem_compare$delta_psi_pos - mean(fem_compare$delta_psi_pos))^2), "\n")
cat("(PGV,N)-based model R2:", 1 - sum(fem_compare$resid_PN^2)/sum((fem_compare$delta_psi_pos - mean(fem_compare$delta_psi_pos))^2), "\n")

# -----------------------------------------------------------------------------
# 7. Functions for arbitrary event sequences
# -----------------------------------------------------------------------------

#' Calculate cumulative intensity for a sequence of earthquakes
#'
#' @param pgv_sequence Vector of PGV values (mm/s) for each event
#' @param alpha Cumulative intensity exponent (default: estimated)
#' @param pgv_ref Reference PGV (default: 10 mm/s)
#' @return Cumulative intensity S after all events
calculate_cumulative_S <- function(pgv_sequence, alpha = best_alpha, pgv_ref = PGV_ref) {
  sum((pgv_sequence / pgv_ref)^alpha)
}

#' Simulate damage accumulation over earthquake sequence
#'
#' @param pgv_sequence Vector of PGV values for each event
#' @param initial_psi Initial damage state
#' @param material Material strength (0.7, 1.0, or 1.3)
#' @param facade_type "A" or "B"
#' @param soil_profile "A" or "B"
#' @param n_sims Number of Monte Carlo simulations
#' @param m1 Stage 1 model (hurdle)
#' @param m2 Stage 2 model (severity)
#' @return List with simulated final damage states
simulate_damage_sequence <- function(
  pgv_sequence,
  initial_psi = 0,
  material = 1.0,
  facade_type = "A",
  soil_profile = "A",
  n_sims = 1000,
  m1 = m1_S,
  m2 = m2_S,
  alpha = best_alpha
) {
  n_events <- length(pgv_sequence)
  psi_trajectories <- matrix(NA, nrow = n_sims, ncol = n_events + 1)
  psi_trajectories[, 1] <- initial_psi

  # Gamma shape from model
  gamma_shape <- 1 / m2$scale

  for (t in seq_along(pgv_sequence)) {
    current_psi <- psi_trajectories[, t]

    # Calculate cumulative S up to and including this event
    S_cumulative <- calculate_cumulative_S(pgv_sequence[1:t], alpha)

    # Prepare prediction data
    pred_data <- tibble(
      log_S = rep(log(S_cumulative + 0.01), n_sims),
      InitialPsi = current_psi,
      Material = material,
      FacadeType = facade_type,
      SoilProfile = soil_profile,
      EarthquakeType = "ZN"  # Placeholder (will integrate over)
    )

    # Predict P(damage) and E[damage|damage>0]
    p_damage <- predict(m1, newdata = pred_data, type = "response")
    mu_damage <- predict(m2, newdata = pred_data, type = "response")

    # Simulate damage increment
    damage_occurs <- rbinom(n_sims, 1, p_damage)
    damage_amount <- rgamma(n_sims, shape = gamma_shape, rate = gamma_shape / mu_damage)
    delta_psi <- damage_occurs * damage_amount

    # Update cumulative damage
    psi_trajectories[, t + 1] <- current_psi + delta_psi
  }

  return(list(
    trajectories = psi_trajectories,
    final_psi = psi_trajectories[, n_events + 1],
    pgv_sequence = pgv_sequence
  ))
}

# -----------------------------------------------------------------------------
# 8. Example: Simulate a realistic earthquake sequence
# -----------------------------------------------------------------------------

cat("\n=== Example Simulation ===\n")

# Example: 10 events with varying PGV
set.seed(42)
example_sequence <- c(3, 5, 8, 4, 12, 6, 15, 3, 7, 10)  # mm/s

result <- simulate_damage_sequence(
  pgv_sequence = example_sequence,
  initial_psi = 0,
  material = 1.0,
  n_sims = 1000
)

cat("Example sequence PGV:", example_sequence, "\n")
cat("Cumulative S:", calculate_cumulative_S(example_sequence), "\n")
cat("Final Psi - Mean:", mean(result$final_psi), "\n")
cat("Final Psi - Median:", median(result$final_psi), "\n")
cat("Final Psi - 95% CI:", quantile(result$final_psi, c(0.025, 0.975)), "\n")
cat("P(Psi >= 1):", mean(result$final_psi >= 1), "\n")

# Plot trajectories
traj_df <- as_tibble(result$trajectories) |>
  mutate(sim = row_number()) |>
  pivot_longer(cols = -sim, names_to = "event", values_to = "psi") |>
  mutate(event = as.integer(gsub("V", "", event)) - 1)

p_traj <- traj_df |>
  filter(sim <= 100) |>  # Subset for clarity
  ggplot(aes(x = event, y = psi, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed",
             linewidth = 1) +
  annotate("text", x = 8, y = 1.1, label = "Visible damage (Psi = 1)",
           color = "red", hjust = 0) +
  labs(
    title = "Simulated Damage Trajectories",
    subtitle = paste("PGV sequence:", paste(example_sequence, collapse = ", ")),
    x = "Event Number",
    y = expression(Psi)
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "accumulation_example_trajectories.png"), p_traj,
       width = 10, height = 6, dpi = 150)
cat("Saved: accumulation_example_trajectories.png\n")

# -----------------------------------------------------------------------------
# 9. Save models and functions
# -----------------------------------------------------------------------------

# Save models
saveRDS(m1_S, here::here("outputs", "models", "accumulation_stage1_S.rds"))
saveRDS(m2_S, here::here("outputs", "models", "accumulation_stage2_S.rds"))

# Save parameters
accum_params <- list(
  alpha = best_alpha,
  pgv_ref = PGV_ref,
  gamma_shape = 1 / m2_S$scale
)
saveRDS(accum_params, here::here("outputs", "models", "accumulation_params.rds"))

# Save functions
save(
  calculate_cumulative_S,
  simulate_damage_sequence,
  file = here::here("outputs", "models", "accumulation_functions.RData")
)

cat("\nAccumulation models and functions saved.\n")
cat("Key parameters:\n")
cat("  alpha =", best_alpha, "\n")
cat("  PGV_ref =", PGV_ref, "mm/s\n")

cat("\nAccumulation analysis complete.\n")
