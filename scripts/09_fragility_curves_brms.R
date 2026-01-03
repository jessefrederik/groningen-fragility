# =============================================================================
# 09_fragility_curves_brms.R
#
# Derive fragility curves P(Psi >= threshold | PGV) using brms model
# with full posterior uncertainty propagation.
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(brms)
library(patchwork)

# Source helper functions
source(here::here("scripts", "helpers", "brms_predict.R"))

# Load brms model
fit_brms <- readRDS(here::here("outputs", "models", "brms_hurdle_gamma_A.rds"))

fig_dir <- here::here("outputs", "figures")

cat("Loaded brms model\n")

# -----------------------------------------------------------------------------
# 1. Generate fragility curves for single event
# -----------------------------------------------------------------------------

cat("\n=== Generating Single-Event Fragility Curves ===\n")

# PGV range for fragility curves (mm/s)
pgv_seq <- c(2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 96, 112, 128)

# Scenarios
scenarios <- expand_grid(
  pgv = pgv_seq,
  initial_psi = c(0, 0.5, 1.0),
  material = 1.0,
  facade_type = "A",
  soil_profile = "A"
)

# Get predictions with full posterior uncertainty
cat("Computing predictions for", nrow(scenarios), "scenario combinations...\n")

n_draws <- 1000
results_list <- list()

pb <- txtProgressBar(min = 0, max = nrow(scenarios), style = 3)
for (i in 1:nrow(scenarios)) {
  sc <- scenarios[i, ]

  pred <- predict_damage_brms(
    fit_brms = fit_brms,
    pgv = sc$pgv,
    n = 1,
    initial_psi = sc$initial_psi,
    material = sc$material,
    facade_type = sc$facade_type,
    soil_profile = sc$soil_profile,
    ndraws = n_draws,
    summary = FALSE
  )

  # Calculate exceedance probabilities
  results_list[[i]] <- tibble(
    pgv = sc$pgv,
    initial_psi = sc$initial_psi,
    material = sc$material,
    facade_type = sc$facade_type,
    soil_profile = sc$soil_profile,
    p_damage = mean(pred > 0),
    p_visible = mean(pred >= 1),
    p_moderate = mean(pred >= 2),
    mean_delta_psi = mean(pred),
    sd_delta_psi = sd(pred),
    q10_delta_psi = quantile(pred, 0.10),
    q50_delta_psi = quantile(pred, 0.50),
    q90_delta_psi = quantile(pred, 0.90)
  )

  setTxtProgressBar(pb, i)
}
close(pb)

fragility_single <- bind_rows(results_list)

cat("Single-event fragility curves computed.\n")

# -----------------------------------------------------------------------------
# 2. Visualize single-event fragility curves
# -----------------------------------------------------------------------------

cat("\n=== Visualizing Fragility Curves ===\n")

# P(Psi >= 1 | PGV) by initial damage
p_frag_visible <- fragility_single |>
  mutate(initial_label = paste0("Initial Psi = ", initial_psi)) |>
  ggplot(aes(x = pgv, y = p_visible, color = initial_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 120, y = 0.12, label = "10% exceedance", color = "gray50", hjust = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Fragility Curve: P(Visible Damage | PGV) - brms Model",
    subtitle = "Single event, Material = 1.0, FacadeType = A, SoilProfile = A",
    x = "PGV (mm/s)",
    y = expression(P(Delta*Psi >= 1)),
    color = "Initial State"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_visible_damage_brms.png"), p_frag_visible,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_visible_damage_brms.png\n")

# P(any damage | PGV)
p_frag_any <- fragility_single |>
  filter(initial_psi == 0) |>
  ggplot(aes(x = pgv, y = p_damage)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = p_damage - 0.1, ymax = pmin(p_damage + 0.1, 1)),
              alpha = 0.2, fill = "steelblue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Probability of Any Damage | PGV - brms Model",
    subtitle = "Single event, undamaged building",
    x = "PGV (mm/s)",
    y = expression(P(Delta*Psi > 0))
  ) +
  theme_minimal()

# Expected damage by PGV
p_expected <- fragility_single |>
  filter(initial_psi == 0) |>
  ggplot(aes(x = pgv)) +
  geom_ribbon(aes(ymin = q10_delta_psi, ymax = q90_delta_psi),
              alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = mean_delta_psi), linewidth = 1.2, color = "steelblue") +
  geom_line(aes(y = q50_delta_psi), linewidth = 1, color = "steelblue", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  annotate("text", x = 5, y = 1.1, label = "Visible damage", color = "red", hjust = 0) +
  scale_x_log10() +
  labs(
    title = "Expected Damage vs PGV - brms Model",
    subtitle = "Blue ribbon: 10-90% posterior interval",
    x = "PGV (mm/s)",
    y = expression(E[Delta*Psi])
  ) +
  theme_minimal()

combined_single <- p_frag_any + p_expected
ggsave(file.path(fig_dir, "fragility_single_event_brms.png"), combined_single,
       width = 14, height = 5, dpi = 150)
cat("Saved: fragility_single_event_brms.png\n")

# -----------------------------------------------------------------------------
# 3. Extract 10% exceedance thresholds
# -----------------------------------------------------------------------------

cat("\n=== 10% Exceedance Thresholds (brms) ===\n")

# Interpolate to find PGV where P(Psi >= 1) = 0.10
find_threshold <- function(data, target = 0.10) {
  if (max(data$p_visible) < target) return(NA)
  if (min(data$p_visible) > target) return(min(data$pgv))

  # Linear interpolation
  approx(data$p_visible, data$pgv, xout = target)$y
}

thresholds_brms <- fragility_single |>
  group_by(initial_psi) |>
  summarise(
    pgv_10pct = find_threshold(cur_data()),
    .groups = "drop"
  )

cat("\n10% exceedance thresholds (brms model):\n")
print(thresholds_brms)

cat("\nComparison with Korswagen:\n")
cat("  Undamaged (Psi0=0): brms =", round(thresholds_brms$pgv_10pct[1], 1),
    "mm/s (Korswagen ~13 mm/s)\n")
cat("  Pre-damaged (Psi0=1): brms =", round(thresholds_brms$pgv_10pct[3], 1),
    "mm/s (Korswagen ~6 mm/s)\n")

# -----------------------------------------------------------------------------
# 4. Fit lognormal fragility function
# -----------------------------------------------------------------------------

cat("\n=== Fitting Lognormal Fragility Function (brms) ===\n")

# For undamaged buildings
frag_undamaged <- fragility_single |>
  filter(initial_psi == 0, p_visible > 0, p_visible < 1)

fit_lognormal <- function(pgv, prob) {
  # Fit: P(Psi >= 1) = Phi((log(PGV) - log(median)) / beta)
  # Using probit link
  if (length(pgv) < 3) return(list(median = NA, beta = NA))

  z <- qnorm(pmax(0.001, pmin(0.999, prob)))
  log_pgv <- log(pgv)

  fit <- lm(z ~ log_pgv)

  beta <- 1 / coef(fit)[2]
  log_median <- -coef(fit)[1] * beta

  list(
    median = exp(log_median),
    beta = abs(beta)
  )
}

lognormal_params <- fit_lognormal(frag_undamaged$pgv, frag_undamaged$p_visible)

cat("\nLognormal fragility parameters (undamaged, brms):\n")
cat("  Median:", round(lognormal_params$median, 1), "mm/s\n")
cat("  Beta:", round(lognormal_params$beta, 2), "\n")

# Plot fitted vs empirical
pgv_fine <- seq(2, 128, length.out = 100)
p_fitted <- pnorm((log(pgv_fine) - log(lognormal_params$median)) / lognormal_params$beta)

p_lognormal <- tibble(pgv = pgv_fine, p_fitted = p_fitted) |>
  ggplot(aes(x = pgv, y = p_fitted)) +
  geom_line(color = "red", linewidth = 1.2) +
  geom_point(data = filter(fragility_single, initial_psi == 0),
             aes(x = pgv, y = p_visible), size = 3) +
  scale_x_log10() +
  labs(
    title = "Lognormal Fragility Fit - brms Model",
    subtitle = sprintf("Median = %.1f mm/s, β = %.2f",
                       lognormal_params$median, lognormal_params$beta),
    x = "PGV (mm/s)",
    y = expression(P(Delta*Psi >= 1))
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "fragility_lognormal_fit_brms.png"), p_lognormal,
       width = 8, height = 5, dpi = 150)
cat("Saved: fragility_lognormal_fit_brms.png\n")

# -----------------------------------------------------------------------------
# 5. Sensitivity to building parameters (brms posterior)
# -----------------------------------------------------------------------------

cat("\n=== Parameter Sensitivity Analysis ===\n")

# Vary facade type and soil profile
sensitivity_scenarios <- expand_grid(
  pgv = pgv_seq,
  initial_psi = 0,
  material = 1.0,
  facade_type = c("A", "B"),
  soil_profile = c("A", "B")
)

cat("Computing sensitivity analysis...\n")
sens_results <- list()

pb <- txtProgressBar(min = 0, max = nrow(sensitivity_scenarios), style = 3)
for (i in 1:nrow(sensitivity_scenarios)) {
  sc <- sensitivity_scenarios[i, ]

  pred <- predict_damage_brms(
    fit_brms = fit_brms,
    pgv = sc$pgv,
    n = 1,
    initial_psi = sc$initial_psi,
    material = sc$material,
    facade_type = sc$facade_type,
    soil_profile = sc$soil_profile,
    ndraws = 500,
    summary = FALSE
  )

  sens_results[[i]] <- tibble(
    pgv = sc$pgv,
    facade_type = sc$facade_type,
    soil_profile = sc$soil_profile,
    p_visible = mean(pred >= 1),
    mean_damage = mean(pred)
  )

  setTxtProgressBar(pb, i)
}
close(pb)

sensitivity_df <- bind_rows(sens_results)

# Plot sensitivity
p_sens <- sensitivity_df |>
  mutate(config = paste0("Facade ", facade_type, ", Soil ", soil_profile)) |>
  ggplot(aes(x = pgv, y = p_visible, color = config)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_log10() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Fragility Sensitivity to Building Parameters - brms",
    x = "PGV (mm/s)",
    y = expression(P(Delta*Psi >= 1)),
    color = "Configuration"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_sensitivity_brms.png"), p_sens,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_sensitivity_brms.png\n")

# -----------------------------------------------------------------------------
# 6. Multi-event cumulative fragility
# -----------------------------------------------------------------------------

cat("\n=== Multi-Event Cumulative Fragility ===\n")

# Simulate cumulative damage over multiple events at same PGV
n_events_seq <- c(1, 2, 4, 8)
pgv_test <- c(8, 16, 32)

cumulative_results <- list()

for (pgv_val in pgv_test) {
  for (n_ev in n_events_seq) {
    # Simulate n_ev events at this PGV
    n_sims <- 500
    final_psi <- numeric(n_sims)

    for (sim in 1:n_sims) {
      psi <- 0
      for (ev in 1:n_ev) {
        delta <- simulate_damage_increment_brms(
          fit_brms = fit_brms,
          pgv = pgv_val,
          current_psi = psi,
          ndraws = 1
        )
        psi <- psi + delta
      }
      final_psi[sim] <- psi
    }

    cumulative_results[[length(cumulative_results) + 1]] <- tibble(
      pgv = pgv_val,
      n_events = n_ev,
      p_visible = mean(final_psi >= 1),
      mean_psi = mean(final_psi),
      sd_psi = sd(final_psi)
    )
  }
}

cumulative_df <- bind_rows(cumulative_results)

cat("\nCumulative damage after multiple events:\n")
print(cumulative_df)

p_cumulative <- cumulative_df |>
  mutate(pgv_label = paste0("PGV = ", pgv, " mm/s")) |>
  ggplot(aes(x = n_events, y = p_visible, color = pgv_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_viridis_d() +
  labs(
    title = "Cumulative Fragility: P(Psi >= 1) vs Number of Events - brms",
    x = "Number of Events (same PGV)",
    y = expression(P(Psi >= 1)),
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_cumulative_brms.png"), p_cumulative,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_cumulative_brms.png\n")

# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------

saveRDS(fragility_single, here::here("outputs", "models", "fragility_curves_brms.rds"))

fragility_params_brms <- list(
  thresholds = thresholds_brms,
  lognormal_undamaged = lognormal_params,
  cumulative = cumulative_df
)
saveRDS(fragility_params_brms, here::here("outputs", "models", "fragility_parameters_brms.rds"))

# Summary table
cat("\n=== Fragility Curve Summary (brms) ===\n")
cat("\n10% Exceedance Thresholds:\n")
print(thresholds_brms)
cat("\nLognormal Parameters (undamaged):\n")
cat("  Median:", round(lognormal_params$median, 1), "mm/s\n")
cat("  Beta:", round(lognormal_params$beta, 2), "\n")

cat("\nFragility curves with brms complete.\n")
