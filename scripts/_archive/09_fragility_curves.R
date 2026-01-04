# =============================================================================
# 09_fragility_curves.R
#
# Derive fragility curves P(Psi >= threshold | PGV) from fitted models.
# Compares to Korswagen's published curves.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(mgcv)

# Load models and parameters
m1_init <- readRDS(here::here("outputs", "models", "baseline_stage1_init.rds"))
m2_severity <- readRDS(here::here("outputs", "models", "baseline_stage2_severity.rds"))
m1_S <- readRDS(here::here("outputs", "models", "accumulation_stage1_S.rds"))
m2_S <- readRDS(here::here("outputs", "models", "accumulation_stage2_S.rds"))
accum_params <- readRDS(here::here("outputs", "models", "accumulation_params.rds"))

fig_dir <- here::here("outputs", "figures")

# -----------------------------------------------------------------------------
# 1. Define fragility curve computation
# -----------------------------------------------------------------------------

#' Compute fragility curve via Monte Carlo simulation
#'
#' @param pgv_values Vector of PGV values to evaluate
#' @param threshold Damage threshold (e.g., 1 for visible damage)
#' @param initial_psi Initial damage state
#' @param material Material strength
#' @param facade_type "A" or "B"
#' @param soil_profile "A" or "B"
#' @param n_sims Number of Monte Carlo samples
#' @param m1 Stage 1 model
#' @param m2 Stage 2 model
#' @param n_events Number of equivalent events (for accumulation)
#' @return Tibble with PGV, P(exceed), and confidence bounds
compute_fragility_curve <- function(
  pgv_values = seq(2, 100, by = 2),
  threshold = 1,
  initial_psi = 0,
  material = 1.0,
  facade_type = "A",
  soil_profile = "A",
  n_sims = 5000,
  m1 = m1_init,
  m2 = m2_severity,
  n_events = 1
) {
  gamma_shape <- 1 / m2$scale
  n_pgv <- length(pgv_values)

  # Results storage
  p_exceed <- numeric(n_pgv)
  p_exceed_lower <- numeric(n_pgv)
  p_exceed_upper <- numeric(n_pgv)

  for (i in seq_along(pgv_values)) {
    pgv <- pgv_values[i]

    # Prepare prediction data
    pred_data <- tibble(
      log_pgv = log(pgv),
      log_n = log(n_events),
      InitialPsi = initial_psi,
      Material = material,
      FacadeType = facade_type,
      SoilProfile = soil_profile,
      EarthquakeType = "ZN"  # Use one record; could integrate over
    )

    # Get predictions
    p_damage <- predict(m1, newdata = pred_data, type = "response")[1]
    mu_damage <- predict(m2, newdata = pred_data, type = "response")[1]

    # Monte Carlo simulation
    damage_occurs <- rbinom(n_sims, 1, p_damage)
    damage_amount <- rgamma(n_sims, shape = gamma_shape, rate = gamma_shape / mu_damage)
    final_psi <- initial_psi + damage_occurs * damage_amount

    # Proportion exceeding threshold
    p_exceed[i] <- mean(final_psi >= threshold)

    # Bootstrap confidence interval
    boot_props <- replicate(500, {
      idx <- sample(n_sims, n_sims, replace = TRUE)
      mean(final_psi[idx] >= threshold)
    })
    p_exceed_lower[i] <- quantile(boot_props, 0.025)
    p_exceed_upper[i] <- quantile(boot_props, 0.975)
  }

  tibble(
    pgv = pgv_values,
    p_exceed = p_exceed,
    p_lower = p_exceed_lower,
    p_upper = p_exceed_upper,
    threshold = threshold,
    initial_psi = initial_psi,
    material = material
  )
}

# -----------------------------------------------------------------------------
# 2. Generate fragility curves for key scenarios
# -----------------------------------------------------------------------------

cat("\n=== Generating Fragility Curves ===\n")

pgv_grid <- seq(2, 80, by = 2)

# Scenario 1: Undamaged, standard material, visible damage (Psi >= 1)
cat("Computing: Undamaged, standard material, Psi >= 1...\n")
frag_undamaged_std <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 1,
  initial_psi = 0,
  material = 1.0
)

# Scenario 2: Pre-damaged (Psi_0 = 0.5), standard material
cat("Computing: Pre-damaged (0.5), standard material, Psi >= 1...\n")
frag_predamaged_std <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 1,
  initial_psi = 0.5,
  material = 1.0
)

# Scenario 3: Undamaged, weak material
cat("Computing: Undamaged, weak material, Psi >= 1...\n")
frag_undamaged_weak <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 1,
  initial_psi = 0,
  material = 0.7
)

# Scenario 4: Undamaged, strong material
cat("Computing: Undamaged, strong material, Psi >= 1...\n")
frag_undamaged_strong <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 1,
  initial_psi = 0,
  material = 1.3
)

# Scenario 5: Multiple events (N=4)
cat("Computing: Undamaged, N=4 events, Psi >= 1...\n")
frag_n4_std <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 1,
  initial_psi = 0,
  material = 1.0,
  n_events = 4
)

# Scenario 6: Moderate damage threshold (Psi >= 2)
cat("Computing: Undamaged, Psi >= 2...\n")
frag_ds2 <- compute_fragility_curve(
  pgv_values = pgv_grid,
  threshold = 2,
  initial_psi = 0,
  material = 1.0
)

# Combine all curves
all_curves <- bind_rows(
  frag_undamaged_std |> mutate(scenario = "Undamaged, m=1.0"),
  frag_predamaged_std |> mutate(scenario = "Pre-damaged (Psi0=0.5), m=1.0"),
  frag_undamaged_weak |> mutate(scenario = "Undamaged, m=0.7 (weak)"),
  frag_undamaged_strong |> mutate(scenario = "Undamaged, m=1.3 (strong)"),
  frag_n4_std |> mutate(scenario = "Undamaged, N=4 events"),
  frag_ds2 |> mutate(scenario = "Undamaged, Psi>=2 (DS2)")
)

# -----------------------------------------------------------------------------
# 3. Plot fragility curves
# -----------------------------------------------------------------------------

cat("\n=== Generating Plots ===\n")

# Main fragility curve comparison
p_frag_main <- all_curves |>
  filter(threshold == 1) |>
  ggplot(aes(x = pgv, y = p_exceed, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Fragility Curves: P(Visible Damage | PGV)",
    subtitle = expression("Threshold: " * Psi >= 1 * " (DS1, cracks > 0.1mm)"),
    x = "PGV (mm/s, log scale)",
    y = expression(P(Psi >= 1)),
    color = "Scenario",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_curves_main.png"), p_frag_main,
       width = 10, height = 7, dpi = 150)
cat("Saved: fragility_curves_main.png\n")

# DS1 vs DS2 comparison
p_frag_ds <- all_curves |>
  filter(scenario %in% c("Undamaged, m=1.0", "Undamaged, Psi>=2 (DS2)")) |>
  mutate(ds = if_else(threshold == 1, "DS1 (Psi>=1)", "DS2 (Psi>=2)")) |>
  ggplot(aes(x = pgv, y = p_exceed, color = ds, fill = ds)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() +
  scale_color_manual(values = c("DS1 (Psi>=1)" = "steelblue", "DS2 (Psi>=2)" = "coral")) +
  scale_fill_manual(values = c("DS1 (Psi>=1)" = "steelblue", "DS2 (Psi>=2)" = "coral")) +
  labs(
    title = "Fragility Curves: DS1 vs DS2",
    subtitle = "Standard material, undamaged initial condition",
    x = "PGV (mm/s, log scale)",
    y = "Exceedance Probability",
    color = "Damage State",
    fill = "Damage State"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "fragility_curves_ds_comparison.png"), p_frag_ds,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_curves_ds_comparison.png\n")

# Effect of initial damage
p_frag_init <- all_curves |>
  filter(scenario %in% c("Undamaged, m=1.0", "Pre-damaged (Psi0=0.5), m=1.0")) |>
  ggplot(aes(x = pgv, y = p_exceed, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() +
  labs(
    title = "Effect of Initial Damage on Fragility",
    subtitle = "Pre-existing damage increases vulnerability",
    x = "PGV (mm/s, log scale)",
    y = expression(P(Psi >= 1)),
    color = "Initial Condition",
    fill = "Initial Condition"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "fragility_curves_initial_damage.png"), p_frag_init,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_curves_initial_damage.png\n")

# -----------------------------------------------------------------------------
# 4. Extract key PGV values at specified exceedance probabilities
# -----------------------------------------------------------------------------

cat("\n=== Key PGV Values ===\n")

# Function to interpolate PGV at given probability
find_pgv_at_prob <- function(curve, target_prob) {
  if (max(curve$p_exceed) < target_prob) return(NA)
  if (min(curve$p_exceed) > target_prob) return(min(curve$pgv))

  idx <- which.min(abs(curve$p_exceed - target_prob))

  # Linear interpolation around closest point
  if (idx == 1) return(curve$pgv[1])
  if (idx == nrow(curve)) return(curve$pgv[nrow(curve)])

  if (curve$p_exceed[idx] > target_prob) {
    x1 <- curve$pgv[idx-1]; x2 <- curve$pgv[idx]
    y1 <- curve$p_exceed[idx-1]; y2 <- curve$p_exceed[idx]
  } else {
    x1 <- curve$pgv[idx]; x2 <- curve$pgv[idx+1]
    y1 <- curve$p_exceed[idx]; y2 <- curve$p_exceed[idx+1]
  }

  # Interpolate
  x1 + (target_prob - y1) * (x2 - x1) / (y2 - y1)
}

# Extract PGV at 10% and 50% exceedance
key_values <- all_curves |>
  filter(threshold == 1) |>
  group_by(scenario) |>
  summarise(
    pgv_10pct = find_pgv_at_prob(cur_data(), 0.10),
    pgv_50pct = find_pgv_at_prob(cur_data(), 0.50),
    .groups = "drop"
  )

cat("\nPGV at 10% and 50% exceedance probability (Psi >= 1):\n")
print(key_values)

# Compare to Korswagen's values
cat("\n=== Comparison to Korswagen Paper ===\n")
cat("Korswagen findings (from abstract):\n")
cat("- 10% exceedance of visible damage at 13 mm/s (undamaged)\n")
cat("- 10% exceedance at 6 mm/s (if pre-damaged)\n")
cat("- 1% exceedance of moderate damage at >15 mm/s\n")

cat("\nOur model estimates:\n")
cat("- 10% exceedance (undamaged, m=1.0):",
    round(key_values$pgv_10pct[key_values$scenario == "Undamaged, m=1.0"], 1), "mm/s\n")
cat("- 10% exceedance (pre-damaged):",
    round(key_values$pgv_10pct[key_values$scenario == "Pre-damaged (Psi0=0.5), m=1.0"], 1), "mm/s\n")

# -----------------------------------------------------------------------------
# 5. Fragility function parameterization (for external use)
# -----------------------------------------------------------------------------

cat("\n=== Fitting Parametric Fragility Functions ===\n")

# Fit lognormal fragility curve: P = Phi((log(PGV) - log(median)) / beta)
fit_lognormal_fragility <- function(curve) {
  # Use optimization to fit lognormal CDF
  nll <- function(params) {
    median_pgv <- params[1]
    beta <- params[2]

    p_pred <- pnorm((log(curve$pgv) - log(median_pgv)) / beta)
    -sum(dbinom(round(curve$p_exceed * 5000), 5000, p_pred, log = TRUE))
  }

  opt <- optim(c(20, 0.5), nll, method = "L-BFGS-B",
               lower = c(1, 0.1), upper = c(100, 2))

  list(
    median = opt$par[1],
    beta = opt$par[2]
  )
}

# Fit to main scenarios
params_undamaged <- fit_lognormal_fragility(frag_undamaged_std)
params_predamaged <- fit_lognormal_fragility(frag_predamaged_std)

cat("\nLognormal fragility parameters (Psi >= 1):\n")
cat("Undamaged: median =", round(params_undamaged$median, 1), "mm/s, beta =",
    round(params_undamaged$beta, 3), "\n")
cat("Pre-damaged: median =", round(params_predamaged$median, 1), "mm/s, beta =",
    round(params_predamaged$beta, 3), "\n")

# Plot fitted vs empirical
p_fit <- frag_undamaged_std |>
  mutate(
    p_fitted = pnorm((log(pgv) - log(params_undamaged$median)) / params_undamaged$beta)
  ) |>
  ggplot(aes(x = pgv)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_upper), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = p_exceed, color = "Simulated"), linewidth = 1) +
  geom_line(aes(y = p_fitted, color = "Lognormal fit"), linetype = "dashed", linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "Fragility Curve: Simulated vs Lognormal Fit",
    subtitle = sprintf("Median = %.1f mm/s, beta = %.3f",
                       params_undamaged$median, params_undamaged$beta),
    x = "PGV (mm/s, log scale)",
    y = expression(P(Psi >= 1)),
    color = ""
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "fragility_lognormal_fit.png"), p_fit,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_lognormal_fit.png\n")

# -----------------------------------------------------------------------------
# 6. Save fragility curves and parameters
# -----------------------------------------------------------------------------

# Save all curves
saveRDS(all_curves, here::here("outputs", "models", "fragility_curves.rds"))

# Save parametric fits
fragility_params <- list(
  undamaged_ds1 = params_undamaged,
  predamaged_ds1 = params_predamaged,
  key_values = key_values
)
saveRDS(fragility_params, here::here("outputs", "models", "fragility_parameters.rds"))

# Create summary table for export
fragility_summary <- all_curves |>
  filter(threshold == 1, pgv %in% c(5, 10, 15, 20, 30, 50)) |>
  select(scenario, pgv, p_exceed, p_lower, p_upper) |>
  pivot_wider(
    names_from = pgv,
    values_from = c(p_exceed, p_lower, p_upper),
    names_glue = "PGV{pgv}_{.value}"
  )

write_csv(fragility_summary, here::here("outputs", "fragility_summary.csv"))
cat("Saved: fragility_summary.csv\n")

cat("\n=== Fragility Analysis Complete ===\n")
cat("Key outputs:\n")
cat("- Fragility curves: outputs/models/fragility_curves.rds\n")
cat("- Lognormal parameters: outputs/models/fragility_parameters.rds\n")
cat("- Summary table: outputs/fragility_summary.csv\n")
cat("- Plots: outputs/figures/fragility_*.png\n")
