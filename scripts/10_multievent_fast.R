# =============================================================================
# 10_multievent_fast.R
#
# Fast multi-event cumulative fragility using pre-extracted posterior samples.
# Avoids repeated brms::posterior_predict() calls.
#
# @author Generated with Claude Code
# @date 2026-01-04
# =============================================================================

library(tidyverse)
library(brms)

# Load brms model
fit_brms <- readRDS(here::here("outputs", "models", "brms_hurdle_gamma_A.rds"))
std_params <- readRDS(here::here("outputs", "models", "standardization_params.rds"))

fig_dir <- here::here("outputs", "figures")

cat("=== Fast Multi-Event Cumulative Fragility ===\n")

# -----------------------------------------------------------------------------
# 1. Pre-extract posterior samples (do this ONCE)
# -----------------------------------------------------------------------------

cat("Extracting posterior samples...\n")

# Get posterior draws as data frame
post <- as_draws_df(fit_brms)
n_post <- nrow(post)

cat("Extracted", n_post, "posterior draws\n")

# Key parameters for hurdle-gamma
# Severity (mu): b_Intercept + bsp_mopgv_ord * pgv_effect + ...
# Hurdle (hu): b_hu_Intercept + bsp_hu_mopgv_ord * pgv_effect + ...

# Extract fixed effects
b_intercept <- post$b_Intercept
b_initial_psi <- post$b_initial_psi_c
bsp_pgv <- post$bsp_mopgv_ord
bsp_n <- post$bsp_mon_ord
shape <- post$shape

hu_intercept <- post$b_hu_Intercept
hu_initial_psi <- post$b_hu_initial_psi_c
hu_bsp_pgv <- post$bsp_hu_mopgv_ord
hu_bsp_n <- post$bsp_hu_mon_ord

# Simplex weights for monotonic effects (cumulative)
# PGV has 7 increments (8 levels - 1)
pgv_simplex <- matrix(NA, n_post, 7)
for (i in 1:7) {
  pgv_simplex[, i] <- post[[paste0("simo_mopgv_ord1[", i, "]")]]
}
pgv_cumsum <- t(apply(pgv_simplex, 1, cumsum))  # Cumulative for each level

hu_pgv_simplex <- matrix(NA, n_post, 7)
for (i in 1:7) {
  hu_pgv_simplex[, i] <- post[[paste0("simo_hu_mopgv_ord1[", i, "]")]]
}
hu_pgv_cumsum <- t(apply(hu_pgv_simplex, 1, cumsum))

# N simplex (4 increments)
n_simplex <- matrix(NA, n_post, 4)
for (i in 1:4) {
  n_simplex[, i] <- post[[paste0("simo_mon_ord1[", i, "]")]]
}
n_cumsum <- t(apply(n_simplex, 1, cumsum))

cat("Posterior parameters extracted.\n")

# -----------------------------------------------------------------------------
# 2. Fast prediction function using pre-extracted posteriors
# -----------------------------------------------------------------------------

# Map PGV to level index (1-8)
pgv_to_idx <- function(pgv) {
  pgv_levels <- c(2, 4, 8, 16, 32, 64, 96, 128)
  sapply(pgv, function(p) which.min(abs(pgv_levels - p)))
}

#' Fast damage prediction using pre-extracted posteriors
#' @param pgv PGV value (scalar)
#' @param initial_psi Initial damage (scalar)
#' @param n_draws Number of posterior draws to use
#' @return Vector of damage samples
fast_predict_damage <- function(pgv, initial_psi, n_draws = 500) {
  # Sample posterior indices
  idx <- sample(n_post, n_draws, replace = TRUE)

  # Get PGV level index (1-8)
  pgv_idx <- pgv_to_idx(pgv)

  # Calculate monotonic PGV effect (0 for level 1, cumsum for others)
  if (pgv_idx == 1) {
    pgv_effect <- rep(0, n_draws)
    hu_pgv_effect <- rep(0, n_draws)
  } else {
    pgv_effect <- pgv_cumsum[idx, pgv_idx - 1]
    hu_pgv_effect <- hu_pgv_cumsum[idx, pgv_idx - 1]
  }

  # N effect (single event = level 1 = 0)
  n_effect <- rep(0, n_draws)

  # Center initial_psi
  psi_c <- initial_psi - std_params$initial_psi_mean

  # Linear predictor for severity (log scale)
  eta_mu <- b_intercept[idx] +
            bsp_pgv[idx] * pgv_effect +
            bsp_n[idx] * n_effect +
            b_initial_psi[idx] * psi_c

  # Linear predictor for hurdle (logit scale)
  eta_hu <- hu_intercept[idx] +
            hu_bsp_pgv[idx] * hu_pgv_effect +
            hu_bsp_n[idx] * n_effect +
            hu_initial_psi[idx] * psi_c

  # Transform
  mu <- exp(eta_mu)
  p_zero <- plogis(eta_hu)  # P(Y = 0)

  # Sample from hurdle-gamma
  is_zero <- rbinom(n_draws, 1, p_zero)
  gamma_sample <- rgamma(n_draws, shape = shape[idx], rate = shape[idx] / mu)

  # Hurdle: zero or gamma
  damage <- ifelse(is_zero == 1, 0, gamma_sample)

  return(damage)
}

# Test it
cat("\nTesting fast prediction...\n")
system.time({
  test <- fast_predict_damage(pgv = 32, initial_psi = 0, n_draws = 1000)
})
cat("Mean damage at PGV=32:", mean(test), "\n")
cat("P(damage >= 1):", mean(test >= 1), "\n")

# -----------------------------------------------------------------------------
# 3. Fast multi-event simulation
# -----------------------------------------------------------------------------

cat("\n=== Running Multi-Event Simulation ===\n")

# Parameters
pgv_values <- c(8, 16, 32, 64)
n_events_seq <- c(1, 2, 4, 8, 16)
n_sims <- 1000

results <- list()
counter <- 0
total <- length(pgv_values) * length(n_events_seq)

for (pgv_val in pgv_values) {
  for (n_ev in n_events_seq) {
    counter <- counter + 1
    cat(sprintf("\r[%d/%d] PGV=%d, N=%d", counter, total, pgv_val, n_ev))

    # Simulate n_sims trajectories
    final_psi <- numeric(n_sims)

    for (sim in 1:n_sims) {
      psi <- 0
      for (ev in 1:n_ev) {
        # Single posterior draw per event
        delta <- fast_predict_damage(pgv_val, psi, n_draws = 1)
        psi <- psi + delta
      }
      final_psi[sim] <- psi
    }

    results[[length(results) + 1]] <- tibble(
      pgv = pgv_val,
      n_events = n_ev,
      mean_psi = mean(final_psi),
      sd_psi = sd(final_psi),
      p_visible = mean(final_psi >= 1),
      p_moderate = mean(final_psi >= 2),
      q10 = quantile(final_psi, 0.10),
      q50 = quantile(final_psi, 0.50),
      q90 = quantile(final_psi, 0.90)
    )
  }
}
cat("\n")

cumulative_df <- bind_rows(results)

cat("\nMulti-event cumulative damage results:\n")
print(cumulative_df, n = 25)

# -----------------------------------------------------------------------------
# 4. Visualize results
# -----------------------------------------------------------------------------

cat("\n=== Generating Plots ===\n")

# P(visible damage) vs number of events
p_cumulative <- cumulative_df |>
  mutate(pgv_label = paste0("PGV = ", pgv, " mm/s")) |>
  ggplot(aes(x = n_events, y = p_visible, color = pgv_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  scale_color_viridis_d(option = "C") +
  scale_x_continuous(breaks = n_events_seq) +
  labs(
    title = "Cumulative Fragility: P(Visible Damage) vs Number of Events",
    subtitle = "1000 Monte Carlo simulations per scenario",
    x = "Number of Events (same PGV)",
    y = expression(P(Psi >= 1)),
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_cumulative_fast.png"), p_cumulative,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_cumulative_fast.png\n")

# Expected damage vs events
p_expected <- cumulative_df |>
  mutate(pgv_label = paste0("PGV = ", pgv, " mm/s")) |>
  ggplot(aes(x = n_events, y = mean_psi, color = pgv_label)) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = pgv_label), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 14, y = 1.1, label = "Visible damage", color = "red") +
  scale_color_viridis_d(option = "C") +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(breaks = n_events_seq) +
  labs(
    title = "Expected Cumulative Damage vs Number of Events",
    subtitle = "Ribbon: 10-90% interval",
    x = "Number of Events",
    y = expression(E[Psi]),
    color = NULL, fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fragility_cumulative_expected.png"), p_expected,
       width = 10, height = 6, dpi = 150)
cat("Saved: fragility_cumulative_expected.png\n")

# Heatmap of P(visible)
p_heatmap <- cumulative_df |>
  ggplot(aes(x = factor(n_events), y = factor(pgv), fill = p_visible)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p_visible)), color = "white", size = 4) +
  scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
  labs(
    title = "P(Visible Damage) by PGV and Number of Events",
    x = "Number of Events",
    y = "PGV (mm/s)",
    fill = expression(P(Psi >= 1))
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "fragility_cumulative_heatmap.png"), p_heatmap,
       width = 8, height = 6, dpi = 150)
cat("Saved: fragility_cumulative_heatmap.png\n")

# -----------------------------------------------------------------------------
# 5. Save results
# -----------------------------------------------------------------------------

saveRDS(cumulative_df, here::here("outputs", "models", "cumulative_fragility_fast.rds"))

cat("\n=== Summary ===\n")
cat("Events needed for 10% exceedance (P >= 0.1):\n")
cumulative_df |>
  filter(p_visible >= 0.1) |>
  group_by(pgv) |>
  summarise(min_events = min(n_events), .groups = "drop") |>
  print()

cat("\nMulti-event analysis complete.\n")
