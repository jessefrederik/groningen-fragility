# =============================================================================
# 03_fit_baseline.R
#
# Fit baseline GAM models using mgcv for fast iteration.
# Uses hurdle structure: P(damage) + E[damage|damage>0]
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(mgcv)

# -----------------------------------------------------------------------------
# 1. Load prepared data
# -----------------------------------------------------------------------------

fem <- readRDS(here::here("outputs", "models", "fem_prepared.rds"))

cat("Data loaded:", nrow(fem), "observations\n")
cat("Positive damage:", sum(fem$has_damage), "observations\n")

# -----------------------------------------------------------------------------
# 2. Fit Stage 1: Probability of damage (hurdle)
# -----------------------------------------------------------------------------

cat("\n=== Fitting Stage 1: P(damage > 0) ===\n")

# Note on monotonicity:
# mgcv's "mpi" basis is monotonic increasing P-spline
# However, it's not in base mgcv - we use "tp" with penalty and check shape

# Standard GAM for initiation probability
m1_init <- gam(
  has_damage ~ s(log_pgv, k = 6, bs = "tp") +
               s(log_n, k = 4, bs = "tp") +
               s(InitialPsi, k = 5, bs = "tp") +
               Material +
               FacadeType +
               SoilProfile +
               EarthquakeType,
  family = binomial(link = "logit"),
  data = fem,
  method = "REML"
)

cat("\nStage 1 Summary:\n")
print(summary(m1_init))

# Check effective degrees of freedom
cat("\nEffective degrees of freedom:\n")
cat("log_pgv:", round(sum(m1_init$edf[1]), 2), "\n")
cat("log_n:", round(sum(m1_init$edf[2]), 2), "\n")
cat("InitialPsi:", round(sum(m1_init$edf[3]), 2), "\n")

# -----------------------------------------------------------------------------
# 3. Fit Stage 2: Damage amount | damage > 0
# -----------------------------------------------------------------------------

cat("\n=== Fitting Stage 2: E[damage | damage > 0] ===\n")

# Subset to positive damage only
fem_pos <- fem |> filter(has_damage == 1)
cat("Observations with damage:", nrow(fem_pos), "\n")

# Gamma GAM with log link
m2_severity <- gam(
  delta_psi_pos ~ s(log_pgv, k = 6, bs = "tp") +
                  s(log_n, k = 4, bs = "tp") +
                  s(InitialPsi, k = 5, bs = "tp") +
                  Material +
                  FacadeType +
                  SoilProfile +
                  EarthquakeType,
  family = Gamma(link = "log"),
  data = fem_pos,
  method = "REML"
)

cat("\nStage 2 Summary:\n")
print(summary(m2_severity))

# Gamma shape parameter
gamma_scale <- m2_severity$scale
gamma_shape <- 1 / gamma_scale
cat("\nGamma shape parameter:", round(gamma_shape, 3), "\n")
cat("Gamma scale (dispersion):", round(gamma_scale, 3), "\n")

# -----------------------------------------------------------------------------
# 4. Visualize smooth effects
# -----------------------------------------------------------------------------

cat("\n=== Generating Effect Plots ===\n")

fig_dir <- here::here("outputs", "figures")

# Stage 1: Initiation probability
png(file.path(fig_dir, "baseline_stage1_effects.png"), width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2))
plot(m1_init, select = 1, shade = TRUE, main = "Effect of log(PGV) on P(damage)")
plot(m1_init, select = 2, shade = TRUE, main = "Effect of log(N) on P(damage)")
plot(m1_init, select = 3, shade = TRUE, main = "Effect of InitialPsi on P(damage)")
dev.off()
cat("Saved: baseline_stage1_effects.png\n")

# Stage 2: Severity
png(file.path(fig_dir, "baseline_stage2_effects.png"), width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2))
plot(m2_severity, select = 1, shade = TRUE, main = "Effect of log(PGV) on E[damage]")
plot(m2_severity, select = 2, shade = TRUE, main = "Effect of log(N) on E[damage]")
plot(m2_severity, select = 3, shade = TRUE, main = "Effect of InitialPsi on E[damage]")
dev.off()
cat("Saved: baseline_stage2_effects.png\n")

# -----------------------------------------------------------------------------
# 5. Check monotonicity
# -----------------------------------------------------------------------------

cat("\n=== Checking Monotonicity ===\n")

# Create prediction grid for log_pgv
pgv_grid <- tibble(
  log_pgv = seq(min(fem$log_pgv), max(fem$log_pgv), length.out = 100),
  log_n = mean(fem$log_n),
  InitialPsi = mean(fem$InitialPsi),
  Material = 1.0,
  FacadeType = "A",
  SoilProfile = "A",
  EarthquakeType = "ZN"
)

# Predictions
pgv_grid <- pgv_grid |>
  mutate(
    pred_init = predict(m1_init, newdata = pgv_grid, type = "response"),
    pred_severity = predict(m2_severity, newdata = pgv_grid, type = "response"),
    pred_combined = pred_init * pred_severity  # Expected damage (hurdle)
  )

# Check for monotonicity violations
init_monotone <- all(diff(pgv_grid$pred_init) >= -1e-6)
sev_monotone <- all(diff(pgv_grid$pred_severity) >= -1e-6)

cat("P(damage) monotonic in PGV:", init_monotone, "\n")
cat("E[damage|>0] monotonic in PGV:", sev_monotone, "\n")

if (!init_monotone || !sev_monotone) {
  cat("WARNING: Consider using monotonic splines (scam package) or brms mo() terms\n")
}

# Plot predicted curves
p_mono <- pgv_grid |>
  pivot_longer(cols = starts_with("pred_"), names_to = "component", values_to = "value") |>
  mutate(component = case_when(
    component == "pred_init" ~ "P(damage > 0)",
    component == "pred_severity" ~ "E[damage | damage > 0]",
    component == "pred_combined" ~ "E[damage] (hurdle)"
  )) |>
  ggplot(aes(x = exp(log_pgv), y = value, color = component)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  facet_wrap(~component, scales = "free_y", ncol = 1) +
  labs(
    title = "Baseline GAM: Predicted Damage vs PGV",
    subtitle = "Holding other predictors at mean/reference levels",
    x = "PGV (mm/s, log scale)",
    y = "Predicted Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "baseline_monotonicity_check.png"), p_mono,
       width = 8, height = 10, dpi = 150)
cat("Saved: baseline_monotonicity_check.png\n")

# -----------------------------------------------------------------------------
# 6. Model diagnostics
# -----------------------------------------------------------------------------

cat("\n=== Model Diagnostics ===\n")

# Stage 1 diagnostics
png(file.path(fig_dir, "baseline_stage1_diagnostics.png"), width = 1000, height = 800, res = 150)
par(mfrow = c(2, 2))
gam.check(m1_init)
dev.off()
cat("Saved: baseline_stage1_diagnostics.png\n")

# Stage 2 diagnostics
png(file.path(fig_dir, "baseline_stage2_diagnostics.png"), width = 1000, height = 800, res = 150)
par(mfrow = c(2, 2))
gam.check(m2_severity)
dev.off()
cat("Saved: baseline_stage2_diagnostics.png\n")

# -----------------------------------------------------------------------------
# 7. Coefficient summaries
# -----------------------------------------------------------------------------

cat("\n=== Parametric Coefficients ===\n")

cat("\nStage 1 (logit scale):\n")
coef_init <- coef(m1_init)
coef_init <- coef_init[!grepl("^s\\(", names(coef_init))]
print(round(coef_init, 3))

cat("\nStage 2 (log scale):\n")
coef_sev <- coef(m2_severity)
coef_sev <- coef_sev[!grepl("^s\\(", names(coef_sev))]
print(round(coef_sev, 3))

# Interpret Material effect
cat("\n=== Material Effect Interpretation ===\n")
cat("Stage 1: Material coefficient =", round(coef_init["Material"], 3), "\n")
cat("  (negative = stronger material reduces P(damage))\n")
cat("Stage 2: Material coefficient =", round(coef_sev["Material"], 3), "\n")
cat("  (negative = stronger material reduces damage severity)\n")

# -----------------------------------------------------------------------------
# 8. Save models
# -----------------------------------------------------------------------------

saveRDS(m1_init, here::here("outputs", "models", "baseline_stage1_init.rds"))
saveRDS(m2_severity, here::here("outputs", "models", "baseline_stage2_severity.rds"))

cat("\nModels saved to outputs/models/\n")

# -----------------------------------------------------------------------------
# 9. Create prediction function
# -----------------------------------------------------------------------------

#' Predict damage from baseline hurdle model
#'
#' @param newdata Data frame with predictors
#' @param m1 Stage 1 model (initiation)
#' @param m2 Stage 2 model (severity)
#' @param n_sims Number of simulations for uncertainty
#' @return Data frame with predictions
predict_hurdle_baseline <- function(newdata, m1 = m1_init, m2 = m2_severity, n_sims = 0) {
  # Point predictions
  p_init <- predict(m1, newdata = newdata, type = "response")
  mu_sev <- predict(m2, newdata = newdata, type = "response")

  result <- newdata |>
    mutate(
      p_damage = p_init,
      mu_if_damage = mu_sev,
      expected_damage = p_init * mu_sev
    )

  if (n_sims > 0) {
    # Simulate from hurdle distribution
    gamma_shape <- 1 / m2$scale
    sims <- matrix(NA, nrow = nrow(newdata), ncol = n_sims)

    for (i in 1:n_sims) {
      # Stage 1: Bernoulli
      damage_occurs <- rbinom(nrow(newdata), 1, p_init)
      # Stage 2: Gamma
      damage_amount <- rgamma(nrow(newdata), shape = gamma_shape, rate = gamma_shape / mu_sev)
      sims[, i] <- damage_occurs * damage_amount
    }

    result <- result |>
      mutate(
        median_damage = apply(sims, 1, median),
        q10_damage = apply(sims, 1, quantile, 0.1),
        q90_damage = apply(sims, 1, quantile, 0.9)
      )
  }

  return(result)
}

# Save prediction function
save(predict_hurdle_baseline, file = here::here("outputs", "models", "predict_hurdle_baseline.RData"))
cat("Prediction function saved.\n")

cat("\nBaseline model fitting complete.\n")
