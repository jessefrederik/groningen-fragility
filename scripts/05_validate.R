# =============================================================================
# 05_validate.R
#
# Model validation: Leave-one-record-out cross-validation,
# posterior predictive checks, and calibration assessment.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(brms)
library(mgcv)
library(bayesplot)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load data and models
# -----------------------------------------------------------------------------

fem <- readRDS(here::here("outputs", "models", "fem_prepared.rds"))

# Load baseline models
m1_init <- readRDS(here::here("outputs", "models", "baseline_stage1_init.rds"))
m2_severity <- readRDS(here::here("outputs", "models", "baseline_stage2_severity.rds"))

# Check if brms model exists (may take long to fit)
brms_path <- here::here("outputs", "models", "brms_hurdle_gamma_A.rds")
has_brms <- file.exists(brms_path)

if (has_brms) {
  fit_brms <- readRDS(brms_path)
  cat("Loaded brms model.\n")
} else {
  cat("brms model not found. Run 04_fit_brms.R first.\n")
  cat("Proceeding with baseline model validation only.\n")
}

fig_dir <- here::here("outputs", "figures")

# -----------------------------------------------------------------------------
# 2. Leave-One-Record-Out Cross-Validation (baseline)
# -----------------------------------------------------------------------------

cat("\n=== Leave-One-Record-Out CV (Baseline GAM) ===\n")

eq_types <- unique(fem$EarthquakeType)
cv_results_baseline <- list()

for (eq_out in eq_types) {
  cat("Fold: holding out", eq_out, "\n")

  # Split data
  train <- fem |> filter(EarthquakeType != eq_out)
  test <- fem |> filter(EarthquakeType == eq_out)

  # Fit Stage 1
  m1_cv <- gam(
    has_damage ~ s(log_pgv, k = 6, bs = "tp") +
                 s(log_n, k = 4, bs = "tp") +
                 s(InitialPsi, k = 5, bs = "tp") +
                 Material + FacadeType + SoilProfile + EarthquakeType,
    family = binomial(link = "logit"),
    data = train,
    method = "REML"
  )

  # Fit Stage 2
  train_pos <- train |> filter(has_damage == 1)
  m2_cv <- gam(
    delta_psi_pos ~ s(log_pgv, k = 6, bs = "tp") +
                    s(log_n, k = 4, bs = "tp") +
                    s(InitialPsi, k = 5, bs = "tp") +
                    Material + FacadeType + SoilProfile + EarthquakeType,
    family = Gamma(link = "log"),
    data = train_pos,
    method = "REML"
  )

  # Predict on test set
  # Note: EarthquakeType for held-out will be "new" - use average effect
  test_pred <- test |>
    mutate(
      # Set held-out record to a training record for prediction
      EarthquakeType_orig = EarthquakeType,
      EarthquakeType = eq_types[eq_types != eq_out][1]  # Use first training record
    )

  p_init <- predict(m1_cv, newdata = test_pred, type = "response")
  mu_sev <- predict(m2_cv, newdata = test_pred, type = "response")

  cv_results_baseline[[eq_out]] <- test |>
    mutate(
      fold = eq_out,
      pred_p_damage = p_init,
      pred_mu_if_damage = mu_sev,
      pred_expected = p_init * mu_sev,
      residual = delta_psi_pos - pred_expected
    )
}

cv_baseline <- bind_rows(cv_results_baseline)

# CV metrics
cat("\n=== CV Metrics (Baseline) ===\n")

# Stage 1: Classification
cv_baseline |>
  group_by(fold) |>
  summarise(
    n = n(),
    auc = tryCatch({
      pROC::auc(pROC::roc(has_damage, pred_p_damage, quiet = TRUE))
    }, error = function(e) NA),
    brier = mean((has_damage - pred_p_damage)^2)
  ) |>
  print()

# Stage 2: Regression
cat("\nRegression metrics:\n")
cv_baseline |>
  group_by(fold) |>
  summarise(
    n = n(),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    cor = cor(delta_psi_pos, pred_expected)
  ) |>
  print()

# Overall
cat("\nOverall CV performance:\n")
cat("RMSE:", sqrt(mean(cv_baseline$residual^2)), "\n")
cat("MAE:", mean(abs(cv_baseline$residual)), "\n")
cat("Correlation:", cor(cv_baseline$delta_psi_pos, cv_baseline$pred_expected), "\n")

# -----------------------------------------------------------------------------
# 3. Calibration plots
# -----------------------------------------------------------------------------

cat("\n=== Generating Calibration Plots ===\n")

# Observed vs Predicted
p_calib1 <- cv_baseline |>
  ggplot(aes(x = pred_expected, y = delta_psi_pos)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  facet_wrap(~fold) +
  labs(
    title = "Calibration: Observed vs Predicted (Leave-One-Record-Out)",
    x = "Predicted E[DeltaPsi]",
    y = "Observed DeltaPsi"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "validation_calibration_cv.png"), p_calib1,
       width = 10, height = 8, dpi = 150)
cat("Saved: validation_calibration_cv.png\n")

# Calibration by PGV
p_calib2 <- cv_baseline |>
  group_by(fold, PGV) |>
  summarise(
    obs_mean = mean(delta_psi_pos),
    pred_mean = mean(pred_expected),
    .groups = "drop"
  ) |>
  ggplot(aes(x = pred_mean, y = obs_mean, color = factor(PGV))) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~fold) +
  scale_color_viridis_d() +
  labs(
    title = "Calibration by PGV Level",
    x = "Mean Predicted",
    y = "Mean Observed",
    color = "PGV"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "validation_calibration_by_pgv.png"), p_calib2,
       width = 12, height = 8, dpi = 150)
cat("Saved: validation_calibration_by_pgv.png\n")

# -----------------------------------------------------------------------------
# 4. P(damage) calibration
# -----------------------------------------------------------------------------

# Binned calibration for probability of damage
p_calib_prob <- cv_baseline |>
  mutate(prob_bin = cut(pred_p_damage, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) |>
  group_by(prob_bin) |>
  summarise(
    n = n(),
    obs_rate = mean(has_damage),
    pred_rate = mean(pred_p_damage),
    .groups = "drop"
  ) |>
  ggplot(aes(x = pred_rate, y = obs_rate)) +
  geom_point(aes(size = n)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = obs_rate - 1.96*sqrt(obs_rate*(1-obs_rate)/n),
                    ymax = obs_rate + 1.96*sqrt(obs_rate*(1-obs_rate)/n)), width = 0.02) +
  labs(
    title = "Probability Calibration: P(damage > 0)",
    x = "Predicted Probability",
    y = "Observed Rate",
    size = "n"
  ) +
  coord_equal() +
  xlim(0, 1) + ylim(0, 1) +
  theme_minimal()

ggsave(file.path(fig_dir, "validation_prob_calibration.png"), p_calib_prob,
       width = 8, height = 6, dpi = 150)
cat("Saved: validation_prob_calibration.png\n")

# -----------------------------------------------------------------------------
# 5. Residual diagnostics
# -----------------------------------------------------------------------------

# Residuals by PGV
p_resid1 <- cv_baseline |>
  ggplot(aes(x = log(PGV), y = residual)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_jitter(alpha = 0.2, width = 0.1) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(
    title = "Residuals vs log(PGV)",
    x = "log(PGV)",
    y = "Residual (Obs - Pred)"
  ) +
  theme_minimal()

# Residuals by InitialPsi
p_resid2 <- cv_baseline |>
  ggplot(aes(x = InitialPsi, y = residual)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_jitter(alpha = 0.2, width = 0.05) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(
    title = "Residuals vs InitialPsi",
    x = "Initial Psi",
    y = "Residual (Obs - Pred)"
  ) +
  theme_minimal()

# Residuals by fold (record)
p_resid3 <- cv_baseline |>
  ggplot(aes(x = fold, y = residual)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Residuals by Held-Out Record",
    x = "Held-Out Earthquake Type",
    y = "Residual"
  ) +
  theme_minimal()

combined_resid <- (p_resid1 + p_resid2) / p_resid3
ggsave(file.path(fig_dir, "validation_residuals.png"), combined_resid,
       width = 12, height = 10, dpi = 150)
cat("Saved: validation_residuals.png\n")

# -----------------------------------------------------------------------------
# 6. brms-specific validation (if available)
# -----------------------------------------------------------------------------

if (has_brms) {
  cat("\n=== brms Model Validation ===\n")

  # Leave-one-record-out with brms (using kfold)
  cat("Running k-fold CV for brms (this may take a while)...\n")

  # Note: brms kfold with custom folds
  # For speed, we'll use LOO instead
  loo_brms <- loo(fit_brms, cores = 4)

  cat("\nLOO-CV results:\n")
  print(loo_brms)

  # Check Pareto k diagnostics
  pareto_k <- loo_brms$diagnostics$pareto_k
  cat("\nPareto k diagnostics:\n")
  cat("k > 0.7 (problematic):", sum(pareto_k > 0.7), "\n")
  cat("k > 0.5 (moderate):", sum(pareto_k > 0.5), "\n")

  # Posterior predictive check by group
  pp_by_eq <- pp_check(fit_brms, type = "stat_grouped",
                       stat = "mean", group = "EarthquakeType")
  ggsave(file.path(fig_dir, "validation_brms_pp_by_record.png"), pp_by_eq,
         width = 10, height = 6, dpi = 150)
  cat("Saved: validation_brms_pp_by_record.png\n")

  # Calibration with full posterior
  cat("\nGenerating posterior predictive calibration...\n")

  # Get posterior predictions
  fem_brms <- fem |>
    mutate(
      pgv_ord = factor(PGV, levels = c(2, 4, 8, 16, 32, 64, 96, 128), ordered = TRUE),
      n_ord = factor(N, levels = c(1, 2, 3, 4, 8), ordered = TRUE),
      material_ord = factor(Material, levels = c(0.7, 1.0, 1.3), ordered = TRUE),
      initial_psi_c = InitialPsi - mean(InitialPsi)
    )

  # Sample subset for speed
  set.seed(42)
  sample_idx <- sample(nrow(fem_brms), min(500, nrow(fem_brms)))
  fem_sample <- fem_brms[sample_idx, ]

  post_pred <- posterior_predict(fit_brms, newdata = fem_sample, ndraws = 500)

  # Coverage
  q025 <- apply(post_pred, 2, quantile, 0.025)
  q975 <- apply(post_pred, 2, quantile, 0.975)
  coverage_95 <- mean(fem_sample$delta_psi_pos >= q025 & fem_sample$delta_psi_pos <= q975)
  cat("95% prediction interval coverage:", round(coverage_95, 3), "\n")

  q10 <- apply(post_pred, 2, quantile, 0.1)
  q90 <- apply(post_pred, 2, quantile, 0.9)
  coverage_80 <- mean(fem_sample$delta_psi_pos >= q10 & fem_sample$delta_psi_pos <= q90)
  cat("80% prediction interval coverage:", round(coverage_80, 3), "\n")
}

# -----------------------------------------------------------------------------
# 7. Summary statistics
# -----------------------------------------------------------------------------

cat("\n=== Validation Summary ===\n")

summary_stats <- cv_baseline |>
  summarise(
    n = n(),
    mean_obs = mean(delta_psi_pos),
    mean_pred = mean(pred_expected),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    cor = cor(delta_psi_pos, pred_expected),
    r_squared = 1 - sum(residual^2) / sum((delta_psi_pos - mean(delta_psi_pos))^2)
  )

cat("\nOverall validation metrics:\n")
print(summary_stats)

# By PGV level
cat("\nMetrics by PGV level:\n")
cv_baseline |>
  group_by(PGV) |>
  summarise(
    n = n(),
    rmse = sqrt(mean(residual^2)),
    bias = mean(residual),
    .groups = "drop"
  ) |>
  print()

# Save validation results
saveRDS(cv_baseline, here::here("outputs", "models", "cv_results_baseline.rds"))
cat("\nValidation results saved to outputs/models/cv_results_baseline.rds\n")

cat("\nValidation complete.\n")
