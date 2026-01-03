# =============================================================================
# 04_fit_brms.R
#
# Fit hierarchical Bayesian hurdle-Gamma model using brms.
# Includes monotonic effects for PGV and partial pooling across records.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(brms)
library(bayesplot)

# Set Stan options for better performance
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# -----------------------------------------------------------------------------
# 1. Load prepared data
# -----------------------------------------------------------------------------

fem <- readRDS(here::here("outputs", "models", "fem_prepared.rds"))

cat("Data loaded:", nrow(fem), "observations\n")

# -----------------------------------------------------------------------------
# 2. Prepare data for brms monotonic effects
# -----------------------------------------------------------------------------

# For monotonic effects, brms uses ordered factors via mo()
# We'll create an ordered PGV category

# Check current PGV levels
cat("PGV levels:", sort(unique(fem$PGV)), "\n")

# Create ordered factor for monotonic effect
fem_brms <- fem |>
  mutate(
    # Ordered factor for monotonic PGV effect
    pgv_ord = factor(PGV, levels = c(2, 4, 8, 16, 32, 64, 96, 128), ordered = TRUE),

    # N as ordered factor for monotonic effect
    n_ord = factor(N, levels = c(1, 2, 3, 4, 8), ordered = TRUE),

    # Material as ordered factor (higher = stronger = less damage)
    material_ord = factor(Material, levels = c(0.7, 1.0, 1.3), ordered = TRUE),

    # Keep continuous versions for alternative formulation
    log_pgv_c = log_pgv - mean(log_pgv),  # centered
    log_n_c = log_n - mean(log_n),
    initial_psi_c = InitialPsi - mean(InitialPsi),
    material_c = Material - mean(Material)
  )

# -----------------------------------------------------------------------------
# 3. Define prior distributions
# -----------------------------------------------------------------------------

# Priors informed by:
# - Physical constraints (PGV effect should be positive)
# - Scale of the data (DeltaPsi typically 0-3)
# - Korswagen paper's findings (σ_record ~ 0.2-0.4)

priors <- c(
  # Fixed effects (weakly informative)
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 2), class = "Intercept"),

  # Monotonic simplex prior (default Dirichlet(1,...,1) is uniform)
  # We keep default for now

  # Hierarchical SD priors
  prior(exponential(2), class = "sd"),  # Expect small random effects

  # Hurdle probability intercept
  prior(normal(0, 2), class = "Intercept", dpar = "hu"),

  # Gamma shape parameter
  prior(gamma(2, 0.5), class = "shape")  # Expect shape ~ 1-4
)

# -----------------------------------------------------------------------------
# 4. Fit Model A: Full hierarchical hurdle-Gamma
# -----------------------------------------------------------------------------

cat("\n=== Fitting Model A: Hierarchical Hurdle-Gamma ===\n")
cat("This may take 15-30 minutes...\n")

# Formula with monotonic effects and hierarchical structure
# mo() = monotonic effect for ordered factors
# (1|group) = random intercept

formula_A <- bf(
  # Gamma component (mean)
  delta_psi_pos ~ mo(pgv_ord) + mo(n_ord) + mo(material_ord) + initial_psi_c +
                  (1 | FacadeType:SoilProfile) + (1 | EarthquakeType),
  # Hurdle component (P(y=0))
  hu ~ mo(pgv_ord) + mo(n_ord) + mo(material_ord) + initial_psi_c +
       (1 | FacadeType:SoilProfile) + (1 | EarthquakeType),
  family = hurdle_gamma()
)

fit_A <- brm(
  formula_A,
  data = fem_brms,
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = here::here("outputs", "models", "brms_hurdle_gamma_A")
)

cat("\nModel A fitted.\n")

# -----------------------------------------------------------------------------
# 5. Model summary and diagnostics
# -----------------------------------------------------------------------------

cat("\n=== Model A Summary ===\n")
print(summary(fit_A))

# Check convergence
cat("\n=== Convergence Diagnostics ===\n")
rhat_vals <- rhat(fit_A)
cat("Max Rhat:", max(rhat_vals, na.rm = TRUE), "\n")
cat("Any Rhat > 1.01:", any(rhat_vals > 1.01, na.rm = TRUE), "\n")

ess_vals <- neff_ratio(fit_A)
cat("Min ESS ratio:", min(ess_vals, na.rm = TRUE), "\n")

# Check for divergences
np <- nuts_params(fit_A)
divergent <- sum(subset(np, Parameter == "divergent__")$Value)
cat("Divergent transitions:", divergent, "\n")

# -----------------------------------------------------------------------------
# 6. Posterior predictive checks
# -----------------------------------------------------------------------------

cat("\n=== Generating Posterior Predictive Checks ===\n")

fig_dir <- here::here("outputs", "figures")

# Overall density overlay
pp_dens <- pp_check(fit_A, type = "dens_overlay", ndraws = 100)
ggsave(file.path(fig_dir, "brms_pp_check_density.png"), pp_dens,
       width = 8, height = 5, dpi = 150)
cat("Saved: brms_pp_check_density.png\n")

# By groups
pp_grouped <- pp_check(fit_A, type = "stat_grouped", stat = "mean", group = "pgv_ord")
ggsave(file.path(fig_dir, "brms_pp_check_by_pgv.png"), pp_grouped,
       width = 10, height = 6, dpi = 150)
cat("Saved: brms_pp_check_by_pgv.png\n")

# Proportion of zeros
pp_zeros <- pp_check(fit_A, type = "stat", stat = function(y) mean(y == 0))
ggsave(file.path(fig_dir, "brms_pp_check_zeros.png"), pp_zeros,
       width = 8, height = 5, dpi = 150)
cat("Saved: brms_pp_check_zeros.png\n")

# -----------------------------------------------------------------------------
# 7. Extract and visualize effects
# -----------------------------------------------------------------------------

cat("\n=== Extracting Effects ===\n")

# Conditional effects for PGV
ce_pgv <- conditional_effects(fit_A, effects = "pgv_ord")
p_ce_pgv <- plot(ce_pgv, ask = FALSE)[[1]] +
  labs(
    title = "Marginal Effect of PGV on Expected Damage",
    x = "PGV (mm/s)",
    y = expression(E[Delta * Psi])
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "brms_effect_pgv.png"), p_ce_pgv,
       width = 8, height = 5, dpi = 150)
cat("Saved: brms_effect_pgv.png\n")

# Conditional effects for Material
ce_mat <- conditional_effects(fit_A, effects = "material_ord")
p_ce_mat <- plot(ce_mat, ask = FALSE)[[1]] +
  labs(
    title = "Marginal Effect of Material Strength on Expected Damage",
    x = "Material Strength",
    y = expression(E[Delta * Psi])
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "brms_effect_material.png"), p_ce_mat,
       width = 8, height = 5, dpi = 150)
cat("Saved: brms_effect_material.png\n")

# -----------------------------------------------------------------------------
# 8. Extract random effects (record-to-record variability)
# -----------------------------------------------------------------------------

cat("\n=== Random Effects (Record-to-Record Variability) ===\n")

# Get random effect estimates
ranef_est <- ranef(fit_A)

cat("\nEarthquakeType random effects:\n")
print(ranef_est$EarthquakeType)

cat("\nFacadeType:SoilProfile random effects:\n")
print(ranef_est$`FacadeType:SoilProfile`)

# SD of random effects
cat("\nRandom effect standard deviations:\n")
posterior_summary(fit_A, pars = "sd_")

# -----------------------------------------------------------------------------
# 9. Alternative Model B: Continuous log(PGV) with spline
# -----------------------------------------------------------------------------

cat("\n=== Fitting Model B: Continuous spline (for comparison) ===\n")

# Alternative with smooth term (less interpretable but more flexible)
# Note: brms doesn't enforce monotonicity on splines, so we use this for comparison

formula_B <- bf(
  delta_psi_pos ~ s(log_pgv_c, k = 5) + log_n_c + material_c + initial_psi_c +
                  (1 | FacadeType:SoilProfile) + (1 | EarthquakeType),
  hu ~ s(log_pgv_c, k = 5) + log_n_c + material_c + initial_psi_c +
       (1 | FacadeType:SoilProfile) + (1 | EarthquakeType),
  family = hurdle_gamma()
)

fit_B <- brm(
  formula_B,
  data = fem_brms,
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 42,
  control = list(adapt_delta = 0.95),
  file = here::here("outputs", "models", "brms_hurdle_gamma_B")
)

cat("\nModel B fitted.\n")

# -----------------------------------------------------------------------------
# 10. Model comparison
# -----------------------------------------------------------------------------

cat("\n=== Model Comparison (LOO-CV) ===\n")

loo_A <- loo(fit_A, cores = 4)
loo_B <- loo(fit_B, cores = 4)

cat("\nModel A (monotonic):\n")
print(loo_A)

cat("\nModel B (spline):\n")
print(loo_B)

cat("\nLOO comparison:\n")
print(loo_compare(loo_A, loo_B))

# -----------------------------------------------------------------------------
# 11. Create prediction function
# -----------------------------------------------------------------------------

#' Predict damage from brms hurdle-gamma model
#'
#' @param newdata Data frame with predictors
#' @param model Fitted brms model
#' @param n_sims Number of posterior draws to use
#' @return Data frame with predictions and uncertainty
predict_brms_hurdle <- function(newdata, model = fit_A, n_sims = 1000) {
  # Get posterior predictions
  post_pred <- posterior_predict(model, newdata = newdata, ndraws = n_sims)

  # Summary statistics
  result <- newdata |>
    mutate(
      mean_damage = colMeans(post_pred),
      median_damage = apply(post_pred, 2, median),
      sd_damage = apply(post_pred, 2, sd),
      q05_damage = apply(post_pred, 2, quantile, 0.05),
      q25_damage = apply(post_pred, 2, quantile, 0.25),
      q75_damage = apply(post_pred, 2, quantile, 0.75),
      q95_damage = apply(post_pred, 2, quantile, 0.95),
      p_no_damage = apply(post_pred, 2, function(x) mean(x == 0))
    )

  return(result)
}

# Save prediction function and fitted model reference
save(predict_brms_hurdle, file = here::here("outputs", "models", "predict_brms_hurdle.RData"))

# -----------------------------------------------------------------------------
# 12. Summary of key findings
# -----------------------------------------------------------------------------

cat("\n=== Key Findings ===\n")

# Extract fixed effects
fixef_A <- fixef(fit_A)
cat("\nFixed effects (Model A):\n")
print(round(fixef_A, 3))

# Monotonic effect increments for PGV
cat("\nMonotonic increments for PGV (Model A, Gamma component):\n")
simo_pgv <- as_draws_df(fit_A) |>
  select(starts_with("simo_mopgv_ord")) |>
  summarise(across(everything(), list(mean = mean, sd = sd)))
print(simo_pgv)

# Gamma shape
cat("\nGamma shape parameter:\n")
print(posterior_summary(fit_A, pars = "shape"))

cat("\nBayesian model fitting complete.\n")
cat("Primary model saved to: outputs/models/brms_hurdle_gamma_A.rds\n")
