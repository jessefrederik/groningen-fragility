# =============================================================================
# 02_explore_data.R
#
# Exploratory data analysis for Korswagen FEM fragility data.
# Visualizes relationships and checks distributional assumptions.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(patchwork)  # For combining plots

# -----------------------------------------------------------------------------
# 1. Load prepared data
# -----------------------------------------------------------------------------

fem <- readRDS(here::here("outputs", "models", "fem_prepared.rds"))

# Create output directory for figures
fig_dir <- here::here("outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. DeltaPsi vs PGV relationship
# -----------------------------------------------------------------------------

# Main relationship: DeltaPsi vs PGV
p1 <- fem |>
  ggplot(aes(x = PGV, y = delta_psi_pos)) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.1, height = 0)) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
  scale_x_log10() +
  labs(
    title = "Damage Increment vs Peak Ground Velocity",
    x = "PGV (mm/s, log scale)",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# By Material strength
p2 <- fem |>
  mutate(Material = factor(Material, labels = c("Weak (0.7)", "Standard (1.0)", "Strong (1.3)"))) |>
  ggplot(aes(x = PGV, y = delta_psi_pos, color = Material)) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.1, height = 0)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_log10() +
  scale_color_viridis_d(option = "C") +
  labs(
    title = "By Material Strength",
    x = "PGV (mm/s, log scale)",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# By Initial damage
p3 <- fem |>
  mutate(InitialDamage = factor(DesignInitialPsi,
    labels = c("None (0)", "Slight (0.5)", "Visible (1.0)", "Moderate (1.5)"))) |>
  ggplot(aes(x = PGV, y = delta_psi_pos, color = InitialDamage)) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.1, height = 0)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_log10() +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "By Initial Damage State",
    x = "PGV (mm/s, log scale)",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# By N (number of events)
p4 <- fem |>
  mutate(N = factor(N)) |>
  ggplot(aes(x = PGV, y = delta_psi_pos, color = N)) +
  geom_point(alpha = 0.2, position = position_jitter(width = 0.1, height = 0)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_log10() +
  scale_color_viridis_d() +
  labs(
    title = "By Number of Events (N)",
    x = "PGV (mm/s, log scale)",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# Combine plots
combined_pgv <- (p1 + p2) / (p3 + p4)
ggsave(file.path(fig_dir, "eda_deltapsi_vs_pgv.png"), combined_pgv,
       width = 12, height = 10, dpi = 150)
cat("Saved: eda_deltapsi_vs_pgv.png\n")

# -----------------------------------------------------------------------------
# 3. DeltaPsi by categorical factors
# -----------------------------------------------------------------------------

# By FacadeType and SoilProfile
p5 <- fem |>
  ggplot(aes(x = FacadeType, y = delta_psi_pos, fill = SoilProfile)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Damage by Facade Type and Soil Profile",
    x = "Facade Type",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# By EarthquakeType (record-to-record variability)
p6 <- fem |>
  ggplot(aes(x = EarthquakeType, y = delta_psi_pos)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Damage by Earthquake Record",
    subtitle = "Record-to-record variability (only 4 records)",
    x = "Earthquake Type",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

# Faceted by all categorical
p7 <- fem |>
  ggplot(aes(x = PGV, y = delta_psi_pos)) +
  geom_jitter(alpha = 0.2, width = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  scale_x_log10() +
  facet_grid(FacadeType ~ SoilProfile, labeller = label_both) +
  labs(
    title = "Damage vs PGV by Facade and Soil",
    x = "PGV (mm/s, log scale)",
    y = expression(Delta * Psi)
  ) +
  theme_minimal()

combined_cat <- (p5 + p6) / p7
ggsave(file.path(fig_dir, "eda_deltapsi_by_factors.png"), combined_cat,
       width = 12, height = 10, dpi = 150)
cat("Saved: eda_deltapsi_by_factors.png\n")

# -----------------------------------------------------------------------------
# 4. Distribution of DeltaPsi (check Gamma assumption)
# -----------------------------------------------------------------------------

# Overall distribution
p8 <- fem |>
  filter(has_damage == 1) |>  # Only positive values
  ggplot(aes(x = delta_psi_pos)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(
    title = expression("Distribution of " * Delta * Psi * " | " * Delta * Psi > 0),
    subtitle = "Right-skewed, consistent with Gamma distribution",
    x = expression(Delta * Psi),
    y = "Density"
  ) +
  theme_minimal()

# Log-scale to check for log-normal
p9 <- fem |>
  filter(has_damage == 1) |>
  ggplot(aes(x = log(delta_psi_pos))) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "coral", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  stat_function(fun = dnorm,
                args = list(mean = mean(log(fem$delta_psi_pos[fem$has_damage == 1])),
                            sd = sd(log(fem$delta_psi_pos[fem$has_damage == 1]))),
                color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = expression("Distribution of log(" * Delta * Psi * ") | " * Delta * Psi > 0),
    subtitle = "Dashed blue: fitted Normal (log-normal assumption)",
    x = expression(log(Delta * Psi)),
    y = "Density"
  ) +
  theme_minimal()

# By PGV level
p10 <- fem |>
  filter(has_damage == 1) |>
  mutate(pgv_group = cut(PGV, breaks = c(0, 8, 32, 128), labels = c("Low (2-8)", "Mid (16-32)", "High (64-128)"))) |>
  ggplot(aes(x = delta_psi_pos, fill = pgv_group)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(
    title = expression("Distribution of " * Delta * Psi * " by PGV Level"),
    x = expression(Delta * Psi),
    y = "Density",
    fill = "PGV Group"
  ) +
  theme_minimal()

# QQ plot for Gamma
pos_data <- fem$delta_psi_pos[fem$has_damage == 1]
shape_est <- mean(pos_data)^2 / var(pos_data)
rate_est <- mean(pos_data) / var(pos_data)

p11 <- tibble(observed = sort(pos_data)) |>
  mutate(
    theoretical = qgamma(ppoints(n()), shape = shape_est, rate = rate_est)
  ) |>
  ggplot(aes(x = theoretical, y = observed)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Q-Q Plot: Gamma Distribution",
    subtitle = sprintf("Shape = %.2f, Rate = %.2f", shape_est, rate_est),
    x = "Theoretical Gamma Quantiles",
    y = "Observed Quantiles"
  ) +
  theme_minimal()

combined_dist <- (p8 + p9) / (p10 + p11)
ggsave(file.path(fig_dir, "eda_deltapsi_distribution.png"), combined_dist,
       width = 12, height = 10, dpi = 150)
cat("Saved: eda_deltapsi_distribution.png\n")

# -----------------------------------------------------------------------------
# 5. Probability of damage (hurdle component)
# -----------------------------------------------------------------------------

# P(damage) vs PGV
p12 <- fem |>
  group_by(PGV) |>
  summarise(
    p_damage = mean(has_damage),
    n = n(),
    se = sqrt(p_damage * (1 - p_damage) / n)
  ) |>
  ggplot(aes(x = PGV, y = p_damage)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0, p_damage - 1.96*se),
                    ymax = pmin(1, p_damage + 1.96*se)), width = 0.1) +
  geom_line(linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Probability of Any Damage by PGV",
    x = "PGV (mm/s, log scale)",
    y = expression(P(Delta * Psi > 0))
  ) +
  theme_minimal()

# P(damage) by InitialPsi
p13 <- fem |>
  group_by(DesignInitialPsi, PGV) |>
  summarise(p_damage = mean(has_damage), .groups = "drop") |>
  mutate(InitialDamage = factor(DesignInitialPsi)) |>
  ggplot(aes(x = PGV, y = p_damage, color = InitialDamage)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "P(Damage) by Initial Damage State",
    x = "PGV (mm/s, log scale)",
    y = expression(P(Delta * Psi > 0)),
    color = expression(Psi[0])
  ) +
  theme_minimal()

combined_hurdle <- p12 + p13
ggsave(file.path(fig_dir, "eda_damage_probability.png"), combined_hurdle,
       width = 12, height = 5, dpi = 150)
cat("Saved: eda_damage_probability.png\n")

# -----------------------------------------------------------------------------
# 6. Check for truncation (missing low-PGV observations)
# -----------------------------------------------------------------------------

# Heatmap of observation counts
p14 <- fem |>
  count(PGV, N, DesignInitialPsi) |>
  ggplot(aes(x = factor(PGV), y = factor(N), fill = n)) +
  geom_tile() +
  facet_wrap(~DesignInitialPsi, labeller = label_both) +
  scale_fill_viridis_c(option = "B") +
  labs(
    title = "Observation Counts by PGV, N, and Initial Damage",
    subtitle = "Missing cells indicate potential left-truncation",
    x = "PGV (mm/s)",
    y = "N (events)",
    fill = "Count"
  ) +
  theme_minimal()

ggsave(file.path(fig_dir, "eda_truncation_pattern.png"), p14,
       width = 12, height = 8, dpi = 150)
cat("Saved: eda_truncation_pattern.png\n")

# -----------------------------------------------------------------------------
# 7. Summary statistics by group
# -----------------------------------------------------------------------------

cat("\n=== Mean DeltaPsi by PGV ===\n")
fem |>
  group_by(PGV) |>
  summarise(
    n = n(),
    mean_delta = mean(delta_psi_pos),
    sd_delta = sd(delta_psi_pos),
    p_damage = mean(has_damage),
    mean_if_damage = mean(delta_psi_pos[has_damage == 1])
  ) |>
  print()

cat("\n=== Mean DeltaPsi by Material ===\n")
fem |>
  group_by(Material) |>
  summarise(
    n = n(),
    mean_delta = mean(delta_psi_pos),
    p_damage = mean(has_damage)
  ) |>
  print()

cat("\n=== Mean DeltaPsi by EarthquakeType ===\n")
fem |>
  group_by(EarthquakeType) |>
  summarise(
    n = n(),
    mean_delta = mean(delta_psi_pos),
    sd_delta = sd(delta_psi_pos),
    p_damage = mean(has_damage)
  ) |>
  print()

cat("\nExploratory analysis complete. Check outputs/figures/ for plots.\n")
