# =============================================================================
# 07_building_mapping.R
#
# Map building covariates (bouwjaar, oppervlakte, postcode) to FEM parameters.
# Creates heuristic priors for vulnerability parameters.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)
library(sf)
library(truncnorm)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load building and site data
# -----------------------------------------------------------------------------

cat("Loading building data...\n")
buildings <- st_read(
  here::here("datafiles", "bag_selectie.gpkg"),
  quiet = TRUE
)
cat("Loaded", nrow(buildings), "buildings\n")

# Load Vs30 data
vs30 <- read_delim(
  here::here("datafiles", "vs30.csv"),
  delim = ";",
  locale = locale(decimal_mark = ","),
  col_names = c("postcode", "vs30", "x1", "x2", "x3"),
  skip = 1,
  show_col_types = FALSE
) |>
  select(postcode = postcode, vs30 = vs30) |>
  mutate(postcode = as.character(postcode))

cat("Loaded Vs30 for", nrow(vs30), "postcodes\n")

# -----------------------------------------------------------------------------
# 2. Explore building characteristics
# -----------------------------------------------------------------------------

cat("\n=== Building Characteristics ===\n")

# Bouwjaar distribution
cat("\nBouwjaar (construction year) summary:\n")
print(summary(buildings$bouwjaar))

# Oppervlakte distribution
cat("\nOppervlakte (floor area) summary:\n")
print(summary(buildings$oppervlakte))

# Gebruiksdoel
cat("\nGebruiksdoel (use type) counts:\n")
print(table(buildings$gebruiksdoel))

# -----------------------------------------------------------------------------
# 3. Define heuristic mapping functions
# -----------------------------------------------------------------------------

#' Sample material strength based on building age
#'
#' Older buildings may have weaker masonry due to:
#' - Degradation over time
#' - Different construction standards
#' - Mortar weathering
#'
#' Based on Korswagen: m ~ N(1, 0.3) truncated to [0.7, 1.3]
#'
#' Heuristic: shift mean down for older buildings
#'
#' @param bouwjaar Construction year
#' @param n Number of samples
#' @return Vector of material strength samples
sample_material_strength <- function(bouwjaar, n = 1) {
  # Base parameters from Korswagen
  base_mean <- 1.0
  base_sd <- 0.3

  # Age-based adjustment (heuristic)
  # Pre-1920: lower mean (older mortar, more degradation)
  # 1920-1970: slight reduction
  # Post-1970: standard

  age_adjustment <- case_when(
    bouwjaar < 1920 ~ -0.15,
    bouwjaar < 1970 ~ -0.05,
    TRUE ~ 0
  )

  adjusted_mean <- base_mean + age_adjustment

  # Sample from truncated normal
  rtruncnorm(n, a = 0.7, b = 1.3, mean = adjusted_mean, sd = base_sd)
}

#' Assign facade type probability based on building characteristics (vectorized)
#'
#' Korswagen facade types:
#' - A: Flexible facade with opening (window), more susceptible to cracking
#' - B: Stiffer shear wall, less susceptible
#'
#' Heuristic:
#' - Larger buildings more likely to have structural walls (B)
#' - Residential typically has windows (A)
#'
#' @param oppervlakte Floor area (m2)
#' @param gebruiksdoel Use type
#' @return Probability of FacadeType = "A"
prob_facade_A <- function(oppervlakte, gebruiksdoel) {
  # Base probability
  p_A <- rep(0.5, length(oppervlakte))

  # Residential buildings more likely to have window-heavy facades
  is_woon <- !is.na(gebruiksdoel) & grepl("woon", gebruiksdoel, ignore.case = TRUE)
  p_A <- p_A + ifelse(is_woon, 0.2, 0)

  # Larger buildings slightly more likely to have structural walls
  p_A <- p_A + case_when(
    is.na(oppervlakte) ~ 0,
    oppervlakte > 200 ~ -0.1,
    oppervlakte < 80 ~ 0.1,
    TRUE ~ 0
  )

  # Clamp to [0.2, 0.8] - always some uncertainty
  pmax(0.2, pmin(0.8, p_A))
}

#' Assign soil profile based on Vs30 (vectorized)
#'
#' Korswagen soil profiles:
#' - A: Stiffer soil (sand-dominated)
#' - B: Softer soil (peat/clay-dominated)
#'
#' Typical Groningen:
#' - Lower Vs30 (< 200 m/s) indicates softer soil
#'
#' @param vs30 Shear wave velocity (m/s)
#' @return "A" or "B"
assign_soil_profile <- function(vs30) {
  n <- length(vs30)
  # Probabilistic assignment for intermediate values
  random_vals <- runif(n)
  prob_A <- (vs30 - 180) / 70  # Probability of A in intermediate range

  result <- case_when(
    is.na(vs30) ~ "A",           # Default to stiffer if unknown
    vs30 < 180 ~ "B",            # Soft soil
    vs30 > 250 ~ "A",            # Stiff soil
    random_vals < prob_A ~ "A",  # Intermediate - probabilistic
    TRUE ~ "B"
  )
  result
}

#' Sample initial damage state (Psi_0)
#'
#' Initial damage highly uncertain. Korswagen uses discrete scenarios.
#' Heuristic: older buildings more likely to have pre-existing damage
#' from settlement, thermal cycling, etc.
#'
#' @param bouwjaar Construction year
#' @param n Number of samples
#' @return Vector of initial Psi values
sample_initial_psi <- function(bouwjaar, n = 1) {
  # Discrete scenarios from Korswagen
  psi_levels <- c(0, 0.5, 1.0, 1.5)

  # Age-based probabilities
  # Older buildings more likely to have pre-damage
  probs <- case_when(
    bouwjaar < 1920 ~ c(0.2, 0.3, 0.3, 0.2),  # High prior damage
    bouwjaar < 1950 ~ c(0.3, 0.35, 0.25, 0.1),
    bouwjaar < 1980 ~ c(0.5, 0.3, 0.15, 0.05),
    TRUE ~ c(0.7, 0.2, 0.08, 0.02)  # Newer: mostly undamaged
  )

  sample(psi_levels, n, replace = TRUE, prob = probs)
}

# -----------------------------------------------------------------------------
# 4. Apply mappings to building dataset
# -----------------------------------------------------------------------------

cat("\n=== Applying Heuristic Mappings ===\n")

# Extract PC4 from postcode
buildings <- buildings |>
  mutate(
    pc4 = substr(postcode, 1, 4)
  )

# Join Vs30
buildings <- buildings |>
  left_join(vs30, by = c("pc4" = "postcode"))

cat("Vs30 coverage:", sum(!is.na(buildings$vs30)), "/", nrow(buildings), "\n")

# For buildings without Vs30, use regional average
mean_vs30 <- mean(buildings$vs30, na.rm = TRUE)
buildings <- buildings |>
  mutate(vs30 = if_else(is.na(vs30), mean_vs30, vs30))

# Assign soil profile
set.seed(42)
buildings <- buildings |>
  mutate(
    soil_profile = assign_soil_profile(vs30),
    prob_facade_a = prob_facade_A(oppervlakte, gebruiksdoel)
  )

# Summary of mappings
cat("\nSoil profile distribution:\n")
print(table(buildings$soil_profile))

cat("\nP(Facade=A) summary:\n")
print(summary(buildings$prob_facade_a))

# -----------------------------------------------------------------------------
# 5. Generate vulnerability parameter samples for subset
# -----------------------------------------------------------------------------

cat("\n=== Generating Vulnerability Samples (Example) ===\n")

# Take a sample of buildings for demonstration
set.seed(42)
sample_idx <- sample(nrow(buildings), min(1000, nrow(buildings)))
buildings_sample <- buildings[sample_idx, ]

# Generate vulnerability parameters
n_sims <- 100  # Samples per building

vulnerability_samples <- buildings_sample |>
  st_drop_geometry() |>
  select(identificatie, bouwjaar, oppervlakte, vs30, soil_profile, prob_facade_a) |>
  rowwise() |>
  mutate(
    # Material strength samples
    material_samples = list(sample_material_strength(bouwjaar, n_sims)),
    # Facade type samples
    facade_samples = list(ifelse(runif(n_sims) < prob_facade_a, "A", "B")),
    # Initial Psi samples
    initial_psi_samples = list(sample_initial_psi(bouwjaar, n_sims))
  ) |>
  ungroup()

# Summarize distributions
cat("\nMaterial strength by bouwjaar range:\n")
vulnerability_samples |>
  mutate(
    bouwjaar_cat = cut(bouwjaar,
      breaks = c(0, 1920, 1950, 1980, 2025),
      labels = c("<1920", "1920-50", "1950-80", ">1980"))
  ) |>
  group_by(bouwjaar_cat) |>
  summarise(
    n = n(),
    mean_material = mean(sapply(material_samples, mean)),
    sd_material = mean(sapply(material_samples, sd))
  ) |>
  print()

# -----------------------------------------------------------------------------
# 6. Create mapping functions for simulation pipeline
# -----------------------------------------------------------------------------

#' Generate vulnerability parameters for a building
#'
#' @param bouwjaar Construction year
#' @param oppervlakte Floor area
#' @param gebruiksdoel Use type
#' @param vs30 Shear wave velocity
#' @param n Number of Monte Carlo samples
#' @return Data frame with n rows of sampled parameters
generate_building_vulnerability <- function(
  bouwjaar = 1970,
  oppervlakte = 100,
  gebruiksdoel = "woonfunctie",
  vs30 = 200,
  n = 100
) {
  tibble(
    material = sample_material_strength(bouwjaar, n),
    facade_type = ifelse(
      runif(n) < prob_facade_A(oppervlakte, gebruiksdoel),
      "A", "B"
    ),
    soil_profile = assign_soil_profile(rep(vs30, n)),
    initial_psi = sample_initial_psi(bouwjaar, n)
  )
}

#' Batch generate vulnerability for building data frame
#'
#' @param buildings_df Data frame with building attributes
#' @param n_sims Number of Monte Carlo samples per building
#' @return List of vulnerability parameter matrices
generate_vulnerability_batch <- function(buildings_df, n_sims = 100) {
  n_buildings <- nrow(buildings_df)

  # Pre-allocate matrices
  material <- matrix(NA, nrow = n_buildings, ncol = n_sims)
  facade <- matrix(NA, nrow = n_buildings, ncol = n_sims)
  soil <- matrix(NA, nrow = n_buildings, ncol = n_sims)
  psi0 <- matrix(NA, nrow = n_buildings, ncol = n_sims)

  for (i in 1:n_buildings) {
    row <- buildings_df[i, ]

    material[i, ] <- sample_material_strength(row$bouwjaar, n_sims)
    facade[i, ] <- ifelse(
      runif(n_sims) < prob_facade_A(row$oppervlakte, row$gebruiksdoel),
      "A", "B"
    )
    soil[i, ] <- assign_soil_profile(rep(row$vs30, n_sims))
    psi0[i, ] <- sample_initial_psi(row$bouwjaar, n_sims)
  }

  list(
    material = material,
    facade_type = facade,
    soil_profile = soil,
    initial_psi = psi0
  )
}

# -----------------------------------------------------------------------------
# 7. Visualize heuristic relationships
# -----------------------------------------------------------------------------

cat("\n=== Generating Visualization ===\n")

fig_dir <- here::here("outputs", "figures")

# Material strength vs bouwjaar
p_material <- tibble(
  bouwjaar = rep(seq(1850, 2020, by = 10), each = 100),
  material = unlist(lapply(seq(1850, 2020, by = 10), function(y) sample_material_strength(y, 100)))
) |>
  ggplot(aes(x = bouwjaar, y = material)) +
  geom_jitter(alpha = 0.1, width = 3) +
  geom_smooth(method = "loess", color = "red") +
  labs(
    title = "Material Strength Prior by Construction Year",
    subtitle = "Heuristic: older buildings have lower mean strength",
    x = "Bouwjaar",
    y = "Material Strength (m)"
  ) +
  theme_minimal()

# Initial Psi vs bouwjaar
p_psi0 <- tibble(
  bouwjaar = rep(seq(1850, 2020, by = 10), each = 100),
  psi0 = unlist(lapply(seq(1850, 2020, by = 10), function(y) sample_initial_psi(y, 100)))
) |>
  ggplot(aes(x = bouwjaar, y = psi0)) +
  geom_jitter(alpha = 0.1, width = 3, height = 0.05) +
  geom_smooth(method = "loess", color = "red") +
  labs(
    title = "Initial Damage Prior by Construction Year",
    subtitle = "Heuristic: older buildings more likely to have pre-damage",
    x = "Bouwjaar",
    y = expression(Psi[0])
  ) +
  theme_minimal()

# Soil profile vs Vs30
p_soil <- tibble(
  vs30 = seq(150, 320, by = 1)
) |>
  mutate(
    p_soft = sapply(vs30, function(v) mean(assign_soil_profile(rep(v, 1000)) == "B"))
  ) |>
  ggplot(aes(x = vs30, y = p_soft)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(180, 250), linetype = "dashed", color = "red") +
  labs(
    title = "Soil Profile Assignment by Vs30",
    subtitle = "Soft soil (B) more likely below 180 m/s",
    x = "Vs30 (m/s)",
    y = "P(SoilProfile = B)"
  ) +
  theme_minimal()

combined_mapping <- (p_material + p_psi0) / p_soil
ggsave(file.path(fig_dir, "building_mapping_heuristics.png"), combined_mapping,
       width = 12, height = 10, dpi = 150)
cat("Saved: building_mapping_heuristics.png\n")

# -----------------------------------------------------------------------------
# 8. Save mapping functions and sample data
# -----------------------------------------------------------------------------

# Save functions
save(
  sample_material_strength,
  prob_facade_A,
  assign_soil_profile,
  sample_initial_psi,
  generate_building_vulnerability,
  generate_vulnerability_batch,
  file = here::here("outputs", "models", "building_mapping_functions.RData")
)

# Save buildings with Vs30 joined
saveRDS(buildings, here::here("outputs", "models", "buildings_with_vs30.rds"))

cat("\nBuilding mapping functions and data saved.\n")
cat("\nBuilding mapping complete.\n")
