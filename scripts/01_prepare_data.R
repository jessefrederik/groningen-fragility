# =============================================================================
# 01_prepare_data.R
#
# Load and transform Korswagen FEM data for fragility modeling.
# Creates analysis-ready dataset with derived variables.
#
# @author Generated with Claude Code
# @date 2026-01-03
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1. Load raw data
# -----------------------------------------------------------------------------

fem_raw <- read_csv(
  here::here("datafiles", "korswagen_fem_summary.csv"),
  col_types = cols(
    PGV = col_integer(),
    N = col_integer(),
    DesignInitialPsi = col_double(),
    InitialPsi = col_double(),
    FinalPsi = col_double(),
    DeltaPsi = col_double(),
    Material = col_double(),
    FacadeType = col_character(),
    SoilProfile = col_character(),
    EarthquakeType = col_character()
  )
)

cat("Loaded", nrow(fem_raw), "rows from korswagen_fem_summary.csv\n")

# -----------------------------------------------------------------------------
# 2. Create derived variables
# -----------------------------------------------------------------------------

fem <- fem_raw |>
  mutate(
    # Log-transform predictors
    log_pgv = log(PGV),
    log_n = log(N),

    # Clip negative DeltaPsi (numerical noise) to 0
    delta_psi_pos = pmax(DeltaPsi, 0),

    # Binary indicator for damage occurrence
    has_damage = as.integer(delta_psi_pos > 0),

    # Create ordered factor for PGV (for monotonic effects in brms)
    pgv_cat = factor(PGV, levels = c(2, 4, 8, 16, 32, 64, 96, 128), ordered = TRUE),

    # Standardize continuous predictors for model stability
    log_pgv_z = scale(log_pgv)[, 1],
    log_n_z = scale(log_n)[, 1],
    material_z = scale(Material)[, 1],
    initial_psi_z = scale(InitialPsi)[, 1],

    # Create combined factor for facade x soil interaction
    facade_soil = interaction(FacadeType, SoilProfile, drop = TRUE),

    # Earthquake record as factor
    eq_record = factor(EarthquakeType)
  )

# -----------------------------------------------------------------------------
# 3. Summary statistics
# -----------------------------------------------------------------------------

cat("\n=== Variable Summary ===\n")
cat("PGV levels:", sort(unique(fem$PGV)), "\n")
cat("N levels:", sort(unique(fem$N)), "\n")
cat("Material levels:", sort(unique(fem$Material)), "\n")
cat("FacadeType:", unique(fem$FacadeType), "\n")
cat("SoilProfile:", unique(fem$SoilProfile), "\n")
cat("EarthquakeType:", unique(fem$EarthquakeType), "\n")
cat("DesignInitialPsi:", sort(unique(fem$DesignInitialPsi)), "\n")

cat("\n=== DeltaPsi Distribution ===\n")
cat("Total rows:", nrow(fem), "\n")
cat("DeltaPsi == 0:", sum(fem$delta_psi_pos == 0),
    sprintf("(%.1f%%)", 100 * mean(fem$delta_psi_pos == 0)), "\n")
cat("DeltaPsi > 0:", sum(fem$delta_psi_pos > 0), "\n")
cat("Original negative:", sum(fem$DeltaPsi < 0), "(clipped to 0)\n")
cat("Mean DeltaPsi (all):", round(mean(fem$delta_psi_pos), 3), "\n")
cat("Mean DeltaPsi (>0):", round(mean(fem$delta_psi_pos[fem$has_damage == 1]), 3), "\n")
cat("Max DeltaPsi:", round(max(fem$delta_psi_pos), 3), "\n")

# -----------------------------------------------------------------------------
# 4. Create train/test splits for cross-validation
# -----------------------------------------------------------------------------

# Leave-one-record-out CV: 4 folds, one per EarthquakeType
cv_folds <- fem |>
  mutate(fold = as.integer(factor(EarthquakeType)))

cat("\n=== Cross-Validation Folds (Leave-One-Record-Out) ===\n")
cv_folds |>
  count(fold, EarthquakeType) |>
  print()

# -----------------------------------------------------------------------------
# 5. Check for missing PGV combinations (truncation)
# -----------------------------------------------------------------------------

cat("\n=== PGV x N combinations (checking truncation) ===\n")
fem |>
  filter(DesignInitialPsi == 0) |>  # Focus on undamaged starting condition
  count(PGV, N) |>
  pivot_wider(names_from = N, values_from = n, names_prefix = "N=") |>
  print(n = 10)

# Note: Low PGV has fewer observations at some N levels due to computational
# triage in original study ("no damage expected" → not run)

# -----------------------------------------------------------------------------
# 6. Save prepared data
# -----------------------------------------------------------------------------

saveRDS(fem, here::here("outputs", "models", "fem_prepared.rds"))
cat("\nSaved prepared data to outputs/models/fem_prepared.rds\n")

# Also save standardization parameters for prediction
std_params <- list(
  log_pgv_mean = mean(fem_raw$PGV |> log()),
  log_pgv_sd = sd(fem_raw$PGV |> log()),
  log_n_mean = mean(fem_raw$N |> log()),
  log_n_sd = sd(fem_raw$N |> log()),
  material_mean = mean(fem_raw$Material),
  material_sd = sd(fem_raw$Material),
  initial_psi_mean = mean(fem_raw$InitialPsi),
  initial_psi_sd = sd(fem_raw$InitialPsi)
)
saveRDS(std_params, here::here("outputs", "models", "standardization_params.rds"))
cat("Saved standardization parameters to outputs/models/standardization_params.rds\n")

# -----------------------------------------------------------------------------
# 7. Quick data quality checks
# -----------------------------------------------------------------------------

cat("\n=== Data Quality Checks ===\n")

# Check for any NA values
na_counts <- colSums(is.na(fem))
if (any(na_counts > 0)) {
  cat("WARNING: Missing values found:\n")
  print(na_counts[na_counts > 0])
} else {
  cat("No missing values found.\n")
}

# Check FinalPsi >= InitialPsi (should always be true if damage accumulates)
violations <- sum(fem$FinalPsi < fem$InitialPsi - 0.001)  # small tolerance
if (violations > 0) {
  cat("WARNING:", violations, "rows where FinalPsi < InitialPsi\n")
} else {
  cat("All FinalPsi >= InitialPsi (damage accumulates correctly).\n")
}

cat("\nData preparation complete.\n")
