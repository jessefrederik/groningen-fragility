# Groningen Building Fragility Model

Statistical model for predicting earthquake-induced building damage in the Groningen gas field region, based on Korswagen's FEM simulation data.

## Overview

This project implements a hierarchical hurdle-Gamma model to estimate fragility curves P(Damage | PGV) for masonry buildings. The model:

- Uses Bayesian inference (brms/Stan) with full posterior uncertainty
- Enforces **monotonicity** (higher PGV → equal or more damage)
- Accounts for record-to-record ground motion variability
- Maps building covariates to FEM vulnerability parameters

## Key Results

| Metric | Value |
|--------|-------|
| 10% exceedance threshold (undamaged) | 15.3 mm/s |
| Lognormal median | 72.9 mm/s |
| Lognormal β | 1.11 |
| Gamma shape | 1.68 |
| Multi-event scaling | Sub-linear (16 events ≈ 4× damage) |

## Project Structure

```
├── datafiles/           # Input data (FEM, buildings, Vs30)
├── helperfuncties/      # Ground motion model, KNMI API
├── scripts/             # Analysis pipeline
│   ├── 01_prepare_data.R
│   ├── 02_explore_data.R
│   ├── 03_fit_baseline.R
│   ├── 04_fit_brms.R        # Main Bayesian model
│   ├── 05_validate.R
│   ├── 06_accumulation.R
│   ├── 07_building_mapping.R
│   ├── 08_spatial_simulation_brms.R
│   ├── 09_fragility_curves_brms.R
│   ├── 10_multievent_fast.R    # Fast cumulative damage simulation
│   └── helpers/
│       └── brms_predict.R
├── outputs/
│   ├── models/          # Fitted models (.rds)
│   └── figures/         # Diagnostic plots
└── METHODOLOGY_REPORT.md
```

## Dependencies

```r
install.packages(c("tidyverse", "brms", "mgcv", "sf", "bayesplot",
                   "patchwork", "truncnorm", "here"))
```

## Usage

Run scripts in order (01-09). Script 04 takes ~50 minutes for MCMC sampling.

```r
source("scripts/01_prepare_data.R")
# ... etc
```

## Data Sources

- **FEM data**: Korswagen (2019) PhD Thesis, TU Delft
- **Ground motion model**: Bommer et al. (2022)
- **Building data**: BAG (Basisregistratie Adressen en Gebouwen)

## References

- Korswagen, P. A. (2019)."; PhD Thesis, TU Delft
- Bommer, J. J., et al. (2022). Ground motion model for Groningen

## License

MIT

---
*Generated with Claude Code*
