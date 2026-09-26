# Data

This directory contains protected source and cached inputs used by the model and its upstream analyses. Do not delete, replace, or exclude user-supplied data during routine cleanup.

Key inputs include:

- `growth_fit_observations.csv` — parental molt-increment observations used to fix parental growth.
- `sea_ice_yearly.csv` — annual sea-ice covariate used by the GAM analyses.
- `cache/crabpack_parental_raw.rds` and `cache/crabpack_hybrid_raw.rds` — complete station-level crabpack extracts used for spatial covariates, abundance centroids, empirical weight-at-size, and future extensions.
- `bsfrf/` — original Snow and Tanner BSFRF survey-efficiency material.
- `catch/` — fishery catch inputs.

The crabpack caches are expensive to recreate because the API is slow. The canonical model and GAM workflows use the saved inputs and do not require an API pull.

Analysis-ready model inputs derived from these sources are stored with the canonical model under `models/model/data/` so that the accepted fit can be inspected and its objective reconstructed.
