# Standard joint population model

This directory is the single current model implementation. It corresponds to the final `bridge_14` fit selected during model development. Earlier implementations, bridge fits, and diagnostics are retained under `../archive/` and are historical rather than runnable.

## Canonical files

- `fits/model_fit.rds`: accepted fitted model object.
- `fits/hessian_diagnostics.csv`: numerical Hessian diagnostics.
- `plots/model_diagnostics/`: complete standard diagnostic figures and tables.
- `data/model_input.rds`: configured model input.
- `data/two_fishery_base_input.rds`: base two-fishery input.
- `run_model.R`: stable fitting entry point; its defaults reproduce the accepted configuration.
- `check_hessian.R` and `plot_model.R`: stable diagnostic entry points.
- `projection/`: projection scripts and local helper functions.

## Accepted assumptions

- Snow, Tanner, and hybrid crab are modeled jointly by sex, with two fisheries.
- Survey selectivity uses the legacy two-fishery formulation retained from the model that best reproduced the survey patterns. BSFRF-derived information remains available upstream in `../../analyses/bsfrf_selectivity/`.
- Directed fishery selectivity is ascending logistic; incidental selectivity may be dome-shaped.
- There is no Snow-fishery mortality in 2022 and the anomalous retained Snow observation is excluded.
- Natural mortality is parameterized as a hazard and reported as annual probability on a 0-1 scale.
- Annual mortality-deviation SD is 1 except immature Snow females, for which it is 0.6.
- Baseline mortality prior SD is 0.10 for mature Snow males, 0.10 for mature Snow females, and 0.05 for immature Snow females.
- Early-era female survey curves are partially pooled toward the current era with SD 3 mm for midpoint offsets and 0.10 for log-slope offsets.
- Parental growth is fixed from `../../data/growth_fit_observations.csv`; hybrid growth is an estimated parental blend.
- Hybrid misclassification follows a power curve below 90 mm and is zero at and above 90 mm.
- Discard survival is 70%.

The accepted fit has objective 56018.7427, maximum absolute gradient 0.000793, and a positive-definite invertible Hessian. See `fits/hessian_diagnostics.csv` for the full checks.

## Run

From the repository root:

```powershell
Rscript models/model/run_model.R
Rscript models/model/check_hessian.R
Rscript models/model/plot_model.R
```

Projection scripts are deliberately kept with the model:

```powershell
Rscript models/model/projection/project_model.R --years=15 --draws=100 --seed=20260901
Rscript models/model/projection/plot_projection_recruitment_gam.R
```

Time-series figures break lines across missing years, fishery size-composition observations are bars with predicted lines, and mortality figures use a 0-1 scale.
