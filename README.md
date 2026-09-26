# Joint Snow, Tanner, and hybrid crab population model

This repository contains the current size-structured rTMB assessment for Snow, Tanner, and hybrid crab. It jointly models male and female population components, two fisheries, survey observations, growth, maturation, recruitment, natural mortality, and hybrid misclassification.

The accepted model is the final bridge-14 configuration. Its saved fit has a positive-definite, invertible Hessian, objective 56018.7427, and maximum absolute gradient 0.000793.

## Repository map

- [`models/model/`](models/model/) — canonical model code, prepared input, accepted fit, diagnostics, and projection workflow.
- [`analyses/GAM/`](analyses/GAM/) — current process GAM scripts and all current GAM outputs.
- [`analyses/bsfrf_selectivity/`](analyses/bsfrf_selectivity/) — protected BSFRF/NMFS survey-efficiency analysis and bootstrap results.
- [`data/`](data/) — source and cached input data. These files must not be removed during routine cleanup.
- [`manuscript/`](manuscript/) — manuscript sources and rendered manuscript files.

Historical model and analysis directories remain available locally under `models/archive/` and `analyses/archive/`, but are intentionally excluded from Git because collaborators do not need the development fits.

## Reproduce current results

Use R 4.4.2 and install the packages listed in [`DEPENDENCIES.md`](DEPENDENCIES.md). Run commands from the repository root.

Missing R packages can be installed with `Rscript install_dependencies.R`.

To regenerate diagnostics and GAM results from the accepted saved fit:

```powershell
Rscript reproduce_results.R
```

Optional flags are `--refit`, `--project`, and `--draws=N`. Refitting is computationally expensive; the default workflow uses the accepted `models/model/fits/model_fit.rds`.

The equivalent individual entry points are:

```powershell
Rscript models/model/check_hessian.R
Rscript models/model/plot_model.R
Rscript analyses/GAM/scripts/mortality_gam.R
Rscript analyses/GAM/scripts/terminal_molt_gam.R
Rscript analyses/GAM/scripts/female_terminal_molt_gam.R
Rscript analyses/GAM/scripts/recruitment_centroid_gam.R
Rscript analyses/GAM/scripts/male_abundance_surface_gam.R
```

To refit or project:

```powershell
Rscript models/model/run_model.R
Rscript models/model/projection/project_model.R --years=15 --draws=100 --seed=20260901
```

See [`models/model/README.md`](models/model/README.md) for assumptions and model details, and [`CLEANUP_LOG.md`](CLEANUP_LOG.md) for repository history.

## Primary outputs

- Accepted fit: `models/model/fits/model_fit.rds`
- Hessian diagnostics: `models/model/fits/hessian_diagnostics.csv`
- Standard model figures: `models/model/plots/model_diagnostics/`
- Current GAM figures and tables: `analyses/GAM/`
- BSFRF selectivity results: `analyses/bsfrf_selectivity/`
