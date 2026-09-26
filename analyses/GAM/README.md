# Current GAM analyses

This is the single current GAM output folder. Scripts are in `scripts/`; their defaults read `models/model/fits/model_fit.rds` when run from the repository root and write directly here.

Output prefixes identify the analysis: `mortality_`, `male_terminal_molt_`, `female_terminal_molt_`, `recruitment_centroid_`, and `hybrid_recruitment_`.

Run from the repository root:

```powershell
Rscript analyses/GAM/scripts/mortality_gam.R
Rscript analyses/GAM/scripts/terminal_molt_gam.R
Rscript analyses/GAM/scripts/female_terminal_molt_gam.R
Rscript analyses/GAM/scripts/recruitment_centroid_gam.R
Rscript analyses/GAM/scripts/male_abundance_surface_gam.R
```

`data/annual_component_centroids.csv` is the retained input needed by the recruitment-centroid GAM.
