# Cleanup record

## 2026-09-26 collaborator-share preparation

Historical model and analysis archives, scratch work, generated projection runs, and duplicate rendered exports were excluded from Git without deleting their local copies. The root documentation, data documentation, dependency list, canonical checksums, manuscript figure paths, and stable model entry points were updated. `reproduce_results.R` now regenerates the Hessian diagnostics, complete standard model plots, and all five current GAM analyses from the accepted fit in one command; optional flags enable refitting and projections.

The default reproduction workflow completed successfully after these changes.

---

## 2026-09-22 consolidation

The final bridge-14 two-fishery model was promoted to `models/model/`; all other model versions and intermediate bridge outputs were moved to `models/archive/`. The model's projection scripts and helper functions were retained inside `models/model/projection/` and verified with a 20-draw smoke projection.

The latest GAM scripts and all current, prefixed GAM outputs were consolidated in `analyses/GAM/`. The full BSFRF/NMFS survey-selectivity analysis was retained separately in `analyses/bsfrf_selectivity/`. Superseded analysis material was moved to `analyses/archive/`.

This consolidation used moves and copies only. No source data or historical analysis files were deleted. Archived scripts are preserved for provenance but are not required to remain runnable.

The canonical fit has an invertible positive-definite Hessian (objective 56018.7427; maximum absolute gradient 0.000793). All five canonical GAM scripts and the canonical projection workflow were run successfully after path updates.

---

Cleanup completed 2026-09-02. The repository was reduced to the current JMR
analysis and its accepted fit.

Removed material included standalone precursor models, J0/J1/JR/JM fits and
plots, exploratory recruitment/spatial/GAM products, legacy ADMB files, obsolete
10- and 15-year projection products, tests, session state, previews, and large
plot-derived tables that can be regenerated.

The JMR runner no longer reads JR or JM fits for starting values. JMR preparation
uses the compact `JMR_base_data.rds` input and explicitly records the recruitment
scaling constants used by the accepted fit. Rebuilding the configured data and
objective reproduced the saved objective exactly, and the Hessian remained
positive definite.

The complete crabpack parental and hybrid caches were subsequently restored for
spatial analyses. They are protected inputs and must not be deleted in future
cleanup unless the user explicitly asks for their removal.
