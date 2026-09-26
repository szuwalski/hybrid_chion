# Software dependencies

The accepted fit was produced on Windows with R 4.4.2. Run scripts from the repository root.

Core model and plotting packages:

- RTMB 1.7
- ggplot2 4.0.3
- patchwork 1.3.0
- ggridges 0.5.6
- GGally
- scales
- rlang

GAM and data-processing packages:

- mgcv
- dplyr
- tidyr

R-recommended packages used by the workflow include `MASS`, `splines`, and `tools`.

The analysis-ready model inputs, component weight-at-size table, annual centroid data, and crabpack caches are bundled in the repository. Reproducing the accepted model and GAMs therefore does not require a new crabpack API query.

For the manuscript, the local workflow additionally requires `rmarkdown`, Pandoc, and a LaTeX installation; see `manuscript/README.md`.
