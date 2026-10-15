# S1

Technical supplement to F1: whether the biopsies cluster, where the variance sits, whether one
group varies more, and what drives the main principal components.

| Panel | Shows | Reads |
|---|---|---|
| a_no_clusters | Biopsy-by-biopsy correlation in the saved average-linkage order, with HR/LR, time point and person strips; beside it the gap statistic for k = 1 to 6 with its standard error and the k the one-standard-error rule picks | `biopsy_correlation`, `biopsy_order`, `gap`, `pc_scores` |
| b_variance_shares | Each protein's variance shares for person, time point, HR/LR, blood index and residual (variancePartition), with medians | `variance_shares` |
| c_diffvar | DiffVar log variance ratio (LR relative to HR) per time point, shaded where p < 0.05; counts per direction against chance (2.5% of tested each way), FDR count, and ratios beyond ±4, which are counted but not drawn. Exploratory: 8 people per group | `diffvar` |
| d_pc_covariates | R² of each covariate alone on PCs 1 to 5: time point, HR/LR, person, blood index, proteins detected, sample weight | `pc_covariates` |

All sheets come from `02_HRvLR/04_Variance/c_data/variance.xlsx`.

`b_reports/S1_pc_scatters.pdf` has every PC-covariate pair behind panel d, one covariate per page;
it is written by `d_pc_covariates.R` and is not part of the composite.

## Run

```sh
Rscript 04_Figures/S1/a_script/S1_composite.R
Rscript 04_Figures/S1/a_script/panels/c_diffvar.R      # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel on its own |
| `b_reports/S1_composite.pdf` | the figure, 180 × 167 mm |
| `b_reports/S1_pc_scatters.pdf` | every PC against every covariate, six pages |
| `c_data/S1.xlsx` | one sheet per panel with the numbers drawn, column meanings, inputs with md5 |

## Read with care

- DiffVar tests only proteins observed in all 16 biopsies of a time point; a protein one group never
  detected has modelled values with no spread, and the test would read detection.
