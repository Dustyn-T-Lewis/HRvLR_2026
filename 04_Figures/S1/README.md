# S1

Technical supplement to F1: whether the biopsies cluster, where the variance sits, whether one
group varies more, and what drives the main principal components.

| Panel | Shows | Reads |
|---|---|---|
| a_no_clusters | Biopsy-by-biopsy correlation in the saved average-linkage order, with HR/LR, time point and person strips; beside it the gap statistic for k = 1 to 6 with its standard error and the k the one-standard-error rule picks | `biopsy_correlation`, `biopsy_order`, `gap`, `pc_scores` |
| b_variance_shares | Each protein's variance shares for person, time point, HR/LR, blood index and residual (variancePartition), with medians | `variance_shares` |
| c_diffvar | DiffVar log variance ratio (LR relative to HR) per time point, shaded where p < 0.05; counts per direction against chance (2.5% of tested each way), FDR count, and ratios beyond ±4, which are counted but not drawn. Exploratory: 8 people per group | `diffvar` |
| d_pc_covariates | R² of each covariate alone on PCs 1 to 5: time point, HR/LR, person, blood index, proteins detected, sample weight | `pc_covariates` |
| e_pca_levels | PCA of the T2 and T3 levels, one point per person, drawn as Figure 1b (`style/pca_grid.R`) with each time point's PERMANOVA, gap k and forced-split ARI | `pca`, `permanova` from `01_Protein_Contrasts`; `separation` |

Sheets come from `02_HRvLR/04_Variance/c_data/variance.xlsx` unless named.

`b_reports/S1_pc_scatters.pdf` has every PC-covariate pair behind panel d, one covariate per page;
`a_script/extra/pc_scatters.R` writes it when the composite runs. It is not a panel.

## Run

```sh
Rscript 04_Figures/S1/a_script/S1_composite.R
Rscript 04_Figures/S1/a_script/panels/c_diffvar.R      # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel on its own |
| `b_reports/S1_composite.pdf` | the figure, 250 × 167 mm landscape |
| `b_reports/S1_pc_scatters.pdf` | every PC against every covariate, six pages |
| `c_data/S1.xlsx` | one sheet per panel with the numbers drawn, column meanings, inputs with md5 |

## Read with care

- DiffVar tests only proteins observed in all 16 biopsies of a time point; a protein one group never
  detected has modelled values with no spread, and the test would read detection.
- HR people vary more than LR people at T1 and T3 (panel c). At p < 0.05, 50 proteins are more
  variable in HR against 11 in LR at T1, and 87 against 5 at T3; chance gives 27 and 17 each way.
  T2 is even (31 against 36, chance 22). No protein passes FDR < 0.05, and each group has 8
  people.
