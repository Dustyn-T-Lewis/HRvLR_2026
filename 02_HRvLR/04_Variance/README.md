# 04_Variance

Describes where the proteome's variation lies: variance shares per protein (variancePartition),
whether HR or LR people spread more (DiffVar, exploratory), whether biopsies form natural clusters
(gap statistic), and what drives the main principal components. Feeds supplement S1.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Quantification/c_data/proteins.rds` | limpa protein values and standard errors |
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | blood-corrected values, for clustering and PCA |
| `02_HRvLR/01_Protein_Contrasts/c_data/fit.rds` | sample weights |

## Run

```sh
quarto render 02_HRvLR/04_Variance/a_script/variance.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/variance.xlsx` | variance shares, DiffVar, gap table, separation per time point and per change (gap k, forced-split adjusted Rand index against HR/LR), biopsy correlation and order, PC scores with covariates, PC-by-covariate R² |
| `b_reports/variance.html` | the rendered walkthrough |

## Read with care

- The gap statistic picks k = 1 across all biopsies and at each time point: the biopsies form no
  natural clusters. Forced into two clusters, they match HR/LR no better than chance (adjusted Rand
  index −0.07 to 0.00); at T2 and T3 the second cluster is S28's thin biopsy alone.
- PC1 of the blood-corrected values follows the number of proteins each biopsy detected (R² 0.85):
  the thin S28 biopsies sit at its extreme. PCA, clustering and DiffVar include them unweighted.
- DiffVar runs only on proteins observed in all 16 biopsies of a time point (1,076 at T1, 870 at
  T2, 673 at T3): limpa's modelled values for a protein a group never detected carry no spread, and
  the test then reads detection. No protein reaches FDR < 0.05. HR is more variable at p < 0.05 for
  50 proteins at T1 and 87 at T3, against about 27 and 17 by chance; with 8 people per group these
  are leads, not findings.
