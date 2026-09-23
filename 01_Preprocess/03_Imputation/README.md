# 01_Preprocess / 03_Imputation

Fills the normalised matrix for the three methods that need it complete.

| | |
|---|---|
| **Script** | `a_script/03_impute.qmd` |
| **Reads** | `02_Normalization/c_data/DAList_normalized.rds` |
| **Writes** | `c_data/DAList_imputed.rds`, `c_data/03_impute.xlsx` |

```sh
quarto render 01_Preprocess/03_Imputation/a_script/03_impute.qmd --output-dir ../b_reports
```

missForest with the ranger backend named, 10 iterations, 100 trees, seed 42, rows sorted before
the fit. Out-of-bag NRMSE 0.151.

Read by `03_Pathway_Enrichment/01` (`fry`), `03_Pathway_Enrichment/04` (singscore) and
`04_Network/01` (WGCNA). The differential fit never sees an imputed value.
