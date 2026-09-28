# 03_Imputation

Fills the normalised matrix with missForest for the methods that need it complete.

| | |
|---|---|
| Reads | `02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/DAList_imputed.rds`, `c_data/03_impute.xlsx` |
| Run | `Rscript 01_Preprocess/03_Imputation/a_script/03_impute.R` |

## Method

missForest with the ranger backend, 10 iterations, 100 trees and seed 42, rows sorted before the
fit. Only fry, singscore and WGCNA read the imputed matrix; every limma fit reads the normalised
one.

## Figures

None.

## Workbook

| Sheet | Holds |
|---|---|
| `imputation` | missForest settings, share of values imputed, out-of-bag error. |
