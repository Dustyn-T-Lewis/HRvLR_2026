# 03_Imputation

Fills the normalised matrix with missForest for the methods that need it complete.

| | |
|---|---|
| Reads | `02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/DAList_imputed.rds`, `c_data/03_impute.xlsx` |
| Run | `Rscript 01_Preprocess/03_Imputation/a_script/03_impute.R` |
| Cost | about 30 s |

Read by every step that needs a complete matrix: `fry`, singscore and WGCNA.
