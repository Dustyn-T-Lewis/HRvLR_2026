# 03_Imputation

Fills the normalised matrix with missForest for the methods that need it complete.

| | |
|---|---|
| Reads | `02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/DAList_imputed.rds`, `c_data/03_impute.xlsx` |
| Run | `quarto render 01_Preprocess/03_Imputation/a_script/03_impute.qmd --output-dir ../b_reports` |
| Cost | about 30 s |

Read by `02_Differential` (the `fry` appendix), `03_Pathway_Enrichment/01`, `/04` and `/05` (the
sample sheet only), and `04_Network/01`, `/03` and `/04`.
