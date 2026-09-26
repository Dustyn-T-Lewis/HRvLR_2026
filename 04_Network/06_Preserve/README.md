# 06_Preserve

Builds modules inside each arm and tests their preservation in the other arm.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/06_preserve.xlsx`, `b_reports/06_preserve_figures.pdf` |
| Run | `Rscript 04_Network/06_Preserve/a_script/06_preserve.R` |
| Cost | about 3.5 min, most of it the 200 permutations |

The `arm_networks` sheet holds each arm's samples, soft power and module count. Each arm module is
labelled with its best-overlapping full-cohort module.
