# 03_preserve_modules

Builds modules inside each arm and tests their preservation in the other arm.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `01_build_modules/c_data/modules.rds` |
| Writes | `c_data/03_preserve_modules.xlsx`, `b_reports/03_preserve_modules_figures.pdf` |
| Run | `Rscript 04_Network/03_preserve_modules/a_script/03_preserve_modules.R` |
| Cost | about 3.5 min, most of it the 200 permutations |

The `arm_networks` sheet holds each arm's samples, soft power and module count. Each arm module is
labelled with its best-overlapping full-cohort module.
