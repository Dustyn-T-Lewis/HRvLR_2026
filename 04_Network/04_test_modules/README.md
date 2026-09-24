# 04_test_modules

Tests the modules on the nine contrasts, by eigengene and by fry, and reads each module's
membership against protein-level significance.

| | |
|---|---|
| Reads | `01_build_modules/c_data/modules.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `02_Differential_Expression/02_Differential/c_data/fit.rds`, `02_Differential_Expression/03_Phenotype/c_data/phenotype.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/04_test_modules.xlsx`, `b_reports/04_test_modules_figures.pdf` |
| Run | `Rscript 04_Network/04_test_modules/a_script/04_test_modules.R` |
| Cost | about 5 s |

The figure PDF holds eigengene and fry tiles across the nine contrasts, membership against each
contrast and each phenotype, and kME against moderated t for `Training_Interaction` in every
module.
