# 02_Contrasts

Tests the modules on the nine contrasts, by eigengene and by fry, and reads each module's
membership against protein-level significance.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `02_Differential_Expression/02_Contrasts/c_data/fit.rds`, `02_Differential_Expression/04_Associate/c_data/associate.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 04_Network/02_Contrasts/a_script/02_contrasts.R` |
| Cost | about 5 s |

Membership against phenotype reads the protein-level training change scores, one column per
outcome. The figure PDF holds eigengene and fry tiles across the nine contrasts, membership against
each contrast and each phenotype outcome, and kME against moderated t for `Training_Interaction` in every
module.
