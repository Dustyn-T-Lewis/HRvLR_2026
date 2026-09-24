# 05_classify_and_associate_modules

Puts the eigengenes through the eight classification tasks and the ten phenotypes.

| | |
|---|---|
| Reads | `01_build_modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/05_classify_and_associate_modules.xlsx`, `b_reports/05_classify_and_associate_modules_figures.pdf` |
| Run | `Rscript 04_Network/05_classify_and_associate_modules/a_script/05_classify_and_associate_modules.R` |
| Cost | about 6 s |

The tasks, AUC and tests are those of `02_Differential` and `03_Pathway_Enrichment/05`. Eigengenes
have no missing values, so no observation floor applies. The figure PDF holds nominal hits over
chance, every module and task as tiles, every module and phenotype per window as tiles, an ROC panel
for each nominal module-task pair (floor excluded) and a scatter for each nominal module-outcome
pair. `module_by_arm` holds every within-arm correlation.
