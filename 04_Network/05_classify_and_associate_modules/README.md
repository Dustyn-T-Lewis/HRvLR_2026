# 05_classify_and_associate_modules

The eigengenes on the eight classification tasks and the ten phenotypes.

| | |
|---|---|
| Reads | `modules.rds`, `00_Input/phenotype.csv` |
| Writes | `module_results.rds`, `05_classify_and_associate_modules.xlsx`, 5 figures over 7 pages |

The eight tasks are those of `02_Differential` and `03_Pathway_Enrichment/05`: pROC AUC with
`direction = "<"`, Wilcoxon p, paired within arm. Associations are Spearman, t approximation, in
the training, baseline and acute windows, pooled and within each arm. Eigengenes have no missing
values, so no observation floor applies.

| Analysis | Tests | Expected at p < 0.05 | Nominal | BH < 0.05 |
|---|---:|---:|---:|---:|
| classification | 96 | 4.8 | 6 | 0 |
| association | 360 | 18 | 12 | 0 |

BH over twelve modules is lenient; no module passes it. Figures: every module and task as tiles,
every module and phenotype per window as tiles, an ROC panel for each nominal module-task pair
(floor excluded), a scatter for each nominal module-outcome pair, and nominal hits over chance.
