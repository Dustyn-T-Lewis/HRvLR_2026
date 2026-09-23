# 03_classify_and_associate_modules

The eigengenes put through the nine contrasts, the eight classification tasks and the ten
phenotypes.

| | |
|---|---|
| **Reads** | `modules.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `00_Input/phenotype.csv` |
| **Writes** | `module_results.rds`, `03_classify_and_associate_modules.xlsx`, 5 figures |

**Contrasts.** `lmFit()` on the eigengenes with the protein design and the subject block,
correlation re-estimated on the eigengenes (0.114), `eBayes(robust = TRUE)`, BH within contrast.

**Classification.** The same eight tasks as `02_Differential` and `03_Pathway_Enrichment/05`:
pROC AUC with `direction = "<"`, Wilcoxon p, paired within arm.

**Association.** Spearman, t approximation, in the training, baseline and acute windows.

| Analysis | Tests | Expected at p < 0.05 | Nominal | BH < 0.05 |
|---|---:|---:|---:|---:|
| contrasts | 108 | 5.4 | 5 | 0 |
| classification | 96 | 4.8 | 6 | 0 |
| association | 360 | 18 | 10 | 0 |

With twelve features BH is a weak filter, and still nothing passes it. Every module is drawn in
each figure, one tile per module and column.
