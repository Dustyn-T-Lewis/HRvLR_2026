# 01_Proteins

Can a single protein tell one group from the other?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md)), run through the same function in `R/classify.R`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds` |
| Writes | `c_data/protein_classification.xlsx`, `b_reports/protein_classification_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classification/01_Proteins/a_script/protein_classification.R` |
| Code | `R/classify.R` |

## Key result

- 1,900 proteins tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
- Proteins at p < 0.05 run from 26 (Acute_change_HivLo) to 123 (Training_change_HivLo).

## Sheets

`protein_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per protein and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
