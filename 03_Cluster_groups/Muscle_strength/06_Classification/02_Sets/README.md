# 02_Sets

Can a single gene-set score tell one group from the other?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`02_HRvLR_groups/06_Classification/02_Sets`](../../../../02_HRvLR_groups/06_Classification/02_Sets/README.md)), run through the same function in `00_Setup/Shared_analyses/run_feature_discrimination.R`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/set_classification.xlsx`, `b_reports/set_classification_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classification/02_Sets/a_script/set_classification.R` |
| Code | `00_Setup/Shared_analyses/run_feature_discrimination.R` |

## Key result

- 1,378 gene sets tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
- Gene Sets at p < 0.05 run from 6 (Acute_change_HivLo) to 214 (Baseline_HivLo).

## Sheets

`set_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per set and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
