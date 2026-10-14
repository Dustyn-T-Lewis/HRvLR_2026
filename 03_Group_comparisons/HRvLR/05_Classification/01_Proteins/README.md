# 01_Proteins

Can a single protein tell one group from the other?

| | |
|---|---|
| Reads | `03_Group_comparisons/HRvLR/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds` |
| Writes | `c_data/protein_classification.xlsx`, `b_reports/protein_classification_report.pdf` |
| Run | `Rscript 03_Group_comparisons/HRvLR/05_Classification/01_Proteins/a_script/protein_classification.R` |
| Code | `00_Setup/Shared_analyses/run_feature_discrimination.R`, shared with the cluster groupings |

## Key result

- 1,900 proteins tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
- Proteins at p < 0.05 run from 44 (Baseline_HRvLR) to 100 (Acute_HRvLR).

## Method

- Five tasks: the groups compared at T1, T2 and T3, and on the training and acute changes.
- AUC is the Mann-Whitney W / (n1 n2), above 0.5 when the protein is higher in the first group; Wilcoxon p, BH and Storey q within task.
- DeLong 95% interval for every AUC; stratified bootstrap (2,000 resamples) for the drawn proteins (`run_feature_discrimination()` in `00_Setup/Shared_analyses/run_feature_discrimination.R`).

## Read with care

- The best AUC among many proteins is high even when nothing separates the groups; read its interval.
- With seven or eight people per group the rank-sum p has a floor (`min_p` in `counts`).

## Sheets

`protein_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per protein and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
