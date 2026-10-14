# 03_Modules

Can a single module eigengene tell one group from the other?

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Processing/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/module_classification.xlsx`, `b_reports/module_classification_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/03_Modules/a_script/module_classification.R` |
| Code | `00_Setup/Shared_analyses/run_feature_discrimination.R`, shared with the cluster groupings |

## Key result

- 15 modules tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
- No modules at p < 0.05 in any task.

## Method

- Five tasks: the groups compared at T1, T2 and T3, and on the training and acute changes.
- AUC is the Mann-Whitney W / (n1 n2), above 0.5 when the module is higher in the first group; Wilcoxon p, BH and Storey q within task.
- DeLong 95% interval for every AUC; stratified bootstrap (2,000 resamples) for the drawn modules (`run_feature_discrimination()` in `00_Setup/Shared_analyses/run_feature_discrimination.R`).

## Read with care

- The best AUC among many module eigengenes is high even when nothing separates the groups; read its interval.
- With seven or eight people per group the rank-sum p has a floor (`min_p` in `counts`).

## Sheets

`module_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per module and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
