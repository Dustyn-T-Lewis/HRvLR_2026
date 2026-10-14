# 03_Modules

Can a single module eigengene tell one group from the other?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`03_Group_comparisons/HRvLR/05_Classification/03_Modules`](../../../../03_Group_comparisons/HRvLR/05_Classification/03_Modules/README.md)), run through the same function in `00_Setup/Shared_analyses/run_feature_discrimination.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Group_comparisons/Proteome/01_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/06_Module_Building/c_data/modules.rds` |
| Writes | `<input>/c_data/module_classification.xlsx`, `<input>/b_reports/module_classification_report.pdf`, `c_data/module_classification_all_inputs.xlsx` |
| Run | `Rscript 03_Group_comparisons/Proteome/05_Classification/03_Modules/a_script/module_classification.R` |
| Code | `00_Setup/Shared_analyses/run_feature_discrimination.R` |

## Key result

- `T1`
  - 15 modules tested in each of 5 tasks; 2 at FDR < 0.05 in any task.
  - Modules at p < 0.05 run from 0 (Acute_change_C1vC2) to 6 (Baseline_C1vC2).
- `T2`
  - 15 modules tested in each of 5 tasks; 6 at FDR < 0.05 in any task.
  - Modules at p < 0.05 run from 0 (Acute_C1vC2) to 4 (Trained_C1vC2).
- `T3`
  - 15 modules tested in each of 5 tasks; 9 at FDR < 0.05 in any task.
  - Modules at p < 0.05 run from 0 (Baseline_C1vC2) to 9 (Acute_C1vC2).
- `training`
  - 15 modules tested in each of 5 tasks; 3 at FDR < 0.05 in any task.
  - Modules at p < 0.05 run from 0 (Acute_C1vC2) to 6 (Training_change_C1vC2).
- `acute`
  - 15 modules tested in each of 5 tasks; 8 at FDR < 0.05 in any task.
  - Modules at p < 0.05 run from 0 (Baseline_C1vC2) to 6 (Acute_change_C1vC2).

## Sheets

`module_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per module and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
