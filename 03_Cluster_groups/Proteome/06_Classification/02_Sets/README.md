# 02_Sets

Can a single gene-set score tell one group from the other?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`02_HRvLR_groups/06_Classification/02_Sets`](../../../../02_HRvLR_groups/06_Classification/02_Sets/README.md)), run through the same function in `R/classify.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `<input>/c_data/set_classification.xlsx`, `<input>/b_reports/set_classification_report.pdf`, `c_data/set_classification_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/06_Classification/02_Sets/a_script/set_classification.R` |
| Code | `R/classify.R` |

## Key result

- `T1`
  - 1,378 gene sets tested in each of 5 tasks; 139 at FDR < 0.05 in any task.
  - Gene Sets at p < 0.05 run from 60 (Acute_change_C1vC2) to 402 (Baseline_C1vC2).
- `T2`
  - 1,378 gene sets tested in each of 5 tasks; 9 at FDR < 0.05 in any task.
  - Gene Sets at p < 0.05 run from 14 (Acute_C1vC2) to 207 (Acute_change_C1vC2).
- `T3`
  - 1,378 gene sets tested in each of 5 tasks; 338 at FDR < 0.05 in any task.
  - Gene Sets at p < 0.05 run from 14 (Training_change_C1vC2) to 513 (Acute_C1vC2).
- `training`
  - 1,378 gene sets tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
  - Gene Sets at p < 0.05 run from 10 (Acute_C1vC2) to 134 (Training_change_C1vC2).
- `acute`
  - 1,378 gene sets tested in each of 5 tasks; 9 at FDR < 0.05 in any task.
  - Gene Sets at p < 0.05 run from 10 (Baseline_C1vC2) to 295 (Acute_change_C1vC2).

## Sheets

`set_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per set and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
