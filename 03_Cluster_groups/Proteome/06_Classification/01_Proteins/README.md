# 01_Proteins

Can a single protein tell one group from the other?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md)), run through the same function in `R/classify.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds` |
| Writes | `<input>/c_data/protein_classification.xlsx`, `<input>/b_reports/protein_classification_report.pdf`, `c_data/protein_classification_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/06_Classification/01_Proteins/a_script/protein_classification.R` |
| Code | `R/classify.R` |

## Key result

- `T1`
  - 1,900 proteins tested in each of 5 tasks; 64 at FDR < 0.05 in any task.
  - Proteins at p < 0.05 run from 79 (Acute_change_C1vC2) to 370 (Baseline_C1vC2).
- `T2`
  - 1,900 proteins tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
  - Proteins at p < 0.05 run from 41 (Acute_C1vC2) to 225 (Acute_change_C1vC2).
- `T3`
  - 1,900 proteins tested in each of 5 tasks; 63 at FDR < 0.05 in any task.
  - Proteins at p < 0.05 run from 45 (Training_change_C1vC2) to 424 (Acute_C1vC2).
- `training`
  - 1,900 proteins tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
  - Proteins at p < 0.05 run from 38 (Acute_C1vC2) to 305 (Training_change_C1vC2).
- `acute`
  - 1,900 proteins tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
  - Proteins at p < 0.05 run from 37 (Baseline_C1vC2) to 322 (Acute_change_C1vC2).

## Sheets

`protein_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per protein and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
