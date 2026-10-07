# 02_Sets

Can a single gene-set score tell one group from the other?

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/set_classification.xlsx`, `b_reports/set_classification_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/02_Sets/a_script/set_classification.R` |
| Code | `R/classify.R`, shared with the cluster groupings |

## Key result

- 1,378 gene sets tested in each of 5 tasks; 0 at FDR < 0.05 in any task.
- Gene Sets at p < 0.05 run from 12 (Trained_HRvLR) to 52 (Acute_HRvLR).

## Method

- Five tasks: the groups compared at T1, T2 and T3, and on the training and acute changes.
- AUC is the Mann-Whitney W / (n1 n2), above 0.5 when the set is higher in the first group; Wilcoxon p, BH and Storey q within task, and within collection.
- DeLong 95% interval for every AUC; stratified bootstrap (2,000 resamples) for the drawn sets (`classify_step()` in `R/classify.R`).

## Read with care

- The best AUC among many gene-set scores is high even when nothing separates the groups; read its interval.
- With seven or eight people per group the rank-sum p has a floor (`min_p` in `counts`).

## Sheets

`set_classification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per set and task: AUC with DeLong 95% interval, Wilcoxon p, BH and Storey q; bootstrap 95% interval for the features drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every drawn row. |
