# 06_Associate

What the HR/LR grouping adds to the phenotype links, for proteins, set scores and module
eigengenes. The pooled, label-free associations are in `04_Continuous_no_group/01_Associate`.

| | |
|---|---|
| Reads | `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `01_Trunk/05_Scores/c_data/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/06_Modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/06_associate.xlsx`, `b_reports/06_associate_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Associate/a_script/06_associate.R` |

## Method

Pairings as in the continuous arc, including the T3 level against the T2 trait levels.

- Within group: Spearman between each feature and each outcome in HR and in LR, exact p from
  `cor.test`.
- Group difference: limma on `feature ~ group * outcome`, HR's slope minus LR's.
- Flips: DGCA's Fisher z difference between the HR and LR Spearman correlations, with DGCA's own
  p for that difference (`pValDiff`). DGCA's permutation option was dropped: its q-value step fails
  on the 15-module screen.

`comp_hypertrophy` stays out of the group difference and flips, since the groups are its split.
Each group needs five subjects with the outcome measured. BH within level, window, outcome and
collection.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per test, level, window, outcome and group. |
| `within_group` | Spearman rho and exact p in each group, every feature and outcome. |
| `group_difference` | HR minus LR slope, t, p, BH, Storey q. |
| `flips` | Each group's rho, Fisher z difference, DGCA's p, BH, Storey q, class. |
