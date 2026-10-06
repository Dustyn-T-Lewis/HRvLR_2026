# 07_Associate

A copy of `02_HRvLR_groups/06_Associate` run on the muscle size clusters: groups Hi and Lo from `01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

What the Hi/Lo grouping adds to the phenotype links, for proteins, set scores and module
eigengenes. The pooled, label-free associations are in `04_Continuous_no_group/01_Associate`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/DAList_corrected.rds`, `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/06_associate.xlsx`, `b_reports/06_associate_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/07_Associate/a_script/07_associate.R` |

## Method

Pairings as in the continuous arc, including the T3 level against the T2 trait levels.

- Within group: Spearman between each feature and each outcome in Hi and in Lo, exact p from
  `cor.test`.
- Group difference: limma on `feature ~ group * outcome`, Hi's slope minus Lo's.
- Flips: DGCA's Fisher z difference between the Hi and Lo Spearman correlations, with DGCA's own
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
| `group_difference` | Hi minus Lo slope, t, p, BH, Storey q. |
| `flips` | Each group's rho, Fisher z difference, DGCA's p, BH, Storey q, class. |
