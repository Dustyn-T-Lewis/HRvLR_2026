# 07_Trait_Associations

What the HR/LR grouping adds to the phenotype links, for proteins, set scores and module
eigengenes. The pooled, label-free associations are in `04_Continuous_no_group/01_Trait_Associations`.

| | |
|---|---|
| Reads | `01_Trunk/03_Quantification/c_data/DAList_corrected.rds`, `01_Trunk/05_Set_Scores/c_data/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/06_Module_Building/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/06_associate.xlsx`, `b_reports/06_associate_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/07_Trait_Associations/a_script/06_associate.R` |

## Method

Pairings as in the continuous arc, including the T3 level against the T2 trait levels.

- Within group: Spearman between each feature and each trait in HR and in LR, exact p from
  `cor.test`. Descriptive: each group holds six to eight people, and a group split on a trait
  restricts that trait's range inside it.
- Group difference: limma on `feature ~ group * trait`, HR's slope minus LR's. This is the one test
  of whether a link differs between the groups; a link significant in one group and not the other
  is not, on that basis, a difference.

`comp_hypertrophy` stays out of the group difference, since the groups are its split. Each group
needs five people with the trait measured. BH within level, window, trait, group and collection.

## Figures

One page: the share of features at p < 0.05 per test, window and trait, against the 5% expected by
chance.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per test, level, window, trait and group. |
| `within_group` | Spearman rho and exact p in each group, every feature and trait. |
| `group_difference` | HR minus LR slope, t, p, BH and Storey q. |
