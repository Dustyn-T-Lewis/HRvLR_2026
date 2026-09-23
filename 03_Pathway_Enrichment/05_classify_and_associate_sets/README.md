# 05_classify_and_associate_sets

How well each set's score separates study groups, and whether it tracks the phenotype.

| | |
|---|---|
| **Reads** | `gene_sets.rds`, `set_tests.rds`, `singscore.rds`, `DAList_imputed.rds`, `phenotype.csv` |
| **Writes** | `set_results.rds`, `05_classify_and_associate_sets.xlsx`, figures |

Eight tasks. Training (T1 to T2) and the acute bout (T2 to T3) within each arm are paired, with the
p from the Wilcoxon signed-rank test. HR against LR at T1, at T2, in training change and in acute
change are unpaired, with the rank-sum test. AUC comes from `pROC` with `direction = "<"` fixed.
Associations are Spearman (t approximation) in three windows: training change, T1 level and acute
change, with subjects pooled and again within each arm.

Each collection is read against its own chance count (5% of its sets), with BH within collection
and task. The tables are in `../README.md`.

Figures: ROC curves for the two strongest sets per collection per task (floor excluded), the same
for associations per window, nominal hits over chance by collection, and hit matrices of every
nominal set per task and per window.

```r
sr <- readRDS("03_Pathway_Enrichment/05_classify_and_associate_sets/c_data/set_results.rds")
sr$chance_expectation # nominal against chance, per collection
sr$set_auc            # set x task: AUC, p, BH
sr$set_association    # set x outcome x window, pooled
sr$by_arm             # correlations within each arm
```
