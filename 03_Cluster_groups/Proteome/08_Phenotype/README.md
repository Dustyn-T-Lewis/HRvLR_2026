# 03_Phenotype

Do the proteome clusters differ in size or strength? Phenotype played no part in the clustering,
so these are tests.

| | |
|---|---|
| Reads | `01_Cluster/c_data/01_cluster.xlsx`, `00_Input/phenotype.csv` |
| Writes | `c_data/03_phenotype.xlsx`, `b_reports/03_phenotype_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/03_Phenotype/a_script/03_phenotype.R` |

## Method

Per input and trait (% change over training, composite, volume load, T1 and T2 levels): each
cluster's median and the Wilcoxon rank-sum test (`wilcox.test`), BH within input. The HR/LR
cross-table per input carries Fisher's exact p as a description.

## Workbook

| Sheet | Holds |
|---|---|
| `tests` | Per input, trait kind and trait: n and median per cluster, W, p, BH. |
| `hr_lr` | Per input: HR and LR counts per cluster, Fisher's p. |
