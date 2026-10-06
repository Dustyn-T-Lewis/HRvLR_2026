# 04_Interactions

Two interactions with the proteome clusters, at protein and module level.

| | |
|---|---|
| Reads | `01_Cluster/c_data/01_cluster.xlsx`, `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `01_Trunk/06_Modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/04_interactions.xlsx`, `b_reports/04_interactions_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/04_Interactions/a_script/04_interactions.R` |

## Method

- Cluster x time: clusters made at one input, tested on a change that comes later and was not used
  to make them. T1 clusters on training and acute change; T2 and training-change clusters on acute
  change. limma, `change ~ cluster`.
- Cluster x phenotype: does a feature's training change track a trait differently in the two
  clusters? limma, `change ~ cluster * trait`, the interaction term, once per trait.

The three-way cluster x time x phenotype model is not fitted: with 14 or 15 people it leaves too
few residual degrees of freedom. BH within each screen.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per screen: features tested, at p < 0.05, at BH < 0.05. |
| `cluster_time` | Cluster 2 minus cluster 1 on the later change: effect, t, p, BH. |
| `cluster_phenotype` | The interaction term: cluster 2's slope minus cluster 1's, t, p, BH. |
