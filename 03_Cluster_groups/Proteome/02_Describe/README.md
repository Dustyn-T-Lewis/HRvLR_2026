# 02_Describe

What separates each split, at the input that made it. Descriptive: the clusters were built from
these data, so the p-values say what the split looks like, not that it is real.

| | |
|---|---|
| Reads | `01_Cluster/c_data/01_cluster.xlsx`, `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `01_Trunk/06_Modules/c_data/modules.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_describe.xlsx`, `b_reports/02_describe_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/02_Describe/a_script/02_describe.R` |

## Method

Per input: each module's mean z-scored eigengene per cluster; limma on every protein at that input,
cluster 2 minus cluster 1; and `clusterProfiler::enricher` among the proteins higher and lower in
cluster 2 at p < 0.05, against all measured proteins.

## Workbook

| Sheet | Holds |
|---|---|
| `module_means` | Per input, module and cluster: mean z-scored eigengene. |
| `protein_difference` | Per input and protein: effect, t, p, BH. |
| `enrichment` | Per input and direction: set, overlap, p, BH. |
