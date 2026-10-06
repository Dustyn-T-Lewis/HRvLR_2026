# 01_Cluster

People clustered on their muscle size response before any proteome result is read.

| | |
|---|---|
| Reads | `00_Input/phenotype.csv`; for the relabelled inputs, `01_Trunk/03_Quantify`, `05_Scores` and `06_Modules` outputs |
| Writes | `c_data/01_cluster.xlsx`, `c_data/inputs/` (proteins, matrices, singscore, modules, groups), `b_reports/01_cluster_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/01_Cluster/a_script/01_cluster.R` |

## Method

PAM on the z-scored % changes (% change in fibre area (mixed, type I, type II) and mCSA). K diagnostics: the gap statistic for K = 1 to 3 with
the one-standard-error rule (`cluster::clusGap`, `maxSE(method = "Tibs2001SEmax")`), mean
silhouette for K = 2 and 3, group sizes, and bootstrap Jaccard stability per cluster
(`fpc::clusterboot`, 500 resamples). `imposed_k` and `imposed_reason` at the top of the script fix
K by hand; both K values are recorded.

The clusters become groups `Hi` and `Lo` by the mean of the first trait. The script writes the
trunk's inputs again with `arm` and `group` set to those groups and any unassigned subject's
samples removed, so the later stages run the HR/LR code unchanged.

## Workbook

| Sheet | Holds |
|---|---|
| `choice` | K selected by the gap statistic, imposed K and reason, K used. |
| `diagnostics` | Per K: gap, SE, mean silhouette, sizes, Jaccard per cluster. |
| `stability` | Per K and cluster: mean bootstrap Jaccard. |
| `assignments` | Per subject: the traits, the K = 2 and K = 3 clusters, the cluster used. |
