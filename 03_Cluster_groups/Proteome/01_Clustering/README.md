# 01_Clustering

| | |
|---|---|
| Reads | `01_Trunk/06_Module_Building/c_data/modules.rds`, `00_Input/phenotype.csv` (HR/LR, for description only) |
| Writes | `c_data/01_cluster.xlsx`, `c_data/inputs/<input>/` (relabelled trunk inputs), `b_reports/01_cluster_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/01_Clustering/a_script/01_cluster.R` |

## Method

Per input, PAM on the z-scored module eigengenes. K diagnostics as in the phenotype clusters: the
gap statistic for K = 1 to 3 with the one-standard-error rule (`cluster::clusGap`,
`maxSE(method = "Tibs2001SEmax")`), mean silhouette for K = 2 and 3, group sizes, and bootstrap
Jaccard stability (`fpc::clusterboot`, 500 resamples). The adjusted Rand index against HR/LR
(`mclust::adjustedRandIndex`) is a description. `imposed_k` and `imposed_reason` fix K per input.

For each input, the trunk's protein, limpa, blood-corrected, singscore and module objects are
written to `c_data/inputs/<input>/` with the clusters as groups C1 and C2; people not clustered at
that input drop out. The stages after this one read those copies.

## Workbook

| Sheet | Holds |
|---|---|
| `choice` | Per input: subjects, features, K selected, imposed K and reason, K used. |
| `diagnostics` | Per input and K: gap, SE, mean silhouette, sizes, HR/LR ARI, Jaccard. |
| `stability` | Per input, K and cluster: mean bootstrap Jaccard. |
| `assignments` | Per input and subject: HR/LR, the K = 2 and K = 3 clusters, the cluster used. |
