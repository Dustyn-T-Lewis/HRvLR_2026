# 01_Clustering

Do people fall into groups by how much their muscle grew, and if two groups are imposed, who is in each?

| | |
|---|---|
| Reads | `00_Input/phenotype.csv` |
| Writes | `c_data/clustering.xlsx`, `c_data/groups.rds`, `b_reports/clustering_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/01_Clustering/a_script/clustering.R` |

## Key result

- 16 people clustered; cluster sizes 8 and 8; bootstrap Jaccard 0.89 and 0.88.
- Gap statistic selects K = 1; K = 2 imposed by decision (2026-10-06). It reproduces the HR/LR split exactly, so this grouping repeats that comparison under a size-cluster label.

## Method

- PAM on the z-scored % changes over training in fibre area (mixed, type I, type II) and muscle area, complete cases.
- K from the gap statistic with the one-SE rule (`cluster::clusGap`, K = 1 to 3, 500 references); stability as bootstrap Jaccard per cluster (`fpc::clusterboot`).
- Agreement with HR/LR as the adjusted Rand index (`mclust::adjustedRandIndex`), a description.
- When K = 2 is used, `groups.rds` names the clusters Hi and Lo by mean fibre-area change; `02_Design_Contrasts` reads it.

## Read with care

- The imposed K = 2 is a decision, not a finding; Hi and Lo differ in size change by construction.

## Sheets

`clustering.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `choice` | K selected by the gap statistic, any imposed K and its reason, and the K used downstream. |
| `diagnostics` | Per K: gap, SE, mean silhouette, group sizes, HR/LR adjusted Rand index, bootstrap Jaccard. |
| `stability` | Per K and cluster: mean bootstrap Jaccard (fpc::clusterboot, PAM). |
| `assignments` | Per subject: the traits, the K = 2 and K = 3 clusters, and the cluster used downstream. |
| `space` | Per subject: HR/LR arm, group, cluster, and PC1 and PC2 of the scaled traits. |
| `space_variance` | Share of variance on PC1 and PC2 of the scaled traits. |
