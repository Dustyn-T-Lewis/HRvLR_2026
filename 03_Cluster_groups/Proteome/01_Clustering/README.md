# 01_Clustering

Does the proteome alone, with no phenotype or group label, split these people into groups?

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/clustering.xlsx`, `b_reports/clustering_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/01_Clustering/a_script/clustering.R` |

## Key result

- The gap statistic selects K = 1 at every input; K = 2 is imposed so the two halves can be compared.

| Input | People | C1/C2 | Jaccard C1/C2 |
|---|---|---|---|
| T1 | 15 | 9/6 | 0.75/0.66 |
| T2 | 15 | 7/8 | 0.62/0.65 |
| T3 | 15 | 5/10 | 0.63/0.84 |
| training | 14 | 7/7 | 0.67/0.56 |
| acute | 15 | 4/11 | 0.65/0.86 |

## Method

- PAM on the 15 module eigengenes per person, z-scored, at five inputs: T1, T2 and T3 levels, training change (T2 - T1) and acute change (T3 - T2).
- K from the gap statistic with the one-SE rule; stability as bootstrap Jaccard; agreement with HR/LR as the adjusted Rand index.
- `groups/<input>.rds` holds each input's clusters, C1 and C2 (PAM's numbering, no direction).

## Read with care

- A cluster built from one time point differs in that proteome by construction; later steps tag every comparison circular, caution or clean (`R/inputs.R`).
- C1 at one input is unrelated to C1 at another.

## Sheets

`clustering.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `choice` | Per input: subjects, features, K selected by the gap statistic, any imposed K and reason. |
| `diagnostics` | Per input and K: gap, SE, mean silhouette, sizes, HR/LR adjusted Rand index, Jaccard. |
| `stability` | Per input, K and cluster: mean bootstrap Jaccard (fpc::clusterboot, PAM). |
| `assignments` | Per input and subject: HR/LR (for description), the K = 2 and K = 3 clusters, cluster used. |
| `space` | Per input and subject: HR/LR, cluster, and PC1 and PC2 of the z-scored eigengenes. |
| `space_variance` | Per input: share of variance on PC1 and PC2 of the z-scored eigengenes. |
