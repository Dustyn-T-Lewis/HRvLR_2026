# 04_Profile_Clusters

Proteins grouped by how they relate to phenotype, every subject pooled.

| | |
|---|---|
| Reads | `01_Associate/01_Proteins/c_data/associate.rds`, `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/06_Modules/c_data/modules.rds` |
| Writes | `c_data/04_profile_clusters.xlsx`, `b_reports/04_profile_clusters_figures.pdf` |
| Run | `Rscript 04_Continuous_no_group/04_Profile_Clusters/a_script/04_profile_clusters.R` |

## Method

Profiles: each protein described by its limma t against every trait in every pairing; average
linkage on 1 - Pearson r, cut by `dynamicTreeCut` (groups of 20 or more). Per trait: for training
change and the T1 level, the proteins at p < 0.05 grouped by how they co-vary across people
(groups of 8 or more). Each group gets its size, mean t, main module and the top set from
`clusterProfiler::enricher` against the proteins it was drawn from. The groups come from the same
t values they summarise, so the enrichment describes them and does not test them; the mean
silhouette says how separate they are.

## Workbook

| Sheet | Holds |
|---|---|
| `profile_groups` | Per profile group: size, mean t, share positive, main module, top set and its p and BH. |
| `profile_members` | Each protein's profile group. |
| `profile_means` | Mean t per profile group and screen. |
| `trait_groups` | Per trait and window: co-variation groups with silhouette and ORA. |
