# 03_Profile_Clusters

Which proteins share a pattern of trait links, with every person pooled?

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/associate.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/profile_clusters.xlsx`, `b_reports/profile_clusters_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/03_Profile_Clusters/a_script/profile_clusters.R` |

## Key result

- 46 profile groups hold 1,900 proteins; the largest holds 117.
- 63 within-trait groups over 14 trait and window pairs.

## Method

- Each protein's limma t across the 42 window-and-trait screens; average-linkage tree on 1 - Pearson r, cut by `dynamicTreeCut::cutreeDynamic(minClusterSize = 20)`.
- Within-trait groups: the same on each trait's p < 0.05 proteins in `training` and `baseline_level` (minimum 8).
- Each group's most common module and its best over-represented gene set (`clusterProfiler::enricher()`).

## Read with care

- The 42 screens are not independent; a group can follow a block of related screens.
- `top_set` is the best of many sets tried for that group.

## Sheets

`profile_clusters.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `profile_groups` | Profile groups: size, mean t, main module, top over-represented set; mean silhouette 0.033. |
| `profile_members` | Each protein's profile group. |
| `profile_means` | Mean limma t per profile group and screen. |
| `trait_groups` | Per trait and window: co-variation groups among p < 0.05 proteins, with silhouette and ORA. |
