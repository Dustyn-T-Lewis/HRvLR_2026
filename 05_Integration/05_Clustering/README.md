# 05_Clustering

Looks for structure without the HR/LR labels, and judges every structure against what shuffled data
give.

| | |
|---|---|
| Reads | `01_Preprocess/03_Quantify/c_data/DAList_corrected.rds`, `02_Proteins/04_Associate/c_data/associate.rds`, `03_Pathways/01_Gene_Sets/c_data/gene_sets.rds`, `04_Modules/01_Build/c_data/modules.rds`, `05_Integration/04_Continuous/c_data/shuffled_t.rds`, `00_Input/phenotype.csv`, the STRING v12 files |
| Writes | `c_data/05_clustering.xlsx`, `b_reports/05_clustering_figures.pdf` |
| Run | `Rscript 05_Integration/05_Clustering/a_script/05_clustering.R` |

## Method

- Profiles: every protein described by its limma t against every trait in every pairing; average
  linkage on 1 - Pearson r, cut by `dynamicTreeCut`. Mean silhouette against the same clustering of
  the 100 shuffled profiles.
- Per trait: the proteins at p < 0.05 for one trait (training change and T1 level), grouped by how
  they co-vary across people; each group annotated with its direction, main module, top pathway
  (over-representation, `clusterProfiler::enricher`) and STRING density; silhouette against 50
  shuffled lists of the same length. Per trait family, each protein's t across the pairings, PAM
  with k from 2 to 6, against the 100 shuffles.
- People: the subjects clustered on module eigengenes and on the 500 most variable proteins, per
  input (T1, T2, T3, training change, acute change), PAM with k from 2 to 4 on 1 - Pearson r.
  Silhouette against 200 shuffles of each feature across people; then agreement with HR/LR
  (adjusted Rand index, Fisher) and with each trait (Kruskal-Wallis). The feature-wise shuffle
  breaks all correlation between proteins, so beating it shows only that people differ along
  correlated axes; the agreement columns say whether those axes are response.
- Groups: the subjects split in two within each pathway's and each module's proteins, per input.
  Searching about 1,400 groups makes some split look clean by chance, so the best silhouette and the
  best agreement with HR/LR are judged against the best of 100 shuffled searches.
- Blocks: two-way hierarchical blocks of the profile matrix (mean abs t against the largest block of
  the shuffled matrices) and `biclust`'s plaid model (biclusters found against 20 shuffles).

## Workbook

| Sheet | Holds |
|---|---|
| `summary` | Per analysis and screen: the statistic, its shuffle null and p. |
| `profile_groups`, `profile_members`, `profile_means` | Profile clusters: size, main module, top pathway; members; mean t per screen. |
| `trait_groups` | Per trait, groups of co-varying associated proteins with their annotation. |
| `trait_profiles`, `trait_profile_means` | Per trait family, clusters across pairings: annotation and mean t. |
| `people` | Per feature set, input and k: silhouette, null, group sizes, agreement with HR/LR and the traits, members. |
| `groups`, `group_summary` | Per pathway or module and input: silhouette, ARI, smallest trait p; per input the best of each against its null. |
| `blocks`, `plaid` | Two-way blocks with their p; plaid biclusters. |
