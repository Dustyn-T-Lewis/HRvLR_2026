# 03_Profile_Clusters

Which proteins share a pattern of trait links, with every person pooled?

## The idea

`01_Trait_Associations/01_Proteins` gives every protein a limma t against every trait in every
window: 24 change screens (three change windows by eight change traits) and 18 level screens
(three level windows by six trait levels), 42 in all. Laid side by side, those 42 t values are the
protein's association profile: a fingerprint of which traits it rises with, which it falls with,
and in which windows.

A profile cluster is a group of proteins whose fingerprints have the same shape. Distance between
two proteins is 1 minus the Pearson correlation of their profiles, so two proteins that rise and
fall with the same traits sit close even if one has larger t values throughout. Average-linkage
clustering builds a tree from those distances, and `dynamicTreeCut` cuts the tree into groups of
20 or more proteins, following its branches rather than cutting at one height. Proteins on no
branch large enough stay unassigned (group 0) and leave the tables.

A second grouping works per trait. For the `training` and `baseline_level` windows, it takes each
trait's proteins at p < 0.05 and groups them by how they co-vary across people: distance is
1 minus the Pearson correlation of their per-person values (change T1 to T2, or T1 level), cut the
same way into groups of 8 or more. This splits one trait's hit list into blocks that move
together in the same people.

Each group is then described: its size, mean t, the share of members with positive t, the module
most members belong to (`01_Trunk/06_Module_Building`), its top members by |t|, and the gene set
most over-represented among its members.

The groups are built from the same t values they summarise. Their mean t and their enrichment
describe the groups; they do not test whether the grouping is real. The mean silhouette width
does report how separate the groups are: near 1, each protein sits much closer to its own group
than to the next; near 0, the groups blur together.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Profiles | `readRDS()` of `associate.rds`, `pivot_wider()` | Protein by screen matrix of limma t. Proteins without a t in every screen drop out (`na.omit()`). |
| Profile tree | `hclust(as.dist(1 - cor(t(...))), "average")` | Distance between profiles. |
| Profile cut | `dynamicTreeCut::cutreeDynamic(minClusterSize = 20, deepSplit = 2)` | Group 0 dropped. |
| Separation | `cluster::silhouette()` | Mean width over the assigned proteins; `NA` when fewer than two groups form. |
| Per-trait groups | `subject_change()`, `subject_values()` in `R/helpers.R`; the same tree and cut with `minClusterSize = 8` | On each trait's p < 0.05 proteins in `training` and `baseline_level`, using their per-person values from `DAList_corrected.rds`. |
| Module | `modules.rds` `membership` | `main_module`: the most common module among the group's proteins. |
| Over-representation | `clusterProfiler::enricher(minGSSize = 15, maxGSSize = 500)` | Group members against a universe: every profiled protein for profile groups, every protein tested in the window for per-trait groups. Gene symbols from the representative protein per gene; the set with the smallest p is kept, with its p and BH. |
| Mean t | `summarise()` | Profile groups: the mean over members of each protein's mean t across all 42 screens. Per-trait groups: the mean of members' t against that trait. |

## How to read the outputs

- The 42 screens are not independent. `baseline` and `baseline_level` share the same T1 values;
  the three change windows share the same traits; the composite is built from fibre-area change.
  A profile group can be driven by a block of related screens.
- A low mean silhouette says the groups are cuts through a continuum, not separate families.
- `top_set` is the best of many sets tried for each group. Its `top_set_fdr` is BH across the sets
  tried for that group only, not across groups.
- A profile group's `mean_t` averages across all screens, so a group that rises with one trait
  and falls with another can show a mean near zero. Read `profile_means` for the shape.
- Per-trait groups exist only for `training` and `baseline_level`, and only for traits with enough
  hits to form a group of 8.

## Files

| | |
|---|---|
| Reads | `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/associate.rds` (t and p per protein, window and trait), `01_Trunk/03_Quantification/c_data/DAList_corrected.rds` (per-person values and gene symbols), `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` (sets and protein-to-gene map), `01_Trunk/06_Module_Building/c_data/modules.rds` (module membership) |
| Writes | `c_data/profile_clusters.xlsx`, `b_reports/profile_clusters_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/03_Profile_Clusters/a_script/profile_clusters.R` |

```
03_Profile_Clusters/
├── a_script/profile_clusters.R
├── b_reports/profile_clusters_report.pdf
└── c_data/profile_clusters.xlsx
```

## Workbook: `profile_clusters.xlsx`

`read_me` lists the sheets, with the profile cut's mean silhouette in the `profile_groups` row;
`columns` defines every column.

| Sheet | One row is | Key columns |
|---|---|---|
| `profile_groups` | a profile group | `group`, `n_proteins`, `mean_t`, `share_positive`, `main_module`, `proteins` (up to 15 member gene symbols, largest \|t\| first), `top_set`, `top_set_p`, `top_set_fdr` |
| `profile_members` | a protein | `feature`, `label`, `group` |
| `profile_means` | a profile group and screen | `group`, `screen` (window and trait), `mean_t` |
| `trait_groups` | a per-trait group | `window`, `trait`, `group`, `n_proteins`, `mean_t`, `share_positive`, `main_module`, `proteins`, `top_set`, `top_set_p`, `top_set_fdr`, `silhouette`, `n_hits` (the trait's p < 0.05 proteins) |

## Report: `profile_clusters_report.pdf`

1. Mean t per profile group and screen as a heatmap; the subtitle gives the proteins grouped, the
   number of groups and the mean silhouette.
2. Per-trait groups: one bar per window and trait, one segment per group, filled by the group's
   mean t against the trait.
