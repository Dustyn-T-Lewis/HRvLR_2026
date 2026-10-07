# Muscle_size

Do people grouped by how much their muscle grew differ in their proteome and its response?

## What was clustered

Each person's % change over training (T1 to T2) in four size traits: fibre cross-sectional area
of mixed, type I and type II fibres (`pct_fcsa_mixed`, `pct_fcsa_type1`, `pct_fcsa_type2`) and
muscle cross-sectional area (`pct_mcsa`). Only people with all four are clustered; anyone missing
one has no group, and their samples drop out of every later stage. No proteome value enters the
clustering.

The gap statistic selects K = 1. K = 2 is imposed with this recorded reason: it reproduces the
HR/LR split exactly, so this arc repeats that comparison under a size-cluster label. The cluster with
the higher mean `pct_fcsa_mixed` is `Hi`, the other `Lo`.

The method (PAM, gap statistic, silhouette, bootstrap Jaccard) is explained in
[`../README.md`](../README.md); the details for this arc are in
[`01_Clustering`](01_Clustering/README.md).

## How to read it

- Hi and Lo were defined by size change, so any size-change difference between them is built in.
  `Baseline_HivLo` is the floor, as `Baseline_HRvLR` is in the template.
- With K = 2 imposed over a gap-statistic K = 1, the groups are halves of a spread, not two kinds
  of people. Read the results as exploratory.
- Where the split matches HR/LR, results here repeat `02_HRvLR_groups` under new names.
  `01_Clustering/c_data/clustering.xlsx` (`space`) lists each person's HR/LR arm beside their
  group.

## Stages

| Stage | Asks | Copy of |
|---|---|---|
| [`01_Clustering`](01_Clustering/README.md) | how people split on size change | new: K diagnostics, Hi/Lo, relabelled inputs |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | which comparisons | `02_HRvLR_groups/01_Design_Contrasts` |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | which proteins differ | `02_HRvLR_groups/02_Protein_Contrasts` |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | which gene sets shift | `02_HRvLR_groups/03_Set_Contrasts` |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | which modules shift | `02_HRvLR_groups/04_Module_Contrasts` |
| [`06_Classification/01_Proteins`](06_Classification/01_Proteins/README.md) | can one protein tell Hi from Lo | `02_HRvLR_groups/06_Classification/01_Proteins` |
| [`06_Classification/02_Sets`](06_Classification/02_Sets/README.md) | can one gene-set score tell Hi from Lo | `02_HRvLR_groups/06_Classification/02_Sets` |
| [`06_Classification/03_Modules`](06_Classification/03_Modules/README.md) | can one module tell Hi from Lo | `02_HRvLR_groups/06_Classification/03_Modules` |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by group | `02_HRvLR_groups/07_Trait_Associations` |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the arc on a page | `02_HRvLR_groups/Figures/F1_Overview` |

Stages 02 to 07 read `01_Clustering/c_data/inputs/` in place of the trunk's files. `run_all.R`
runs 01 to 07 in order, then the F1 script.

## Files

```
Muscle_size/
├── 01_Clustering/
│   ├── a_script/clustering.R
│   ├── b_reports/clustering_report.pdf
│   └── c_data/
│       ├── clustering.xlsx
│       └── inputs/              trunk inputs relabelled Hi/Lo, read by 02 to 07
├── 02_Design_Contrasts/
├── 03_Protein_Contrasts/
├── 04_Set_Contrasts/
├── 05_Module_Contrasts/
├── 06_Classification/
│   ├── 01_Proteins/
│   ├── 02_Sets/
│   └── 03_Modules/
├── 07_Trait_Associations/
└── Figures/F1_Overview/
```
