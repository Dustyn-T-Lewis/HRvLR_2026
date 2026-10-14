# 03_Cluster_groups

If the groups came from the data rather than from the study's HR/LR labels, would the results
change?

`02_HRvLR_groups` compares the study's High and Low Responders. Here the data propose the groups,
three ways, and the same steps run on each grouping:

| Folder | People clustered on | Groups |
|---|---|---|
| [`Muscle_size`](Muscle_size/README.md) | % change over training in fibre area (mixed, type I, type II) and muscle area | `Hi`, `Lo` |
| [`Muscle_strength`](Muscle_strength/README.md) | % change over training in leg-press and leg-extension 1RM | `Hi`, `Lo` |
| [`Proteome`](Proteome/README.md) | module eigengenes at T1, T2, T3, training change and acute change | `C1`, `C2` per input |

## Clustering

Each trait or eigengene is z-scored across people, then PAM (`cluster::pam`) picks K real people
as centres and assigns everyone to the nearest. PAM always returns K groups, so three diagnostics
judge the answer:

| Diagnostic | Function | Measures |
|---|---|---|
| Gap statistic, one-SE rule | `cluster::clusGap` (K = 1 to 3, 500 reference sets), `cluster::maxSE(method = "Tibs2001SEmax")` | whether K clusters are tighter than in data with no groups; K = 1 is a possible answer |
| Silhouette | `pam()$silinfo$avg.width` | how much closer people sit to their own cluster than to the next |
| Bootstrap Jaccard | `fpc::clusterboot`, 500 resamples | how often a cluster comes back when people are resampled: 0.75 or more stable, 0.6 to 0.75 doubtful, below 0.6 not to be trusted |

The gap statistic selects K = 1 for size, for strength and for the proteome at every input. K = 2
is then imposed by hand (`imposed_k`), with the reason recorded in each `clustering.xlsx`
(`choice`). An imposed split cuts one spread of people in two; every result downstream of it is
exploratory.

## How the groups reach the shared steps

`01_Clustering` writes `groups.rds` (one per input for `Proteome`): who is in each group, the two
group names and the traits that defined them. `02_Design_Contrasts` reads it and records the
groups in `design.rds`; every later step reads the 01_Processing files through `subset_for_design()` in
`00_Setup/Shared_functions/analysis_reporting_helpers.R`, which keeps the grouped people's samples and puts the cluster names in `arm` and
`group`. The steps are the same functions `02_HRvLR_groups` calls, so the files, sheets and
columns match. The contrasts keep their shape: `Training_Hi`, `Baseline_HivLo`,
`Training_Interaction`.

`05_Module_Preservation` runs only for HR/LR; `08_Cluster_Phenotype` only for `Proteome`.
