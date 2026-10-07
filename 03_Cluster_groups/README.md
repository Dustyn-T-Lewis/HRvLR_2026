# 03_Cluster_groups

If the groups came from the data rather than from the study's median split, would the proteome
story change?

## The idea

`02_HRvLR_groups` compares High and Low Responders, a split the study made at the median of a
composite hypertrophy score. That split is one choice among many. This arc lets the data propose
the groups instead, three ways, and then reruns the HR/LR stages on each grouping:

| Sub-pipeline | People clustered on | Groups |
|---|---|---|
| [`Muscle_size`](Muscle_size/README.md) | % change over training in fibre area (mixed, type I, type II) and mCSA | `Hi`, `Lo` |
| [`Muscle_strength`](Muscle_strength/README.md) | % change over training in leg-press and leg-extension 1RM | `Hi`, `Lo` |
| [`Proteome`](Proteome/README.md) | module eigengenes, separately at T1, T2, T3, training change (T2 - T1) and acute change (T3 - T2) | `C1`, `C2` per input |

The two muscle arcs ask whether a size or a strength response, taken on its own, marks a
different proteome. The proteome arc turns the question round: it groups people by their proteome
with no phenotype or group label, then asks what else differs between those groups.

## Clustering from the ground up

Clustering puts people into K groups so that people in a group are close and people in different
groups are far apart. Each trait or eigengene is first z-scored across people (mean 0, SD 1), so
one with large numbers does not dominate the distance.

PAM (partitioning around medoids, `cluster::pam`) picks K real people as centres, the medoids,
and assigns everyone else to the nearest one. It chooses the medoids that make the total distance
to them smallest. Using real people as centres makes it less sensitive to one extreme person than
k-means.

PAM always returns K groups when asked for K, whether or not the data hold groups. Three
diagnostics judge the answer:

| Diagnostic | Function | What it measures |
|---|---|---|
| Gap statistic | `cluster::clusGap`, 500 reference sets, K = 1 to 3 | How much tighter K clusters are in the real data than in data with no groups. The reference data are drawn uniformly over the real data's range along its principal axes, 500 times; the gap for K is the mean log within-cluster dispersion of the reference minus that of the real data, with a simulation SE. |
| One-SE rule | `cluster::maxSE(method = "Tibs2001SEmax")` | Picks the smallest K whose gap is at least the next K's gap minus the next K's SE (Tibshirani, Walther & Hastie 2001). K = 1 is a candidate, so the rule can say there are no groups. |
| Silhouette | `pam()$silinfo$avg.width`, K = 2 and 3 | For each person, how much closer they sit to their own cluster than to the nearest other one, from -1 to 1; averaged over people. It needs two clusters, so it cannot judge K = 1 against K = 2. |
| Bootstrap Jaccard | `fpc::clusterboot`, PAM, 500 resamples, K = 2 and 3 | Resamples people with replacement, reclusters, and finds each original cluster's best match. Jaccard is the people in both divided by the people in either; the mean over resamples is the cluster's stability (Hennig 2007). |

The fpc guidance reads a mean Jaccard of 0.75 or more as a stable cluster, 0.6 to 0.75 as a
pattern whose exact membership is doubtful, and below 0.6 as a cluster not to be trusted.

## How K is chosen

The gap statistic selects K = 1 for size, for strength and for the proteome at every input: by
that rule these people form one group. Each `01_Clustering` workbook reports this as the result
(`choice`, `k_gap`).

K = 2 is then imposed by hand in each clustering script (`imposed_k`), with the reason written to
`imposed_reason` and carried into `choice`:

| Arc | Recorded reason |
|---|---|
| `Muscle_size` | K = 2 reproduces the HR/LR split exactly, so the arc repeats that comparison under a size-cluster label. |
| `Muscle_strength` | Both clusters are stable by bootstrap Jaccard and the split is independent of HR/LR. Exploratory. |
| `Proteome` | Imposed so the arc can compare the two halves the proteome suggests; read each split with its bootstrap Jaccard. |

An imposed K = 2 cuts a spread of people in two. It does not show that two kinds of people exist.
Read every result downstream of a clustering as exploratory, and read each split with its
Jaccard: a contrast between two unstable clusters compares groups that a resample would redraw.

## How the copies reuse the HR/LR template

`01_Clustering` in each arc writes relabelled copies of the 01_Processing inputs to
`01_Clustering/c_data/inputs/` (one folder per input for `Proteome`): `proteins.rds`,
`DAList_corrected.rds`, `DAList_limpa.rds`, `singscore.rds`, `modules.rds` and `groups.rds`. In
each copy `arm` holds the cluster and `group` the cluster and time point (`Hi_T1`, `C2_T3`), and
the samples of anyone without a cluster are dropped.

The stages from `02_Design_Contrasts` on are copies of the HR/LR stages that read those files
instead of 01_Processing's, with HR and LR renamed. The contrasts keep their shape: `Training_Hi`,
`Baseline_HivLo`, `Training_Interaction` and so on.

| Cluster stage | Copy of |
|---|---|
| `02_Design_Contrasts` | `02_HRvLR_groups/01_Design_Contrasts` |
| `03_Protein_Contrasts` | `02_HRvLR_groups/02_Protein_Contrasts` |
| `04_Set_Contrasts` | `02_HRvLR_groups/03_Set_Contrasts` |
| `05_Module_Contrasts` | `02_HRvLR_groups/04_Module_Contrasts` |
| `06_Classification/01_Proteins`, `02_Sets`, `03_Modules` | `02_HRvLR_groups/06_Classification/...` |
| `07_Trait_Associations` | `02_HRvLR_groups/07_Trait_Associations` |

`05_Module_Preservation` has no copy. The `Proteome` copies wrap the template in a loop over the
five inputs and add `08_Cluster_Phenotype`. Each copy is its own script, so a fix to a template
stage has to be made in every copy by hand.

`Figures/F1_Overview` in each arc draws the HR/LR Figure 1 on the clusters, with a clustering
panel in place of the design panel. `run_all.R` runs every stage in this arc in order (size,
strength, proteome), each arc's F1 script after its last stage.

## Files

```
03_Cluster_groups/
├── Muscle_size/
│   ├── 01_Clustering/           K diagnostics, Hi/Lo, relabelled inputs
│   ├── 02_Design_Contrasts/ ... 07_Trait_Associations/
│   └── Figures/F1_Overview/
├── Muscle_strength/
│   ├── 01_Clustering/
│   ├── 02_Design_Contrasts/ ... 07_Trait_Associations/
│   └── Figures/F1_Overview/
└── Proteome/
    ├── R/inputs.R               the five inputs, validity(), stack_inputs()
    ├── 01_Clustering/           K diagnostics per input, C1/C2, relabelled inputs per input
    ├── 02_Design_Contrasts/ ... 07_Trait_Associations/   one folder per input in each
    ├── 08_Cluster_Phenotype/    do the proteome clusters differ in size or strength
    └── Figures/F1_Overview/     one figure per input
```

The muscle arcs' outputs are tracked in git. The `Proteome` arc's `b_reports/` and `c_data/`
folders are not (`.gitignore`); `run_all.R` and the F1 script regenerate them.
