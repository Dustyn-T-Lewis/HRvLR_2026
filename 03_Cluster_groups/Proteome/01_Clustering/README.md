# 01_Clustering

Does the proteome alone, with no phenotype or group label, split these people into groups, and if
two groups are imposed, who is in each?

## The idea

The features are the 01_Processing module eigengenes (`01_Processing/06_Module_Building`): one summary value
per co-expression module per biopsy. The distance between two people then rests on the modules,
not on every protein, so thousands of noisy single proteins do not set it.

People are clustered five times, once per input:

| Input | Each person's values | Who is clustered |
|---|---|---|
| `T1` | eigengenes at T1 | people with a T1 biopsy |
| `T2` | eigengenes at T2 | people with a T2 biopsy |
| `T3` | eigengenes at T3 | people with a T3 biopsy |
| `training` | T2 minus T1 | people with both (`subject_change()` in `R/helpers.R`) |
| `acute` | T3 minus T2 | people with both |

Within each input every eigengene is z-scored across people, so each module weighs the same in
the distance. HR/LR is never used to form the clusters; the adjusted Rand index against HR/LR
(0 for chance agreement, 1 for the same split) is reported as a description.

The K rule is the one the muscle arcs use, explained in [`../../README.md`](../../README.md):
PAM, the gap statistic for K = 1 to 3 with the one-standard-error rule, mean silhouette for
K = 2 and 3, and bootstrap Jaccard per cluster. The gap statistic selects K = 1 at every input,
and the workbook reports that as the result. K = 2 is imposed at every input (`imposed_k`), with
the recorded reason: so the exploratory arc can compare the two halves the proteome suggests; read
each split with its bootstrap Jaccard.

For each input the script then writes the 01_Processing inputs again with the clusters as groups `C1`
and `C2`, so the later stages run the HR/LR code on them.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Inputs | `subject_values()`, `subject_change()` in `R/helpers.R` | Eigengenes per person at T1, T2, T3, and their T2 - T1 and T3 - T2 changes. |
| Scaling | `scale(t(m))` | People as rows, each eigengene z-scored within the input. |
| Gap statistic | `cluster::clusGap(FUNcluster = pam, K.max = 3, B = 500)`, `cluster::maxSE(method = "Tibs2001SEmax")` | Gap and SE for K = 1 to 3 per input; the one-SE rule picks K. Seed 1. |
| Solutions | `cluster::pam()` for K = 2 and 3 | Cluster per person, mean silhouette width, sizes. |
| Agreement with HR/LR | `mclust::adjustedRandIndex()` | Per input and K. Description only. |
| Stability | `fpc::clusterboot(clustermethod = fpc::claraCBI, usepam = TRUE, B = 500, seed = 1)` | Mean Jaccard per cluster over 500 bootstrap resamples of people, K = 2 and 3. |
| Choice | `coalesce(imposed_k[[input]], k_gap)` | The K used downstream per input. |
| Relabelled inputs | `relabel()` in the script | Per input: `arm` becomes `C1` or `C2` (PAM's cluster number), `group` becomes `C1_T1` and so on, and the samples of anyone not clustered at that input drop out. |
| Space | `prcomp()` on the scaled eigengenes | PC1 and PC2 per person and input, for Figure 1 panel a. A projection, no test. |

## How to read the outputs

- `k_gap` is the result: at each input the proteome forms one spread of people, not groups.
- C1 and C2 are imposed halves. The number is PAM's cluster number and carries no direction:
  `C1` is not the higher or the better group, and C1 at one input is unrelated to C1 at another.
- Read each split with its Jaccard: 0.75 or more stable, 0.6 to 0.75 a pattern with doubtful
  membership, below 0.6 not to be trusted. A later contrast between unstable clusters compares
  groups that a resample would redraw.
- A cluster made from one time point's proteome differs in that proteome by construction. Later
  stages tag every comparison circular, caution or clean against the input
  ([`../README.md`](../README.md)).
- The people differ between inputs: someone missing a biopsy is absent from the inputs that need
  it.

## Files

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds` (features), `00_Input/phenotype.csv` (HR/LR, for description only); for the relabelled inputs, `01_Processing/03_Quantification/c_data/proteins.rds`, `DAList_corrected.rds`, `DAList_limpa.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds` |
| Writes | `c_data/clustering.xlsx`, `c_data/inputs/<input>/` (read by stages 02 to 07 and the F1 script), `b_reports/clustering_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/01_Clustering/a_script/clustering.R` |

The five inputs are listed in the script (`inputs_used`) and again in `../R/inputs.R`, which the
later stages read; change both together. The outputs are not tracked in git.

```
01_Clustering/
├── a_script/clustering.R
├── b_reports/clustering_report.pdf
└── c_data/
    ├── clustering.xlsx
    └── inputs/
        └── <input>/                   one each: T1, T2, T3, training, acute
            ├── proteins.rds           limpa object, metadata relabelled C1/C2
            ├── DAList_corrected.rds   blood-corrected complete matrix
            ├── DAList_limpa.rds       limpa's complete matrix
            ├── singscore.rds          gene-set scores, subset to the kept samples
            ├── modules.rds            module eigengenes, metadata relabelled
            └── groups.rds             group_of (subject to C1/C2), input
```

## Workbook: `clustering.xlsx`

`read_me` lists the sheets; `columns` defines every column. Every result sheet starts with `input`.

| Sheet | One row is | Key columns |
|---|---|---|
| `choice` | an input | `n_subjects`, `n_features`, `k_gap`, `k_imposed`, `imposed_reason`, `k_used`, `grouped` |
| `diagnostics` | an input and K from 1 to 3 | `gap`, `SE.sim`, `mean_silhouette`, `sizes`, `hr_lr_ari`, `min_jaccard`, `jaccard` |
| `stability` | an input and a cluster at K = 2 or 3 | `k`, `cluster`, `jaccard` |
| `assignments` | an input and a person | `arm` (HR/LR, description), `k2`, `k3`, `cluster` |
| `space` | an input and a person | `arm`, `cluster`, `PC1`, `PC2` |
| `space_variance` | an input and a component | `variance_explained` |

`08_Cluster_Phenotype` and the F1 script read `assignments`, `space`, `space_variance`, `choice`
and `stability` from here.

## Report: `clustering_report.pdf`

1. Gap statistic against K with its SE, one panel per input; the dashed line marks the K the
   one-SE rule selects.
