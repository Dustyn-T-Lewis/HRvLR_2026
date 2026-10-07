# 01_Clustering

Do people fall into groups by how much stronger they became, and if two groups are imposed, who
is in each?

## The idea

Each person is described by two numbers: the % change from T1 to T2 in leg-press 1RM and in
leg-extension 1RM, `100 * (T2 - T1) / T1` (`phenotype_tables()` in `R/helpers.R`). Each trait
is z-scored across people so both weigh the same in the distance between two people. Only people
with both traits are clustered. The proteome is not read until the clusters are
fixed.

PAM finds K groups around K real people (medoids). It returns K groups whatever the data look
like, so the script asks separately whether the data hold groups at all: the gap statistic with
the one-standard-error rule compares K = 1, 2 and 3 against data with no groups, and can choose
K = 1. Silhouette and bootstrap Jaccard then describe the K = 2 and K = 3 solutions. The arc
README ([`../../README.md`](../../README.md)) explains each from the ground up.

`imposed_k` at the top of the script overrides the gap statistic, and `imposed_reason` records
why. Both the selected and the imposed K are written to the workbook. The gap statistic selects
K = 1 here; K = 2 is imposed because both clusters are stable by bootstrap Jaccard and the split
is independent of HR/LR. The reason itself calls the arc exploratory.

When the K used is 2, the script writes the 01_Processing inputs again with the clusters as the groups.
The later stages read these copies, so they run the HR/LR code unchanged.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Traits | `phenotype_tables()$pairings$training` in `R/helpers.R`, `tidyr::drop_na()`, `scale()` | % change over training in `pct_leg_press_1rm` and `pct_leg_ext_1rm`, complete cases only, z-scored. |
| Gap statistic | `cluster::clusGap(FUNcluster = pam, K.max = 3, B = 500)`, `cluster::maxSE(method = "Tibs2001SEmax")` | Gap and its SE for K = 1 to 3; the one-SE rule picks K. Seed 1. |
| Solutions | `cluster::pam()` for K = 2 and 3 | Cluster per person and mean silhouette width. |
| Agreement with HR/LR | `mclust::adjustedRandIndex()` | Per K, the clusters against HR/LR: 1 is the same split, 0 no more alike than chance. Description only. |
| Stability | `fpc::clusterboot(clustermethod = fpc::claraCBI, usepam = TRUE, B = 500, seed = 1)` | Mean Jaccard per cluster over 500 bootstrap resamples of people, for K = 2 and 3. |
| Choice | `coalesce(imposed_k, k_gap)` | The K used downstream: the imposed K, or the gap statistic's if none is imposed. |
| Labels | mean `pct_leg_press_1rm` per cluster | The cluster with the higher mean is `Hi`, the other `Lo`. |
| Relabelled inputs | `relabel()` in the script | Only when K used is 2. Rewrites `arm` as Hi/Lo and `group` as `Hi_T1` and so on in the 01_Processing protein, matrix, singscore and module objects, and drops the samples of anyone without a cluster. |
| Space | `prcomp()` on the scaled traits | PC1 and PC2 per person, for Figure 1 panel a. A projection, no test. |

## How to read the outputs

- `k_gap` is the result. K = 1 means the strength changes form one spread with no gap in it.
- The imposed K = 2 is a decision, not a finding. Hi and Lo are halves of that spread; every
  comparison between them downstream is exploratory.
- Read the K = 2 Jaccard per cluster: 0.75 or more is stable, 0.6 to 0.75 a pattern with doubtful
  membership, below 0.6 not to be trusted.
- Hi and Lo differ in strength change by construction. A strength difference between them later
  restates the split.
- Anyone without both traits has no cluster: an empty `cluster` in `assignments`, and no samples
  in `c_data/inputs/`.
- `space` carries each person's HR/LR arm beside their group, so how the split crosses HR/LR can
  be read off directly.

## Files

| | |
|---|---|
| Reads | `00_Input/phenotype.csv`; for the relabelled inputs, `01_Processing/03_Quantification/c_data/proteins.rds`, `DAList_corrected.rds`, `DAList_limpa.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/clustering.xlsx`, `c_data/inputs/` (read by stages 02 to 07 and the F1 script), `b_reports/clustering_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/01_Clustering/a_script/clustering.R` |

`input_manifest` lists `phenotype.csv` only; 01_Processing files are read for relabelling, not for
the clustering.

```
01_Clustering/
├── a_script/clustering.R
├── b_reports/clustering_report.pdf
└── c_data/
    ├── clustering.xlsx
    └── inputs/
        ├── proteins.rds           limpa object, metadata relabelled
        ├── DAList_corrected.rds   blood-corrected complete matrix
        ├── DAList_limpa.rds       limpa's complete matrix
        ├── singscore.rds          gene-set scores, subset to the kept samples
        ├── modules.rds            module eigengenes, metadata relabelled
        └── groups.rds             group_of (subject to Hi/Lo), labels, family, traits
```

## Workbook: `clustering.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | One row is | Key columns |
|---|---|---|
| `choice` | the arc | `n_subjects`, `traits`, `k_gap`, `k_imposed`, `imposed_reason`, `k_used`, `grouped` |
| `diagnostics` | a K from 1 to 3 | `logW`, `E.logW`, `gap`, `SE.sim`, `mean_silhouette`, `sizes`, `hr_lr_ari`, `min_jaccard`, `jaccard` |
| `stability` | a cluster at K = 2 or 3 | `k`, `cluster`, `jaccard` |
| `assignments` | a person | the two traits, `k2`, `k3`, `cluster` (empty for anyone not clustered) |
| `space` | a clustered person | `arm` (HR/LR), `group` (Hi/Lo), `cluster`, `PC1`, `PC2` |
| `space_variance` | a component | `variance_explained` |

## Report: `clustering_report.pdf`

1. Gap statistic against K with its SE; the dashed line marks the K the one-SE rule selects.
2. People on the first two components of the scaled traits, coloured by the K = 2 and K = 3
   clusters.
