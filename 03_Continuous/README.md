# 03_Continuous

Which proteins, gene sets and modules follow muscle size and strength when all 16 people are
pooled and each trait stays a number? No group label enters any model. Traits come from
`00_Setup/Input/phenotype.csv` through `trait_windows()` in `02_HRvLR/helpers.R`: fibre area
(mixed, type I, type II), muscle area and two 1RMs as levels at T1 and T2 and as % change, plus the
composite hypertrophy score and total volume load. T3 has no phenotype.

| Step | Does | Writes |
|---|---|---|
| [`01_Trait_Associations`](01_Trait_Associations/README.md) | sample model, change score, two-predictor and interaction fits for three feature types | `associations.rds`, `trait_associations.xlsx` |
| [`02_Trait_Pathways`](02_Trait_Pathways/README.md) | fgsea on each window-and-trait protein ranking, FDR within database; Jaccard tree selection for the figure | `trait_pathways.rds`, `trait_pathways.xlsx` |
| [`03_Profile_Clusters`](03_Profile_Clusters/README.md) | proteins grouped by their pattern of trait links, and per trait by co-variation | `profile_clusters.rds`, `profile_clusters.xlsx` |
| [`04_Elastic_Net`](04_Elastic_Net/README.md) | elastic net of traits on proteins, modules and Hallmark scores in six pairings: stability selection and held-out Q² | `elastic_net.rds`, `elastic_net.xlsx` |

Six windows pair a feature value with its traits, 42 window-trait screens in all:

| window | feature value | traits |
|---|---|---|
| `baseline_level` | T1 level | trait levels at T1 |
| `baseline` | T1 level | % change, composite score, volume load |
| `training` | change T1 → T2 | % change, composite score, volume load |
| `post_level` | T2 level | trait levels at T2 |
| `acute` | change T2 → T3 | % change, composite score, volume load |
| `post_bout_level` | T3 level | trait levels at T2 |

## Inputs

`01_Processing` outputs (`proteins_corrected.rds`, `singscore.rds`, `gene_sets.rds`,
`modules.rds`), `00_Setup/Input/phenotype.csv` and `02_HRvLR/helpers.R`.

## Run

```sh
Rscript run_all.R 03_Continuous
```

`02` and `03` read the protein change scores `01` saves, so `01` runs first.

## Read with care

Every test rests on 16 people or fewer, and HR and LR are the extremes of 25, so correlations
across these 16 run stronger than they would across all 25. Protein-level calls at FDR < 0.05 are
rare (3 in the change score); set-level calls are many, and fgsea treats proteins as independent.
