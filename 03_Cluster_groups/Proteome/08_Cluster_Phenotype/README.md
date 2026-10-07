# 08_Cluster_Phenotype

Do people the proteome put in different clusters differ in muscle size or strength?

## The idea

The proteome clusters were formed from module eigengenes alone; no size or strength value played
any part. So a size or strength difference between C1 and C2 is not built in, and these are
tests, not descriptions. This is the arc's direct question: does a split the proteome proposes
line up with how people respond to training?

Every trait is tested against each of the five clusterings (T1, T2, T3, training change, acute
change). Three kinds of trait, from `phenotype_tables()` in `R/helpers.R`:

| Kind | Traits |
|---|---|
| `change` | % change over training in fibre area (mixed, type I, type II), mCSA, leg-press 1RM and leg-extension 1RM; the composite hypertrophy score; total training volume load |
| `T1 level` | the six size and strength traits at T1 |
| `T2 level` | the six traits at T2 |

No trait is measured at T3.

The Wilcoxon rank-sum test compares the two clusters' values by ranks, so it assumes no
distribution and one extreme person cannot drive it. Each cluster's median is the matching
summary.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Clusters | `readxl::read_excel(..., "assignments")` | Each person's cluster per input, from `01_Clustering`. |
| Traits | `phenotype_tables()$pairings`: `training`, `baseline_level`, `post_level` | Long table of subject, kind, trait and value; missing values dropped per trait. |
| Test | `wilcox.test(value ~ cluster)` | Per input, kind and trait: n and median per cluster, W and p. Exact where there are no ties, the normal approximation otherwise. |
| Multiple testing | `p.adjust(method = "BH")` | Within input, over every kind and trait together. |
| HR/LR overlap | `fisher.test(table(arm, cluster))` | Per input: HR and LR counts in each cluster with Fisher's exact p. Description only. |

## How to read the outputs

- `n_1`, `median_1` belong to C1 and `n_2`, `median_2` to C2 for that input. C1 at one input is
  unrelated to C1 at another.
- BH runs within each input, never across the five. Five clusterings tested against the same
  traits give five chances at a hit; read a result that appears at one input only with that in
  mind.
- The clusters are imposed halves. Check each input's Jaccard (`01_Clustering`, `stability`)
  first: a test of an unstable split asks about groups that a resample would redraw.
- At most sixteen people are split in two, so the rank-sum test has little power. A large p is
  not evidence of no difference.
- HR and LR are the top and bottom halves of `comp_hypertrophy` in `phenotype.csv`, so a
  difference in that score partly restates how the clusters cross HR/LR. `hr_lr` shows that
  crossing directly; it describes, and does not test the arc's question.

## Files

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/clustering.xlsx` (`assignments`), `00_Input/phenotype.csv` |
| Writes | `c_data/cluster_phenotype.xlsx`, `b_reports/cluster_phenotype_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/08_Cluster_Phenotype/a_script/cluster_phenotype.R` |

The outputs are not tracked in git.

```
08_Cluster_Phenotype/
├── a_script/cluster_phenotype.R
├── b_reports/cluster_phenotype_report.pdf
└── c_data/cluster_phenotype.xlsx
```

## Workbook: `cluster_phenotype.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | One row is | Key columns |
|---|---|---|
| `tests` | an input, kind and trait | `n_1`, `n_2`, `median_1`, `median_2`, `w`, `p`, `fdr` |
| `hr_lr` | an input | `hr_in_1`, `lr_in_1`, `hr_in_2`, `lr_in_2`, `fisher_p` |

## Report: `cluster_phenotype_report.pdf`

1. Each change trait by cluster, one panel per trait and input: box plot with every person as a
   point. The p for each panel is in `tests`.
