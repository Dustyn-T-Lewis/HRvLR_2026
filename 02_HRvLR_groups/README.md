# 02 · HRvLR groups

Do the muscle proteome and its response to training and to one bout differ between High and Low
Responders?

## The groups

The study split its sixteen people at the median of a composite hypertrophy score: eight High
Responders (HR) and eight Low Responders (LR). The composite tracks fibre-area change far more
than strength gain (its formula is not in the repository), so HR and LR differ in fibre-area gain
by construction: a fibre-area difference between them restates the split rather than finding
something. HR and LR are different people, and each person gives up to three biopsies (T1, T2,
T3).

## Stages

| Stage | Asks | Main method |
|---|---|---|
| [`01_Design_Contrasts`](01_Design_Contrasts/README.md) | which comparisons | six cell means, blood covariate, person blocked; nine contrasts |
| [`02_Protein_Contrasts`](02_Protein_Contrasts/README.md) | which proteins change, and which differ by group | limpa and limma, logFC with 95% CI; fry concordance; PCA and PERMANOVA |
| [`03_Set_Contrasts`](03_Set_Contrasts/README.md) | which gene sets shift | fgsea (rank against the other proteins) and fry (change across people) |
| [`04_Module_Contrasts`](04_Module_Contrasts/README.md) | which co-expression modules shift | limma on eigengenes, fry on module members |
| [`05_Module_Preservation`](05_Module_Preservation/README.md) | are each group's modules found in the other | WGCNA `modulePreservation` |
| [`06_Classification`](06_Classification/README.md) | can one protein, set or module tell HR from LR | AUC with DeLong and bootstrap CIs, five tasks |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by group | Spearman within each group; limma `feature ~ group * trait` |
| [`Figures`](Figures/README.md) | the arc on a page | F1 overview; F2 and F3 parked; supplements S01 to S03 |

The nine contrasts are defined in `01_Design_Contrasts` and explained in
`02_Protein_Contrasts/README.md`. The interactions (HR's change minus LR's) are the direct test of
a group difference in response; a feature significant in one group and not the other is not, on
that basis, a group difference.

```
02_HRvLR_groups/
├── 01_Design_Contrasts/      design.rds read by every later stage
├── 02_Protein_Contrasts/     fit.rds read by 03 and 04
├── 03_Set_Contrasts/         set_tests.rds read by the figures
├── 04_Module_Contrasts/
├── 05_Module_Preservation/
├── 06_Classification/
│   ├── 01_Proteins/
│   ├── 02_Sets/
│   └── 03_Modules/
├── 07_Trait_Associations/
└── Figures/
    ├── F1_Overview/  F2_Programmes/  F3_People/
    └── S01_Volcano/  S02_NES_Scatter/  S03_Signatures/
```

## Reading the arc

- Seven or eight people per group, so intervals are wide and FDR is a high bar. Read effect sizes
  with their intervals, and nominal results as leads.
- `Baseline_HRvLR` is the floor, not a negative control: the groups come from the outcome, so a
  difference at T1 shows how far the split separates people before training.
- The cluster arcs (`03_Cluster_groups`) rerun these stages on other groupings, and the continuous
  arc (`04_Continuous_no_group`) asks the same questions with traits as numbers.
