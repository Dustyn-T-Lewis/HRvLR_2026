# Muscle_strength

Do people grouped by how much stronger they became differ in their proteome and its response?

People are clustered on their % change over training in leg-press and leg-extension 1RM
(`pct_leg_press_1rm`, `pct_leg_ext_1rm`); anyone missing one has no group. No proteome value enters
the clustering. The cluster with the higher mean `pct_leg_press_1rm` is `Hi`.

The gap statistic selects K = 1; K = 2 is imposed because both clusters are stable by bootstrap
Jaccard and the split is unrelated to HR/LR (adjusted Rand index -0.07). Results here are
exploratory.

## Steps

| Step | Asks |
|---|---|
| [`01_Clustering`](01_Clustering/README.md) | how people split on strength change; writes `groups.rds` |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | which comparisons |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | which proteins differ |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | which gene sets shift |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | which modules shift |
| [`06_Classification`](06_Classification/README.md) | can one protein, set or module tell Hi from Lo |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by group |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the folder on a page |

Steps 02 to 07 call the same functions as `03_Group_comparisons/HRvLR`, with Hi and Lo in place of HR and LR.

## Read with care

- Hi and Lo were defined by strength change, so a strength-change difference between them is built in.
- K = 2 over a gap-statistic K = 1 cuts one spread of people in two; it does not show two kinds of
  people.
