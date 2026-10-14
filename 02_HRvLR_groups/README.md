# 02_HRvLR_groups

Do the muscle proteome and its response to training and to one bout differ between High and Low
Responders?

HR and LR are the top and bottom eight of 25 trained people on a composite hypertrophy score,
(2 × %ΔmCSA + %ΔfCSA) / 3 (`02_Groups/HRvLR`). Growth defines the groups, so a growth difference
between them restates the selection rather than finding something.

| Step | Asks | Method |
|---|---|---|
| [`01_Design_Contrasts`](01_Design_Contrasts/README.md) | which comparisons | six cell means, blood covariate, person blocked; nine contrasts |
| [`02_Protein_Contrasts`](02_Protein_Contrasts/README.md) | which proteins change, and which differ by group | limpa and limma, logFC with 95% CI; fry concordance; PCA and PERMANOVA |
| [`03_Set_Contrasts`](03_Set_Contrasts/README.md) | which gene sets shift | fgsea and fry |
| [`04_Module_Contrasts`](04_Module_Contrasts/README.md) | which modules shift | limma on eigengenes, fry on members |
| [`05_Module_Preservation`](05_Module_Preservation/README.md) | is each group's module structure found in the other | WGCNA `modulePreservation` |
| [`06_Classification`](06_Classification/README.md) | can one protein, set or module tell HR from LR | AUC with DeLong and bootstrap CIs, five tasks |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by group | Spearman within group; limma `feature ~ group * trait` |
| [`Figures`](Figures/README.md) | the folder on a page | F1; F2 and F3 parked; supplements S01 to S03 |

Steps 01 to 04, 06 and 07 are shared functions in `R/` (named in each README), called here with
HR and LR and in `03_Cluster_groups` with the cluster groups. `05_Module_Preservation` runs only
here.

## Read with care

- Seven or eight people per group: intervals are wide and FDR is a high bar. Read effect sizes
  with their intervals, and nominal results as leads.
- The interactions (HR's change minus LR's) are the test of a group difference in response. A
  feature significant in one group and not the other is not, on that basis, a group difference.
- `Baseline_HRvLR` is the floor, not a negative control.
