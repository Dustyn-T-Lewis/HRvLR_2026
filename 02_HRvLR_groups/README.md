# 02 · HRvLR groups

The study's own split: High and Low Responders, the median split of the composite hypertrophy
score, eight subjects each. The composite tracks fibre-area change (Spearman 0.90) far more than
strength gain (0.11), so a fibre-area difference between the groups is built in.

| Substage | Asks | Method |
|---|---|---|
| [`01_Design`](01_Design/README.md) | which comparisons | six cell means, blood covariate, subject blocked; nine contrasts |
| [`02_DEP`](02_DEP/README.md) | which proteins change, per contrast | limpa `dpcDE`, limma `eBayes(robust)`, `treat`; ANCOVA; PCA and PERMANOVA; fry concordance |
| [`03_Pathways`](03_Pathways/README.md) | which gene sets shift, per contrast | fgsea, cameraPR, fry |
| [`04_Modules`](04_Modules/README.md) | which modules shift, per contrast | eigengene limma, fry on members, ANCOVA |
| [`05_Classify`](05_Classify) | does one feature separate the groups | Wilcoxon and AUC, at protein, set and module level |
| [`06_Associate`](06_Associate/README.md) | does a phenotype link differ by group | Spearman within group; `feature ~ group * trait`; DGCA |
| [`07_Preserve`](07_Preserve/README.md) | do HR and LR share a network | WGCNA `modulePreservation` against random halves |
| [`Figures/`](Figures) | Figure 1 and supplements S01 to S03 | saved outputs only |

The interaction contrasts are the direct test of a group difference in response. A feature
significant in one group and not the other is not, on that basis, a group difference.
