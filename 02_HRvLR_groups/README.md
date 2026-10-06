# 02 · HRvLR groups

![HR/LR groups: who, how they are grouped, what we ask](../docs/concepts/hrvlr_groups.png)

The study's own split: High and Low Responders, the median split of the composite hypertrophy
score, eight subjects each. The composite tracks fibre-area change (Spearman 0.90) far more than
strength gain (0.11), so a fibre-area difference between the groups is built in.

| Substage | Asks | Method |
|---|---|---|
| [`01_Design`](01_Design/README.md) | which comparisons | six cell means, blood covariate, subject blocked; nine contrasts |
| [`02_DEP`](02_DEP/README.md) | which proteins change, per contrast | limpa `dpcDE`, limma `eBayes(robust)`, logFC with 95% CI; PCA and PERMANOVA; fry concordance |
| [`03_Pathways`](03_Pathways/README.md) | which gene sets shift, per contrast | fgsea and fry |
| [`04_Modules/01_Contrasts`](04_Modules/01_Contrasts/README.md) | which modules shift, per contrast | eigengene limma, fry on members |
| [`04_Modules/02_Preserve`](04_Modules/02_Preserve/README.md) | are each group's modules found in the other | WGCNA `modulePreservation` |
| [`05_Classify`](05_Classify) | does one feature separate the groups | Wilcoxon AUC with DeLong and bootstrap CIs, at protein, set and module level |
| [`06_Associate`](06_Associate/README.md) | does a phenotype link differ by group | Spearman within group; limma `feature ~ group * trait` |
| [`Figures/`](Figures) | Figures 1 to 3 and supplements S01 to S03 | saved outputs only |

The interaction contrasts are the direct test of a group difference in response. A feature
significant in one group and not the other is not, on that basis, a group difference.
