# 06_Classify/01_Proteins

A copy of `02_HRvLR_groups/05_Classify/01_Proteins` run on the strength clusters: groups Hi and Lo from `01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

Scores how well each protein separates the groups on eight tasks.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_corrected.rds`, `DAList_limpa.rds` (overview row order) |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classify/01_Proteins/a_script/classify_proteins.R` |

## Method

Eight tasks on limpa's complete blood-corrected log2 abundance; `observed_share` marks how much of
each protein was measured rather than estimated. Within-arm tasks (training and the acute bout, in each arm) pair each
subject's two biopsies: the Wilcoxon signed-rank test gives p and the effect is the median paired
change. Hi-vs-Lo tasks (T1 level, T2 level, training change, acute change) use the rank-sum test,
and AUC = W / (n1 x n2), above 0.5 when higher in Hi. A feature needs three observations per side.
BH within task. `counts` gives each task's smallest attainable p: 2 / 2^n for n pairs, 2 / choose(n1
+ n2, n1) for two groups.

## Figures

An overview page: median paired change and AUC for every protein, each on its own colour scale,
beside a status map. Then every protein at p < 0.05 in every task, twelve panels to a page, ordered
by p, then the larger effect, then ID: a paired grid for within-arm tasks, an ROC grid
(`pROC::ggroc`) for Hi-vs-Lo tasks. Headers give the effect, p, BH q and the protein's own n.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per task: pairing, n per side, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per protein and task: AUC with a DeLong 95% interval (between-arm tasks, none at AUC 0 or 1), median paired change (paired tasks), Wilcoxon p, BH FDR and Storey q. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every p < 0.05 row. |
