# 03_Classify

Scores how well each module eigengene separates the groups on the protein level's eight tasks.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 04_Network/03_Classify/a_script/03_classify.R` |

## Method

Eight tasks on eigengenes. Within-arm tasks (training and the acute bout, in each arm) pair each
subject's two biopsies: the Wilcoxon signed-rank test gives p and the effect is the median paired
change. HR-vs-LR tasks (T1 level, T2 level, training change, acute change) use the rank-sum test,
and AUC = W / (n1 x n2), above 0.5 when higher in HR. A feature needs three observations per side.
BH within task. `counts` gives each task's smallest attainable p: 2 / 2^n for n pairs, 2 / choose(n1
+ n2, n1) for two groups.

## Figures

An overview page: median paired change and AUC for every module, each on its own colour scale,
beside a status map. Then every module at p < 0.05 in every task, twelve panels to a page, ordered
by p, then the larger effect, then ID: a paired grid for within-arm tasks, an ROC grid
(`pROC::ggroc`) for HR-vs-LR tasks. Headers give the effect, p, BH q and the module's own n.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per task: pairing, n per side, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per module and task: AUC with a DeLong 95% interval (between-arm tasks, none at AUC 0 or 1), median paired change (paired tasks), Wilcoxon p, BH FDR and Storey q. |
| `plot_index` | PDF page of each overview, and page and panel position (1 to 12) of every p < 0.05 row. |
