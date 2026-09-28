# 03_Classify

Scores how well each module eigengene separates the groups on the protein level's eight tasks.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 04_Network/03_Classify/a_script/03_classify.R` |
| Cost | about 5 s |

Within-arm tasks (training and the acute bout in each arm) pair each subject's two biopsies. Their
p is the Wilcoxon signed-rank test and their effect the median paired change. Between-arm tasks
(HR against LR on the T1 level, the T2 level, the training change and the acute change) use the
rank-sum test, and AUC is W / (n1·n2), above 0.5 when higher in HR. BH runs within task. A
module needs three observations per side.

`counts` gives each task's pairing, n per side, smallest attainable p (2 / 2^n for n pairs,
2 / choose(n1 + n2, n1) for two groups), and the tested, p < 0.05 and BH < 0.05 counts. HR has 6
training and 7 acute pairs, LR 8 of each, so HR's smallest p is 0.031 and 0.016 against LR's
0.0078. `auc` holds every module and task; its AUC for paired tasks ignores the pairing and is not
drawn.

The figure PDF opens with an overview page: every module's median paired change in the four
within-arm tasks and AUC in the four HR-vs-LR tasks, each on its own colour scale, beside a status
map (not significant, p < 0.05, BH < 0.05, not tested). Rows follow module size. The grids follow:
every module at p < 0.05 in every task, baseline included, twelve panels to a page, ordered within
task by p, then the larger effect, then ID. Within-arm tasks get a paired grid, each subject's two
biopsies joined by a line; between-arm tasks get an ROC grid (`pROC::ggroc`). Each header gives the
module, its effect, p, BH q and its own n (complete pairs, or HR and LR subjects). Each page
subtitle gives the task's counts and each caption the test, the largest n and the smallest
attainable p. `plot_index` gives the PDF page and panel of every drawn row.
