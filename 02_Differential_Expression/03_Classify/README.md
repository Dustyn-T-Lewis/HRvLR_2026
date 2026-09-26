# 03_Classify

Scores how well each protein separates the groups on eight tasks.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds` |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 02_Differential_Expression/03_Classify/a_script/03_classify.R` |
| Cost | about 15 s |

The tasks are training (T1 to T2) and the acute bout (T2 to T3) within each arm, paired, and HR
against LR on the T1 level, the T2 level, the training change and the acute change. AUC is the
rank-sum statistic over n1·n2, so above 0.5 means higher at the later timepoint or in HR. p is
Wilcoxon, signed-rank when paired; BH within task. A protein needs three observations per side.

`chance_expectation` gives each task's group sizes, smallest attainable p, tested, nominal and BH
counts, and the chance expectation (5% of tested). HR has 6 training and 7 acute pairs and LR 8 of
each, so HR's smallest p is 0.031 (training) and 0.016 (acute) against LR's 0.0078. With so few
pairs fewer than 5% of null tests reach 0.05, and a ratio below 1 is not worse than chance.

The figure PDF opens on the chance ratio per task, then draws an ROC curve (`pROC::ggroc`) for
every protein at p < 0.05 on each task, by p, twelve to a page. The floor task gets no curves. Each
page caption repeats the task's counts and smallest attainable p.
