# 05_Classify/03_Modules

Scores how well each module eigengene separates HR from LR on five between-group tasks.

| | |
|---|---|
| Reads | `01_Trunk/06_Modules/c_data/modules.rds` |
| Writes | `c_data/03_modules.xlsx`, `b_reports/03_modules_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/05_Classify/03_Modules/a_script/classify_modules.R` |

## Method

`classify_step()` in `R/classify.R`, on module eigengenes. The five tasks compare the groups: the T1, T2
and T3 levels, and the training (T2 - T1) and acute (T3 - T2) changes. Within-group change is the
subject of `02_DEP`, so no task here pairs a person with themselves.

Nothing is fitted. Each module's AUC is the Mann-Whitney statistic W / (n1 x n2), above 0.5
when higher in HR, with the Wilcoxon rank-sum p, BH and Storey q within task. Every AUC
carries DeLong's 95% interval (`pROC::ci.auc`); with seven or eight people per group it runs
narrow and should be read as approximate. The modules drawn in the report also carry pROC's
stratified bootstrap interval (2,000 resamples, seeded). `counts` gives each task's smallest
attainable p, 2 / choose(n1 + n2, n1). The share of modules at p < 0.05 is a multiple-testing
summary to set against the 5% expected by chance, not a measure of how well the groups separate.

## Figures

An overview page: AUC for every module in the five tasks beside a status map. Then ROC curves
(`pROC::ggroc`) for the modules at p < 0.05 in each task, at most 24 per task, ordered by p,
then distance of AUC from 0.5. Headers give the AUC with its bootstrap interval, p, BH q and the
module's n per group. Every row stays in the workbook.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per module and task: AUC with DeLong 95% interval (`auc_lower`, `auc_upper`), Wilcoxon p, BH FDR and Storey q; bootstrap interval (`boot_lower`, `boot_upper`) for the modules drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position of every drawn row. |
