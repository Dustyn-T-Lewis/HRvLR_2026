# 06_Classification/02_Sets

Scores how well each gene set's singscore separates HR from LR on five between-group tasks.

| | |
|---|---|
| Reads | `01_Trunk/05_Set_Scores/c_data/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/03_Quantification/c_data/DAList_corrected.rds` (sample sheet) |
| Writes | `c_data/02_sets.xlsx`, `b_reports/02_sets_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/02_Sets/a_script/classify_sets.R` |

## Method

`classify_step()` in `R/classify.R`, on singscore. The five tasks compare the groups: the T1, T2
and T3 levels, and the training (T2 - T1) and acute (T3 - T2) changes. Within-group change is the
subject of `02_Protein_Contrasts`, so no task here pairs a person with themselves.

Nothing is fitted. Each set's AUC is the Mann-Whitney statistic W / (n1 x n2), above 0.5
when higher in HR, with the Wilcoxon rank-sum p, BH and Storey q within task and collection. Every AUC
carries DeLong's 95% interval (`pROC::ci.auc`); with seven or eight people per group it runs
narrow and should be read as approximate. The sets drawn in the report also carry pROC's
stratified bootstrap interval (2,000 resamples, seeded). `counts` gives each task's smallest
attainable p, 2 / choose(n1 + n2, n1). The share of sets at p < 0.05 is a multiple-testing
summary to set against the 5% expected by chance, not a measure of how well the groups separate.

## Figures

One overview page per collection: AUC for every set in the five tasks beside a status map. Then ROC curves
(`pROC::ggroc`) for the sets at p < 0.05 in each task, at most 24 per task, ordered by p,
then distance of AUC from 0.5. Headers give the AUC with its bootstrap interval, p, BH q and the
set's n per group. Every row stays in the workbook.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts. |
| `auc` | Per set and task: AUC with DeLong 95% interval (`auc_lower`, `auc_upper`), Wilcoxon p, BH FDR and Storey q; bootstrap interval (`boot_lower`, `boot_upper`) for the sets drawn in the report. |
| `plot_index` | PDF page of each overview, and page and panel position of every drawn row. |
