# 03_Classify

Scores how well each set's singscore separates the groups on the protein level's eight tasks.

| | |
|---|---|
| Reads | `01_Scores/c_data/singscore.rds`, `00_Gene_Sets/c_data/gene_sets.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` (sample sheet) |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/03_Classify/a_script/03_classify.R` |
| Cost | about 15 s |

The tasks, AUC and tests are those of `02_Differential_Expression/03_Classify`. BH and the chance
count run within task and collection. `chance_expectation` gives each task and collection's group
sizes, smallest attainable p (HR training 0.031, HR acute 0.016, LR 0.0078), tested, nominal and BH
counts. With 6 to 8 pairs fewer than 5% of null tests reach 0.05, so a ratio below 1 is not worse
than chance.

The figure PDF opens on the chance ratio per task and collection, then draws an ROC curve
(`pROC::ggroc`) for every set at p < 0.05 on each task, by p, twelve to a page. The floor task gets
no curves.
