# 01_Scores

Scores every sample on every set with singscore.

| | |
|---|---|
| Reads | `00_Gene_Sets/c_data/gene_sets.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/singscore.rds`, `c_data/01_scores.xlsx`, `b_reports/01_scores_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/01_Scores/a_script/01_scores.R` |
| Cost | about 9 s |

`singscore.rds` holds `scores`, 1,378 sets by 45 samples; `03_Classify` and `04_Associate` read it. The figure PDF holds mean score against
mean dispersion per set, and the score distribution by group.
