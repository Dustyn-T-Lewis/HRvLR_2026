# 04_run_singscore

Scores every sample on every set with singscore.

| | |
|---|---|
| Reads | `00_build_gene_sets/c_data/gene_sets.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/singscore.rds`, `c_data/04_run_singscore.xlsx`, `b_reports/04_run_singscore_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.R` |
| Cost | about 9 s |

`singscore.rds` holds `scores`, 1,378 sets by 45 samples. The figure PDF holds mean score against
mean dispersion per set, and the score distribution by group.
