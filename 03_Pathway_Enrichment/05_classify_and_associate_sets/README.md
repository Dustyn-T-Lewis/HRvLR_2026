# 05_classify_and_associate_sets

Tests how well each set score separates the eight tasks and whether it tracks the ten phenotypes.

| | |
|---|---|
| Reads | `00_build_gene_sets/c_data/gene_sets.rds`, `01_run_fgsea_and_fry/c_data/set_tests.rds`, `04_run_singscore/c_data/singscore.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` (sample sheet), `00_Input/phenotype.csv` |
| Writes | `c_data/05_classify_and_associate_sets.xlsx`, `b_reports/05_classify_and_associate_sets_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/05_classify_and_associate_sets/a_script/05_classify_and_associate_sets.R` |
| Cost | about 90 s, most of it drawing |

The figure PDF opens on nominal hits over chance by collection, then an ROC panel for every set
reaching nominal p on each task (floor excluded) and a scatter for every set-outcome pair reaching
nominal p in each window, 12 panels to a page, by collection then p. `set_by_arm` keeps the
within-arm correlations with exact p < 0.05.
