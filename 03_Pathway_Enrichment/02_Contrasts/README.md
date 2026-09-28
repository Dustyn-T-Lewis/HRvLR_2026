# 02_Contrasts

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `00_Gene_Sets/c_data/gene_sets.rds`, `02_Differential_Expression/02_Contrasts/c_data/fit.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/set_tests.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/02_Contrasts/a_script/02_contrasts.R` |
| Cost | about 55 s |

`set_tests.rds` holds `set_tests`, one row per set, contrast and method (`nes`, `leading_edge` and
`main` are fgsea-only), and `protein_results`, the nine contrasts rebuilt from the fit with
`topTable()`. The workbook keeps the sets at FDR < 0.05. The figure PDF holds one dot plot per contrast of every
fgsea collapse survivor at FDR < 0.05, placed by NES, then one per contrast of every set at fry
p < 0.05, placed by signed −log10 p. Forty sets to a page, coloured by collection.
