# 03_enrich_scatter_fgsea

Plots each set's HR NES against its LR NES, once for training and once for the acute bout.

| | |
|---|---|
| Reads | `01_run_fgsea_and_fry/c_data/set_tests.rds` |
| Writes | `c_data/03_enrich_scatter_fgsea.xlsx`, `b_reports/03_enrich_scatter_fgsea_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/03_enrich_scatter_fgsea/a_script/03_enrich_scatter_fgsea.R` |
| Cost | about 9 s |

Four composites: for each pair, every collection with collapse survivors and discordant sets, then
Hallmark and GO Slim with each concordant quadrant rescaled so every significant set is named.
