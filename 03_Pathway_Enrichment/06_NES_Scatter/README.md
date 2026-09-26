# 06_NES_Scatter

Plots each set's HR NES against its LR NES, once for training and once for the acute bout.

| | |
|---|---|
| Reads | `02_Contrasts/c_data/set_tests.rds` |
| Writes | `c_data/06_nes_scatter.xlsx`, `b_reports/06_nes_scatter_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/06_NES_Scatter/a_script/06_nes_scatter.R` |
| Cost | about 9 s |

Four composites: for each pair, every collection with collapse survivors and discordant sets, then
Hallmark and GO Slim with each concordant quadrant rescaled so every significant set is named.
