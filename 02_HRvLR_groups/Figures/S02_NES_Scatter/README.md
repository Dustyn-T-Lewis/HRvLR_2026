# S02_NES_Scatter

Plots each set's HR NES against its LR NES, for training and for the acute bout.

| | |
|---|---|
| Reads | `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` |
| Writes | `c_data/s02_nes_scatter.xlsx`, `b_reports/s02_nes_scatter_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/Figures/S02_NES_Scatter/a_script/s02_nes_scatter.R` |

## Method

A set is discordant when its two NES differ in sign. Spearman rho of HR against LR NES summarises
agreement per population: all sets, collapse survivors, and Hallmark with GO Slim.

## Figures

Four composites: per pair, every collection with collapse survivors and discordant sets, then
Hallmark and GO Slim with each concordant quadrant rescaled so every significant set is named.

## Workbook

| Sheet | Holds |
|---|---|
| `concordance` | Spearman rho of HR against LR NES, and significant and discordant counts, per population. |
| `nes_scatter` | One row per set scored in both arms: both NES, both FDR, significance and discordance. |
