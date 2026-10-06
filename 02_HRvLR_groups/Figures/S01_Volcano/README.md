# S01_Volcano

Draws the protein volcanoes with up to eight collapse-surviving fgsea sets ringed.

| | |
|---|---|
| Reads | `02_Contrasts/c_data/set_tests.rds` |
| Writes | `c_data/s01_volcano.xlsx`, `b_reports/s01_volcano_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/Figures/S01_Volcano/a_script/s01_volcano.R` |

## Method

`enrichVolcano::volcano_ring()` on the protein results from `set_tests.rds`. Point colour reads
protein BH FDR; the ring reads set fgsea FDR. The floor contrast is not drawn.

## Figures

Six volcanoes by FDR (both interactions and the four within-arm contrasts), then the two training
contrasts again with pi-ranked labels.

## Workbook

| Sheet | Holds |
|---|---|
| `ringed_sets` | The fgsea sets ringed on each volcano. Points come from set_tests.rds, protein_results. |
| `labelled_proteins` | The proteins named on each volcano, by FDR or by pi score. |
