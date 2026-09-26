# 05_Volcano

Draws the protein volcanoes with up to eight collapse-surviving fgsea sets ringed.

| | |
|---|---|
| Reads | `02_Contrasts/c_data/set_tests.rds` |
| Writes | `c_data/05_volcano.xlsx`, `b_reports/05_volcano_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/05_Volcano/a_script/05_volcano.R` |
| Cost | about 6 s |

Six FDR panels (both interactions and the four within-arm contrasts), then the two training
contrasts again with pi-ranked labels. The floor is not drawn. Point colour reads protein BH FDR;
the ring reads set fgsea FDR. No protein reaches BH < 0.05, so the FDR panels carry no protein
labels.
