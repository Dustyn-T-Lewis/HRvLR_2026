# 02_enrich_volcano_fgsea

Draws the protein volcanoes with up to eight collapse-surviving fgsea sets ringed.

| | |
|---|---|
| Reads | `01_run_fgsea_and_fry/c_data/set_tests.rds` |
| Writes | `c_data/02_enrich_volcano_fgsea.xlsx`, `b_reports/02_enrich_volcano_fgsea_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.R` |
| Cost | about 6 s |

Six FDR panels (both interactions and the four within-arm contrasts), then the two training
contrasts again with pi-ranked labels. The floor is not drawn. Point colour reads protein BH FDR;
the ring reads set fgsea FDR. No protein reaches BH < 0.05, so the FDR panels carry no protein
labels.
