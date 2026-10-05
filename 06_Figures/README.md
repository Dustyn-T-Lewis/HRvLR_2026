# 06 · Figures

Builds figures from saved outputs; nothing here refits a model. Each figure is one substage with
one script, its PDF and PNG in `b_reports/` and the plotted numbers in `c_data/`. Each script checks
its plotted numbers against the source tables before it writes.

| Figure | Draws |
|---|---|
| [`F01_Overview`](F01_Overview/README.md) | design, proteins per contrast, proteome and phenotype PCA, shared pathways, training and acute signature ranks |
| [`F02_Groups`](F02_Groups/README.md) | interaction pathways, OXPHOS catch-up and acute response, sarcomere-strength slope by arm, classification null |
| [`F03_Continuum`](F03_Continuum/README.md) | screens against chance, sarcomere growth x strength interaction, its drivers and pathways, baseline strength proteins, clustering |
| [`S01_Volcano`](S01_Volcano/README.md) | protein volcano inside a ring of enriched sets, per contrast |
| [`S02_NES_Scatter`](S02_NES_Scatter/README.md) | NES of one contrast against another |
| [`S03_Signatures`](S03_Signatures/README.md) | signature ranks with the Π sets; training sets followed through the bout |
| [`S04_QC`](S04_QC/README.md) | filtering, detection, intensity, dynamic range, missingness, sample correlation, blood, p-value calibration and power |

Main figures are `F0n`, supplementary `S0n`. Two shared files serve every figure:
`R/figure_theme.R` holds the contrast palette (sampled from the F01 design panel), the time and
direction colours, the base theme, panel-letter placement and coloured facet strips;
`R/figure_panels.R` holds the panel builders (PCA, signature barcode, pathway pruning and labels)
and `save_figure()`, which writes a cairo PDF with embedded fonts and a 600 dpi PNG.
