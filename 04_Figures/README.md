# 04_Figures

Main figures F1 to F5 and supplements S1 to S5, drawn from the saved results of `01_Processing`,
`02_HRvLR` and `03_Continuous`. No figure fits a model; each reads tables and draws them.

The analysis steps are exploratory, and their results informed which panels became main figures
and which went to supplements (see the root README).

| Figure | Shows | Supplement |
|---|---|---|
| F1 | design; PCA per time point with PERMANOVA, gap k and forced-split ARI; trait change per group; proteins and gene sets per contrast; each group's training and acute responders tested in the other group (fry) | S1: no natural clusters (correlation heatmap, gap statistic); variance shares; DiffVar; PCs against covariates |
| F2 | proteins each group changes with training or the bout | S2: single-feature AUC; exploratory HR − LR protein candidates |
| F3 | pathways: shared, HR-only and LR-only blocks for training and the bout, then between groups; one set per tree cluster | S3: every candidate set grouped by tree, cuts 0.75 and 0.9, GO Slim overview |
| F4 | traits: people's changes; links against chance; trait pathways | S4: HR slope minus LR slope; protein profile groups |
| F5 | modules: score heatmap (contrasts and trait links), coupling after the bout within each group, preservation | S5: module characterisation with STRING; WGCNA module QC; the smallest-p module person by person. Modules/: per-module hubs and rewiring; coupling, per-person responses, PLS |
| | | S6: elastic net of traits, everyone pooled: every fit, what carries the lead, held-out predictions, HR against LR |

## Layout

```
style/figure_style.R            palette, theme, contrast and trait labels, strips (every panel)
style/pathway_grid.R            the pathway dot-matrix (F3, S3, F4) and its one-set-per-cluster rows
style/protein_grid.R            the protein dot-matrix (F2)
style/share_grid.R              the share-against-chance tiles of the trait screens (F4b, S4a)
style/module_labels.R           module names from their leading set (F5, S2, S5, Modules)
style/score_heatmap.R           the YvO-style score heatmap (F5a, S3e)
style/elastic_net.R             the elastic-net results and the lead S6 draws
style/composite.R               run a figure's panels, tag them, write its workbook
assets/design.png               study design schematic (F1a)
F1/a_script/panels/a_*.R        one script per panel
F1/a_script/F1_composite.R      sources the panels and stitches them
F1/b_reports/panels/a_*.pdf     one PDF per panel
F1/b_reports/F1_composite.pdf   the figure
F1/c_data/F1.xlsx               one sheet per panel, the numbers drawn
Modules/a_script/modules_pages.R  every module page, from 02_HRvLR/08_Module_Detail
Modules/b_reports/<module>/     hubs.pdf and rewiring.pdf per module
Modules/b_reports/overview/     coupling.pdf, per_person.pdf, pls.pdf
```

Each panel script runs on its own from the project root and leaves four objects: `panel` (the
plot), `panel_data` (the plotted numbers, checked against the source table), `holds` (one line on
what the sheet holds) and `inputs` (the files read). The composite runs every panel script with
`run_panels()`, stitches the plots with patchwork and writes the workbook with
`write_figure_workbook()`.

## Run

```sh
Rscript run_all.R 04_Figures
Rscript 04_Figures/F1/a_script/panels/b_pca.R   # one panel
```
