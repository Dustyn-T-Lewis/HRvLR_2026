# 04_Figures

Main figures F1 to F5, each with its supplements, drawn from the saved results of `01_Processing`,
`02_HRvLR` and `03_Continuous`. No figure fits a model; each reads tables and draws them.

The analysis steps are exploratory, and their results informed which panels became main figures
and which went to supplements (see the root README).

| Folder | Main figure | Supplements (in `supp/`) |
|---|---|---|
| `F1/` | design; PCA per time point with PERMANOVA, gap k and forced-split ARI; trait change per group; proteins and gene sets per contrast; each group's training and acute responders tested in the other group (fry) | S1: no natural clusters (correlation heatmap, gap statistic); variance shares; DiffVar; PCs against covariates |
| `F2/` | proteins each group changes with training or the bout | S2: single-feature AUC; exploratory HR − LR protein candidates |
| `F3/` | pathways: shared, HR-only and LR-only blocks for training and the bout, then between groups; one set per tree cluster | S3: every candidate set grouped by tree, cuts 0.75 and 0.9, GO Slim themes, Hallmark score heatmap |
| `F4/` | traits: people's changes; links against chance; trait pathways | S4: HR slope minus LR slope; protein profile groups. S6: elastic net of traits |
| `F5/` | modules: score heatmap (contrasts and trait links), coupling after the bout within each group, preservation | S5: module characterisation with STRING; WGCNA module QC; the smallest-p module person by person. S7: hubs and rewiring per module, coupling, per-person responses, PLS |

## Layout

Every figure folder has the same shape:

```
F4/
├── README.md               what each panel shows, its inputs, how to run it
├── a_script/
│   ├── a_trait_changes.R   one script per panel
│   ├── ...
│   ├── F4_composite.R      runs the panels, stitches the figure, writes the workbook
│   └── supp/               S4_a_slope_counts.R ... S4_composite.R, S6_a_fit_overview.R ... S6_composite.R
├── b_reports/
│   ├── a_trait_changes.pdf ...  F4_composite.pdf
│   └── supp/               S4_*.pdf, S4_composite.pdf, S6_*.pdf, S6_composite.pdf
└── c_data/                 F4.xlsx, S4.xlsx, S6.xlsx: one sheet per panel, the numbers drawn
```

Shared builders:

```
style/figure_style.R    palette, theme, labels, strips, figure_dir() and save_panel()
style/composite.R       run a figure's panels, tag them, save the figure and its workbook
style/pathway_grid.R    the pathway dot-matrix (F3, S3, F4)
style/protein_grid.R    the protein dot-matrix (F2, S2)
style/share_grid.R      the share-against-chance tiles of the trait screens (F4, S4)
style/pca_grid.R        the PCA rows (F1, S1)
style/module_labels.R   module names from their leading set (F5, S2, S5, S7)
style/score_heatmap.R   the YvO-style score heatmap (F5, S3)
style/elastic_net.R     the elastic-net results and the lead S6 draws
assets/design.png       study design schematic (F1 a)
```

Each panel script runs on its own from the project root and leaves `panel` (the plot) and
`panel_data` (the plotted numbers, checked against the source table). The composite runs every
panel script with `run_panels()`, stitches the plots with patchwork and writes the workbook with
`write_figure_workbook()`. S7 draws the step tables of `02_HRvLR/08_Module_Detail` whole, so it
keeps no workbook of its own.

## Run

```sh
Rscript run_all.R 04_Figures                  # every figure and supplement
Rscript run_all.R 04_Figures/F4               # F4, S4 and S6
Rscript 04_Figures/F4/a_script/b_links_vs_chance.R   # one panel
```
