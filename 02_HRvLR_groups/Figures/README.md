# Figures

The HR/LR arc's main figures and supplements, drawn from the saved outputs of stages 01 to 07.

## The figures

| Folder | Status | Asks | Draws from |
|---|---|---|---|
| [`F1_Overview`](F1_Overview/README.md) | active | What does the design look like, and how does the whole proteome respond before any single result? | `02_Protein_Contrasts`, `03_Set_Contrasts`, `phenotype.csv` |
| [`F2_Programmes`](F2_Programmes/README.md) | parked draft | Which co-expression programmes move with training and the bout, and do the groups differ in them? | `01_Design_Contrasts`, `03_Set_Contrasts`, `04_Module_Contrasts`, trunk modules |
| [`F3_People`](F3_People/README.md) | parked draft | Can one feature tell the groups apart, and does a protein's link to strength gain differ between them? | `06_Classification`, `07_Trait_Associations`, the continuous arc |
| [`S01_Volcano`](S01_Volcano/README.md) | active | Which proteins and which fgsea sets move in each within-group and interaction contrast? | `03_Set_Contrasts` |
| [`S02_NES_Scatter`](S02_NES_Scatter/README.md) | active | Do HR and LR shift the same gene sets, set by set? | `03_Set_Contrasts` |
| [`S03_Signatures`](S03_Signatures/README.md) | active | Supplement to F1: the signature ranks with the stricter Π sets, and the training proteins followed through the bout. | `02_Protein_Contrasts`, `03_Set_Contrasts` |

Parked means the script still runs (`run_all.R` includes it) and its outputs are rebuilt, but the
figure is not final and its README describes what it draws, not what it shows.

## What they share

- They read saved outputs only. Each script lists its `inputs <- c(...)`; any computation it adds
  is a description with no new test (F1's phenotype PCA, F3's Spearman across proteins, S02's
  Spearman of NES across sets).
- Each writes `b_reports/<stem>.pdf` and `c_data/<stem>.xlsx`. The workbook holds the numbers each
  panel draws; F2 and F3 add `panel_sources`, naming the file, sheet and selection rule behind
  every panel.
- F1, F2, F3 and S03 check the plotted numbers against their source tables before writing
  anything, and stop on a mismatch.
- Colours and theme come from `R/figure_theme.R`. The contrast palette is sampled from F1's design
  panel: Training blue, Acute orange, HR - LR slate, Interaction yellow. Within Training and Acute
  the first group (HR) takes a darker shade of the family colour. Time points take the colour of
  the window that produced them: T1 grey, T2 blue, T3 orange.
- Panel builders live in `R/figure_panels.R` and, for F1, `R/figure_overview.R`. Both are written
  for any two groups; the cluster arcs' F1 scripts call the same builders with their own group
  names.
- PDFs go through `grDevices::cairo_pdf`, which embeds fonts and keeps the Π glyph the base PDF
  device drops.

## Run

After every analysis stage, in the order `run_all.R` uses (F3 reads the continuous arc, so the
figures run last):

```
Rscript 02_HRvLR_groups/Figures/F1_Overview/a_script/f1_overview.R
Rscript 02_HRvLR_groups/Figures/F2_Programmes/a_script/f2_programmes.R
Rscript 02_HRvLR_groups/Figures/F3_People/a_script/f3_people.R
Rscript 02_HRvLR_groups/Figures/S01_Volcano/a_script/s01_volcano.R
Rscript 02_HRvLR_groups/Figures/S02_NES_Scatter/a_script/s02_nes_scatter.R
Rscript 02_HRvLR_groups/Figures/S03_Signatures/a_script/s03_signatures.R
```

```
Figures/
├── F1_Overview/      a_script/, assets/design.png, b_reports/f1_overview.pdf, c_data/f1_overview.xlsx
├── F2_Programmes/    a_script/, b_reports/f2_programmes.pdf, c_data/f2_programmes.xlsx
├── F3_People/        a_script/, b_reports/f3_people.pdf, c_data/f3_people.xlsx
├── S01_Volcano/      a_script/, b_reports/s01_volcano.pdf, c_data/s01_volcano.xlsx
├── S02_NES_Scatter/  a_script/, b_reports/s02_nes_scatter.pdf, c_data/s02_nes_scatter.xlsx
└── S03_Signatures/   a_script/, b_reports/s03_signatures.pdf, c_data/s03_signatures.xlsx
```
