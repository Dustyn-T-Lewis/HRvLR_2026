# Figures

The HR/LR main figures and supplements, drawn from the saved outputs of steps 01 to 07.

| Folder | Status | Asks | Draws from |
|---|---|---|---|
| [`F1_Overview`](F1_Overview/README.md) | active | What does the design look like, and how does the whole proteome respond before any single result? | `02_Protein_Contrasts`, `03_Set_Contrasts`, `phenotype.csv` |
| [`F2_Programmes`](F2_Programmes/README.md) | parked draft | Which co-expression programmes move with training and the bout, and do the groups differ in them? | `01_Design_Contrasts`, `03_Set_Contrasts`, `04_Module_Contrasts`, 01_Processing modules |
| [`F3_People`](F3_People/README.md) | parked draft | Can one feature tell the groups apart, and does a protein's link to strength gain differ between them? | `06_Classification`, `07_Trait_Associations`, the continuous analysis |
| [`S01_Volcano`](S01_Volcano/README.md) | active | Which proteins and which fgsea sets move in each within-group and interaction contrast? | `03_Set_Contrasts` |
| [`S02_NES_Scatter`](S02_NES_Scatter/README.md) | active | Do HR and LR shift the same gene sets, set by set? | `03_Set_Contrasts` |
| [`S03_Signatures`](S03_Signatures/README.md) | active | Supplement to F1: the signature ranks with the stricter Π sets, and the training proteins followed through the bout. | `02_Protein_Contrasts`, `03_Set_Contrasts` |

Parked means the script runs (`run_all.R` includes it) and its outputs are rebuilt, but the figure
is not final.

- The scripts read saved outputs only; any computation they add is a description with no new test
  (F1's phenotype PCA, F3's Spearman across proteins, S02's Spearman of NES across sets).
- Each writes `b_reports/<name>.pdf` and `c_data/<name>.xlsx` with the numbers each panel draws;
  F2 and F3 add `panel_sources`, naming the file, sheet and rule behind every panel.
- Colours, the contrast palette and panel builders come from `R/figure_theme.R`,
  `R/figure_panels.R` and `R/figure_overview.R`, which the cluster folders' F1 scripts share.
- F3 reads `04_Continuous_no_group`, so `run_all.R` runs these figures last.
