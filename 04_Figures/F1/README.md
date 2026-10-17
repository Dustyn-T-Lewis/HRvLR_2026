# F1

The global view of HR against LR. Row 1: the design and the proteome at each time point. Row 2: how
much each group grew and got stronger, how many proteins each contrast moves, and whether one
group's response shows up in the other.

| Panel | Shows | Reads |
|---|---|---|
| a_design | Study design and the contrast colours used throughout | `assets/design.png`; people and biopsies from `sample_weights` |
| b_pca | PC1 and PC2 within T1, T2 and T3, three in a row (blood-corrected values, proteins scaled), HR and LR as dark circles and light triangles. Each carries the PERMANOVA for HR/LR, the gap statistic's k (1 = no natural clusters) and the adjusted Rand index of a forced two-cluster split against HR/LR (0 = chance); the subtitle gives the PERMANOVA for time across all biopsies | `pca_by_time`, `permanova`; `separation` from `04_Variance` |
| c_trait_effects | Each group's T1 to T2 change in fibre area, muscle area and strength: geometric-mean % with 95% CI on a log-ratio axis, filled at p < 0.05. Muscle area and mixed fibre area formed the groups and sit in their own block without an HR − LR p; the other traits carry Welch's HR − LR p | `trait_effects` from `07_Group_Trait_Slopes` |
| d_counts | Proteins per contrast at p < 0.05 (light bar), Π < 0.05 (dark bar) and FDR < 0.05 (number), up and down; dashed lines mark 2.5% of tested per direction, the 5% expected at p < 0.05 | `protein_by_contrast`, checked against `contrast_summary` |
| e_fry_transfer | Each group's training and acute responders (Π or p < 0.05, up and down) tested with fry in the other group's contrast: bar length −log10 p, colour the direction fry reports, shade the evidence, label set size · FDR. Sets under 5 proteins were not tested and are labelled "too few" | `fry_concordance` |

Sheets come from `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` unless named.

## Run

```sh
Rscript 04_Figures/F1/a_script/F1_composite.R          # every panel, the figure, the workbook
Rscript 04_Figures/F1/a_script/panels/b_pca.R          # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel on its own |
| `b_reports/F1_composite.pdf` | the figure, 180 × 185 mm |
| `c_data/F1.xlsx` | one sheet per panel with the numbers drawn, column meanings, inputs with md5 |

## Read with care

- Panel c: HR and LR were the top and bottom 8 of 25 on muscle and mixed fibre area growth, so they
  differ there by construction. Type I and type II fibre areas were not in the score, but mixed
  fibre area is built from the same fibres, so their HR − LR gap largely follows from the
  selection. Strength gains do not differ (leg press p 0.85, leg extension p 0.43).
- Panel d against F3: LR changes as many single proteins as HR (203 against 156 with training),
  while HR changes far more gene sets.
