# F1

The global view of HR against LR: the design, the proteome at each time point, how many proteins
each contrast moves, and whether one group's response shows up in the other.

| Panel | Shows | Reads |
|---|---|---|
| a_design | Study design and the contrast colours used throughout | `assets/design.png`; people and biopsies from `sample_weights` |
| b_pca | PC1 and PC2 within T1, T2 and T3 (blood-corrected values, proteins scaled), HR and LR as dark circles and light triangles; PERMANOVA for group per time point, and for time across all biopsies | `pca_by_time`, `permanova` |
| c_counts | Proteins per contrast at p < 0.05 (light bar), Π < 0.05 (dark bar) and FDR < 0.05 (number), up and down; dashed lines mark 2.5% of tested per direction, the 5% expected at p < 0.05 | `protein_by_contrast`, checked against `contrast_summary` |
| d_fry_transfer | Each group's training and acute responders (Π or p < 0.05, up and down) tested with fry in the other group's contrast: bar length −log10 p, colour the direction fry reports, shade the evidence, label set size and FDR. Sets under 5 proteins were not tested and are labelled "too few" | `fry_concordance`, set sizes from `protein_by_contrast` |

All sheets come from `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx`.

## Run

```sh
Rscript 04_Figures/F1/a_script/F1_composite.R          # every panel, the figure, the workbook
Rscript 04_Figures/F1/a_script/panels/b_pca.R          # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel on its own |
| `b_reports/F1_composite.pdf` | the figure, 180 × 172 mm |
| `c_data/F1.xlsx` | one sheet per panel with the numbers drawn, column meanings, inputs with md5 |
