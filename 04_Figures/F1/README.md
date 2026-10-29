# F1

The global view of HR against LR. Row 1: the design and the proteome at each time point. Row 2: how
much each group grew and got stronger, and how many proteins (d) and gene sets (e) each contrast
moves. Row 3: whether one group's response shows up in the other.

| Panel | Shows | Reads |
|---|---|---|
| a_design | Study design and the contrast colours used throughout | `assets/design.png`; people and biopsies from `sample_weights` |
| b_pca | Three PCAs in a row, one point per person (blood-corrected values, proteins scaled): T1 levels, the training change (T2 − T1) and the acute change (T3 − T2). HR deep indigo circles, LR periwinkle triangles, each group in its 95% normal ellipse with its centre marked and labelled. Each carries the PERMANOVA for HR/LR (on a change: did the groups respond differently?), the gap statistic's k (1 = no natural clusters) and the adjusted Rand index of a forced two-cluster split against HR/LR (0 = chance); the subtitle gives the PERMANOVA for time across all biopsies | `pca`, `permanova`; `separation` from `04_Variance`; drawn by `style/pca_grid.R` |
| c_trait_effects | Each group's T1 to T2 change in fibre area, muscle area and strength: geometric-mean % with 95% CI on a log-ratio axis, filled at p < 0.05. Muscle area and mixed fibre area formed the groups and sit in their own block without an HR − LR p; the other traits carry Welch's HR − LR p | `trait_effects` from `07_Group_Trait_Slopes` |
| d_counts | Proteins per contrast at p < 0.05 (light bar), Π < 0.05 (dark bar) and FDR < 0.05 (number), up and down; dashed lines mark 2.5% of tested per direction, the 5% expected at p < 0.05 | `protein_by_contrast`, checked against `contrast_summary` |
| e_set_counts | Gene sets at fgsea FDR < 0.05 per contrast, up and down, over Hallmark, Reactome, KEGG and GO BP (1,364 sets, FDR within each database). Beside d: HR and LR change a similar number of proteins, but HR moves 142 sets with training against LR's 31, and 164 against 28 after the bout | `set_tests` from `02_Set_Contrasts` |
| f_fry_transfer | Each group's training and acute responders (Π or p < 0.05, up and down) tested with fry in the other group's contrast: bar length −log10 p, colour the direction fry reports, shade the evidence, label set size · FDR. Sets under 5 proteins were not tested and are labelled "too few" | `fry_concordance` |

Sheets come from `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` unless named.

## Read with care


- Panel c: HR and LR were the top and bottom 8 of 25 on muscle and mixed fibre area growth, so they
  differ there by construction. Type I and type II fibre areas were not in the score, but mixed
  fibre area is built from the same fibres, so their HR − LR gap largely follows from the
  selection. Strength gains do not differ (leg press p 0.85, leg extension p 0.43).
- Panel b: the whole-proteome response does not differ between groups (training R² 0.06, p 0.61;
  acute R² 0.08, p 0.17), though LR's responses sit in a tighter ellipse than HR's.
- Panels d and e: LR changes as many single proteins as HR (203 against 156 with training),
  while HR changes more gene sets (142 against 31).

## Supplementary figure S1

Technical supplement to F1: whether the biopsies cluster, where the variance sits, whether one
group varies more, and what drives the main principal components.

| Panel | Shows | Reads |
|---|---|---|
| a_no_clusters | Biopsy-by-biopsy correlation in the saved average-linkage order, with HR/LR, time point and person strips; beside it the gap statistic for k = 1 to 6 with its standard error and the k the one-standard-error rule picks | `biopsy_correlation`, `biopsy_order`, `gap`, `pc_scores` |
| b_variance_shares | Each protein's variance shares for person, time point, HR/LR, blood index and residual (variancePartition), with medians | `variance_shares` |
| c_diffvar | DiffVar log variance ratio (LR relative to HR) per time point, shaded where p < 0.05; counts per direction against chance (2.5% of tested each way), FDR count, and ratios beyond ±4, which are counted but not drawn. Exploratory: 8 people per group | `diffvar` |
| d_pc_covariates | R² of each covariate alone on PCs 1 to 5: time point, HR/LR, person, blood index, proteins detected, sample weight | `pc_covariates` |
| e_pca_levels | PCA of the T2 and T3 levels, one point per person, drawn as Figure 1b (`style/pca_grid.R`) with each time point's PERMANOVA, gap k and forced-split ARI | `pca`, `permanova` from `01_Protein_Contrasts`; `separation` |

Sheets come from `02_HRvLR/04_Variance/c_data/variance.xlsx` unless named.

`b_reports/supp/S1_pc_scatters.pdf` has every PC-covariate pair behind panel d, one covariate per page;
`a_script/supp/extra/pc_scatters.R` writes it when the composite runs. It is not a panel.

### Read with care


- DiffVar tests only proteins observed in all 16 biopsies of a time point; a protein one group never
  detected has modelled values with no spread, and the test would read detection.
- HR people vary more than LR people at T1 and T3 (panel c). At p < 0.05, 50 proteins are more
  variable in HR against 11 in LR at T1, and 87 against 5 at T3; chance gives 27 and 17 each way.
  T2 is even (31 against 36, chance 22). No protein passes FDR < 0.05, and each group has 8
  people.

## Run

```sh
Rscript run_all.R 04_Figures/F1                       # F1 and S1
Rscript 04_Figures/F1/a_script/F1_composite.R        # F1 alone
```

## Outputs

| File | Holds |
|---|---|
| `a_script/`, `b_reports/` | F1's panel scripts and PDFs (`a_*.R`, `a_*.pdf`, ...) and its composite (`F1_composite.pdf`: the figure, 250 × 235 mm (landscape)) |
| `a_script/supp/`, `b_reports/supp/` | S1: `S1_*.R`, `S1_*.pdf` and `S1_composite.pdf` (the figure, 250 × 167 mm landscape) |
| `c_data/F1.xlsx`, `c_data/S1.xlsx` | one sheet per panel, the numbers drawn |
