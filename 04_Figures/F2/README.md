# F2: proteins

Responders first, then the exploratory HR − LR candidates. The nominal call is Π = p^|logFC|
< 0.05 (Xiao et al. 2014), drawn beside logFC, p and FDR; a Π hit is not called significant.

| Panel | Shows |
|---|---|
| a `a_responders` | Protein dot-matrix over the nine contrasts. Rows: proteins at Π < 0.05 in a group's training or acute contrast, in three blocks (shared by HR and LR, HR only, LR only), the top 10 per block by p (shared proteins rank on the weaker group's p). |
| b `b_hr_lr_candidates` | Same grammar. Rows: proteins at Π < 0.05 in an HR − LR contrast (T1, T2, T3) or an interaction, in the block of the hit with the smallest p, the top 12 per block. Exploratory: the subtitle gives the range of p < 0.05 counts over these five contrasts against the number chance gives, and the count at FDR < 0.05 (none). |
| c `c_trajectories` | Blood-corrected values per person at T1, T2 and T3 for the four candidates with the smallest p, HR and LR shaded, with group means. Each header gives logFC, Π, p and FDR. Open points are biopsies where no precursor was seen and limpa modelled the value. |

Dot-matrix code: colour is logFC (PiYG, capped at the median |logFC| of all Π hits so most dots
stay readable; a and b share the cap and one key; the workbook holds the uncapped values), size is
−log10 p (capped at 4), a filled dot is a Π hit and an open one is not, a black outline marks
FDR < 0.05. Column labels are the group (training, acute), the time point (HR − LR) or the window
(interaction). Each row reads "gene (n/48)", the biopsies where the protein was observed. The grid
and the selection rules are in `04_Figures/style/protein_grid.R`.

## Why rows are ranked by p, not Π

Π decides which proteins may enter; p decides which are drawn. A protein seen in few biopsies takes
most of its values from limpa's detection model, and a modelled floor can give it a large logFC.
Π = p^|logFC| rewards a large logFC, so ranking on Π filled panel b with such proteins (CNN1, 17/48,
led it at p = 0.011). The p from `dpcDE` already carries limpa's quantification uncertainty, which
is how limma ranks (`topTable`). Proteins seen in fewer than 24 of 48 biopsies are 13% of those
tested but a third of all Π hits.

## Read with care

- Panel b is at chance: each HR − LR contrast has about as many proteins at p < 0.05 as chance
  gives, and none passes FDR < 0.05. Its rows are leads.
- "HR only" in panel a means a Π hit in HR and not in LR; it is not a test of HR against LR. The
  interaction columns are that test.
- Sparse proteins still enter when their p is small (RBP1 12/48, EIF4B 11/48); their counts are on
  the labels.

## Inputs

| File | Used for |
|---|---|
| `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` | `protein_by_contrast` (a, b, c), `contrast_summary` (b subtitle) |
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | values and observation counts (c) |

## Run

```sh
Rscript 04_Figures/F2/a_script/F2_composite.R             # every panel, the composite, the workbook
Rscript 04_Figures/F2/a_script/panels/a_responders.R      # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel alone |
| `b_reports/F2_composite.pdf` | the figure, 250 × 140 mm landscape |
| `c_data/F2.xlsx` | one sheet per panel with the numbers drawn, column dictionary, input md5 |
