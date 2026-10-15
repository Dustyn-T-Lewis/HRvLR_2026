# F2: proteins

Responders first, then the exploratory HR − LR candidates. The nominal call is Π = p^|logFC|
< 0.05 (Xiao et al. 2014), drawn beside logFC, p and FDR; a Π hit is not called significant.

| Panel | Shows |
|---|---|
| a `a_responders` | Protein dot-matrix over the nine contrasts. Rows: proteins at Π < 0.05 in a group's training or acute contrast, in three blocks (shared by HR and LR, HR only, LR only), the top 10 per block by Π (shared proteins rank on the weaker group's Π). |
| b `b_hr_lr_candidates` | Same grammar. Rows: proteins at Π < 0.05 in an HR − LR contrast (T1, T2, T3) or an interaction, in the block of the contrast with the smallest Π, the top 12 per block. Exploratory: the subtitle gives the range of p < 0.05 counts over these five contrasts against the number chance gives, and the count at FDR < 0.05 (none). |
| c `c_trajectories` | Blood-corrected values per person at T1, T2 and T3 for the four candidates with the smallest Π in panel b's contrasts, HR and LR shaded, with group means. Each header gives logFC, Π, p and FDR. Open points are biopsies where no precursor was seen and limpa modelled the value. |

Dot-matrix code: colour is logFC (PiYG, capped at the median |logFC| of the drawn Π hits so most
dots stay readable; the workbook holds the uncapped values), size is −log10 p (capped at 4), a
filled dot is a Π hit and an open one is not, a black outline marks FDR < 0.05. Column labels are
the group (training, acute), the time point (HR − LR) or the window (interaction). Each row reads
"gene (n/48)", the biopsies where the protein was observed. The grid is drawn by
`04_Figures/style/protein_grid.R`.

## Read with care

Panel b is led by proteins seen in few biopsies (CNN1 17/48, RBP1 12/48). Most of their values
come from limpa's detection model, and that modelled floor gives them a large logFC, so a tiny Π
on little data. Their p values sit near 0.01 and none passes FDR < 0.05. The rows are kept, with
their counts shown, rather than filtered.

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
| `b_reports/F2_composite.pdf` | the figure, 180 mm wide |
| `c_data/F2.xlsx` | one sheet per panel with the numbers drawn, column dictionary, input md5 |
