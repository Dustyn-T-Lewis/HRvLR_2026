# F2: proteins

The proteins each group changes with training or the bout. The exploratory HR − LR candidates are
in S2 (b, c). The nominal call is Π = p^|logFC|
< 0.05 (Xiao et al. 2014), drawn beside logFC, p and FDR; a Π hit is not called significant.

| Panel | Shows |
|---|---|
| a `a_responders` | Protein dot-matrix over the nine contrasts. Rows: proteins at Π < 0.05 in a group's training or acute contrast, in three blocks (shared by HR and LR, HR only, LR only), the top 10 per block by p (shared proteins rank on the weaker group's p). |

Dot-matrix code: colour is logFC (PiYG, capped at the median |logFC| of all Π hits so most dots
stay readable; the HR − LR candidates in S2 b use the same cap; the workbook holds the uncapped
values), size is
−log10 p (capped at 4), a filled dot is a Π hit and an open one is not, a black outline marks
FDR < 0.05. Column labels are the group (training, acute), the time point (HR − LR) or the window
(interaction). Each row reads "gene (n/48)", the biopsies where the protein was observed. The grid
and the selection rules are in `04_Figures/style/protein_grid.R`.

## Why rows are ranked by p, not Π

Π decides which proteins may enter; p decides which are drawn. A protein seen in few biopsies takes
most of its values from limpa's detection model, and a modelled floor can give it a large logFC.
Π = p^|logFC| rewards a large logFC, so ranking on Π filled the HR − LR candidates (now S2 b) with
such proteins (CNN1, 17/48,
led it at p = 0.011). The p from `dpcDE` already carries limpa's quantification uncertainty, which
is how limma ranks (`topTable`). Proteins seen in fewer than 24 of 48 biopsies are 13% of those
tested but a third of all Π hits.

## Read with care

- "HR only" in panel a means a Π hit in HR and not in LR; it is not a test of HR against LR. The
  interaction columns are that test.

## Inputs

| File | Used for |
|---|---|
| `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` | `protein_by_contrast`, `sample_weights` (biopsy count) |

## Run

```sh
Rscript 04_Figures/F2/a_script/F2_composite.R             # every panel, the composite, the workbook
Rscript 04_Figures/F2/a_script/panels/a_responders.R      # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel alone |
| `b_reports/F2_composite.pdf` | the figure, 250 × 140 mm landscape, one panel |
| `c_data/F2.xlsx` | one sheet per panel with the numbers drawn, column dictionary, input md5 |
