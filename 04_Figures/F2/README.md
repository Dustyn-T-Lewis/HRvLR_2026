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

## Supplementary figure S2: can one feature tell HR from LR?

Panel a: how well one protein, set score or module eigengene separates the groups, from
`02_HRvLR/06_Classification`; no model is fitted and no AUC is measured on held-out people. Panels
b and c: the exploratory HR − LR protein candidates from `02_HRvLR/01_Protein_Contrasts`, drawn with
the F2 dot-matrix (`style/protein_grid.R`).

| Panel | Shows |
|---|---|
| a `a_auc_top` | Per task (T1, T2, T3 levels; T2 − T1 and T3 − T2 changes) and feature type, the three features with the smallest Wilcoxon p (ties broken by AUC furthest from 0.5; the step's `top` table): AUC with DeLong's 95% interval, no interval at AUC 0 or 1. Dashed line: 0.5, chance. Shade: which group is higher. Shape: filled at FDR < 0.05, open at p < 0.05 only, a cross otherwise. Under each cell: features at p < 0.05 against the number chance gives. |
| b `b_hr_lr_candidates` | The F2 dot-matrix. Rows: proteins at Π < 0.05 in an HR − LR contrast (T1, T2, T3) or an interaction, in the block of the hit with the smallest p, the top 12 per block. Exploratory: the subtitle gives the range of p < 0.05 counts over these five contrasts against the number chance gives, and the count at FDR < 0.05 (none). |
| c `c_trajectories` | Blood-corrected values per person at T1, T2 and T3 for the four candidates with the smallest p, HR and LR shaded, with group means. Each header gives logFC, Π, p and FDR. Open points are biopsies where no precursor was seen and limpa modelled the value. |

Set names come from `pathway_label()` (`style/pathway_grid.R`); modules are named by their leading
set with the colour after it (`style/module_labels.R`).

### Inputs


| File | Used for |
|---|---|
| `02_HRvLR/06_Classification/c_data/classification.rds` | `top` (the features drawn) and `counts` |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | module names (a) |
| `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` | `protein_by_contrast`, `contrast_summary`, `sample_weights` (b, c) |
| `01_Processing/03_Blood_Correction/c_data/blood_correction.rds` | values and observation counts (c) |

### Read with care


- No single feature tells the groups apart. Proteins reach p < 0.05 at or below chance (56 to 127
  of 2,059, chance 103), and none passes FDR < 0.05. Set scores at T3 are the exception on count
  (121 of 1,422 against 71 by chance), with none at FDR; one set passes FDR for the training change.
- AUC with 8 people per group is coarse: each step is 1/64, and an AUC near 0 or 1 often rests on
  one or two people.
- Panel b is at chance: each HR − LR contrast has about as many proteins at p < 0.05 as chance
  gives, and none passes FDR < 0.05. Its rows are leads.
- Sparse proteins still enter when their p is small (RBP1 12/48, EIF4B 11/48); their counts are on
  the labels.

## Run

```sh
Rscript run_all.R 04_Figures/F2                       # F2 and S2
Rscript 04_Figures/F2/a_script/F2_composite.R        # F2 alone
```

## Outputs

| File | Holds |
|---|---|
| `a_script/`, `b_reports/` | F2's panel scripts and PDFs (`a_*.R`, `a_*.pdf`, ...) and its composite (`F2_composite.pdf`: the figure, 250 × 140 mm landscape, one panel) |
| `a_script/supp/`, `b_reports/supp/` | S2: `S2_*.R`, `S2_*.pdf` and `S2_composite.pdf` (the figure, 250 × 240 mm: a across the page, b and c below) |
| `c_data/F2.xlsx`, `c_data/S2.xlsx` | one sheet per panel, the numbers drawn |
