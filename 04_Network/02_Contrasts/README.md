# 02_Contrasts

Tests the modules on the nine contrasts, by eigengene and by fry, and reads membership against
protein-level results.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `02_Differential_Expression/02_Contrasts/c_data/fit.rds`, `02_Differential_Expression/04_Associate/c_data/associate.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 04_Network/02_Contrasts/a_script/02_contrasts.R` |

## Method

Eigengenes get the protein design and subject block with their own within-subject correlation,
`eBayes(robust = TRUE)` and BH within contrast. fry tests each module as a protein set on the
imputed matrix, FDR over the twelve modules within contrast. Membership against significance is the
Spearman correlation, within module, between a member's kME and its protein-level moderated t, and
between kME and its training change-score rho with each outcome. Members share a module, so these
correlations are descriptive.

## Figures

Two overview pages: the eigengenes across the 45 biopsies, with arm, time point, type I share, fCSA
mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the biopsy traits are blank at
T3, where none were measured); then eigengene log2 fold change in the nine contrasts beside a status
map. Then eigengene and fry tiles, membership against each contrast and each outcome, and kME
against moderated t for `Training_Interaction` in every module.

## Workbook

| Sheet | Holds |
|---|---|
| `eigengene_tests` | limma on module eigengenes: logFC, moderated t, p and BH FDR per contrast. |
| `fry_tests` | fry on each module as a protein set: direction, p and FDR per contrast. |
| `membership_contrast` | Spearman rho between member kME and protein moderated t, per module and contrast. |
| `membership_phenotype` | Spearman rho between member kME and protein rho with each phenotype over training. |
| `correlation` | Within-subject correlation on the eigengenes and on the imputed matrix. |
