# 05_Modules

A copy of `02_HRvLR_groups/04_Modules` run on the muscle size clusters: groups Hi and Lo from `01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

Tests the modules on the nine contrasts, by eigengene and by fry, and reads membership against
protein-level results.

| | |
|---|---|
| Reads | `01_Build/c_data/modules.rds`, `03_Cluster_groups/Muscle_size/02_Design/c_data/design.rds`, `03_Cluster_groups/Muscle_size/03_DEP/c_data/fit.rds`, `04_Continuous_no_group/01_Associate/01_Proteins/c_data/associate.rds`, `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/05_Modules/a_script/05_modules.R` |

## Method

Eigengenes get the protein design and subject block with their own within-subject correlation,
`eBayes(robust = TRUE)` and BH within contrast. fry tests each module as a protein set on the
limpa's complete matrix, FDR over the modules within contrast. Membership against significance is the
Spearman correlation, within module, between a member's kME and its protein-level moderated t, and
between kME and its training change-score rho with each outcome. Members share a module, so these
correlations are descriptive.

ANCOVA on eigengenes sits beside the interactions: per module and window, `lm(later ~ earlier +
arm)` over subjects with both biopsies, T2 on T1 for training and T3 on T2 for the acute bout, variances
moderated across modules (`squeezeVar`). BH within window.

## Figures

Two overview pages: the eigengenes across the 45 biopsies, with arm, time point, fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the biopsy traits are blank at
T3, where none were measured); then eigengene log2 fold change in the nine contrasts beside a status
map. Then eigengene and fry tiles, membership against each contrast and each outcome, and kME
against moderated t for `Training_Interaction` in every module, and the ANCOVA volcano for both windows.

## Workbook

| Sheet | Holds |
|---|---|
| `eigengene_tests` | limma on module eigengenes: logFC, moderated t, p, BH FDR and Storey q per contrast. |
| `fry_tests` | fry on each module as a protein set: direction, p and FDR per contrast. |
| `membership_contrast` | Spearman rho between member kME and protein moderated t, per module and contrast. |
| `membership_phenotype` | Spearman rho between member kME and protein rho with each phenotype over training. |
| `correlation` | Within-subject correlation on the eigengenes and on limpa's matrix. |
| `ancova` | ANCOVA per module and window on eigengenes: arm effect at equal start, t, p, BH FDR, Storey q. |
