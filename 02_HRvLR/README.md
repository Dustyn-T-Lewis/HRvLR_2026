# 02_HRvLR

How do HR and LR differ? Seven steps, each one method, all on the 16 people (8 HR, 8 LR) at T1, T2
and T3. The nine contrasts and the design (one mean per group and time point, blood index as a
covariate, person as a block) are set once in `01_Protein_Contrasts`; steps 02 and 03 reuse that
fit.

| Step | Does | Writes for later steps |
|---|---|---|
| [`01_Protein_Contrasts`](01_Protein_Contrasts/README.md) | limpa `dpcDE()` with sample weights; nine contrasts with logFC, 95% CI, p, FDR and Π; trajectory F-tests; fry concordance between groups; PCA and PERMANOVA per time point | `protein_contrasts.rds`, `protein_contrasts.xlsx` |
| [`02_Set_Contrasts`](02_Set_Contrasts/README.md) | fgsea and fry on every gene set in the nine contrasts, FDR within database; Jaccard tree selection for the pathway figure | `set_contrasts.rds`, `set_contrasts.xlsx` |
| [`03_Module_Contrasts`](03_Module_Contrasts/README.md) | limma on module eigengenes and fry on modules as protein sets, same design; kME against protein t | `module_contrasts.rds`, `module_contrasts.xlsx` |
| [`04_Variance`](04_Variance/README.md) | variance shares (person, time, group, blood), DiffVar between groups, gap statistic, PCs against covariates | `variance.xlsx` |
| [`05_Module_Preservation`](05_Module_Preservation/README.md) | WGCNA modules built inside each group, preservation of each in the other | `module_preservation.rds`, `module_preservation.xlsx` |
| [`06_Classification`](06_Classification/README.md) | single-feature AUC (Wilcoxon, DeLong interval) for proteins, set scores and modules in five tasks | `classification.rds`, `classification.xlsx` |
| [`07_Group_Trait_Slopes`](07_Group_Trait_Slopes/README.md) | Spearman within each group; limma HR slope minus LR slope per feature, window and trait | `group_trait_slopes.rds`, `group_trait_slopes.xlsx` |
| [`08_Module_Detail`](08_Module_Detail/README.md) | module hubs, coupling within each group, kME rewiring between groups, per-person changes, sPLS-DA and sPLS-R | `module_detail.rds`, `module_detail.xlsx` |

Shared functions live in `00_Setup/Shared_functions/`: `person_traits.R` (per-person values and
trait windows) and `set_tree.R` (the gene-set tree of step 02).

## Inputs

`01_Processing` outputs (`quantification.rds`, `blood_correction.rds`, `gene_sets.rds`,
`set_scores.rds`, `module_building.rds`) and `00_Setup/Input/phenotype.csv`.

## Run

```sh
Rscript run_all.R 02_HRvLR
```

Steps 02 and 03 read `01_Protein_Contrasts/c_data/protein_contrasts.rds`, so 01 runs first; the others read only
`01_Processing`.

## Read with care

No protein, module or fry set test separates HR from LR at FDR < 0.05 in an HR − LR contrast or an
interaction, and one single-feature AUC of 17,455 does (a GO Slim set). fgsea, which treats proteins
as independent, calls 403 sets in the same five contrasts (four databases, GO Slim aside), 263 of
them at T3. Nominal protein counts in the HR − LR contrasts run from 0.8 to 1.4 times chance. HR and
LR are the top and bottom 8 of 25 by the composite hypertrophy score, so differences in growth
between them are built in.
