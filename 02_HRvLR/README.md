# 02_HRvLR

How do HR and LR differ? Seven steps, each one method, all on the 16 people (8 HR, 8 LR) at T1, T2
and T3. The nine contrasts and the design (one mean per group and time point, blood index as a
covariate, person as a block) are set once in `01_Protein_Contrasts`; steps 02 and 03 reuse that
fit.

| Step | Does | Writes for later steps |
|---|---|---|
| [`01_Protein_Contrasts`](01_Protein_Contrasts/README.md) | limpa `dpcDE()` with sample weights; nine contrasts with logFC, 95% CI, p, FDR and Π; trajectory F-tests; fry concordance between groups; PCA and PERMANOVA per time point | `fit.rds`, `protein_contrasts.xlsx` |
| [`02_Set_Contrasts`](02_Set_Contrasts/README.md) | fgsea and fry on every gene set in the nine contrasts, FDR within database; Jaccard tree selection for the pathway figure | `set_tests.rds`, `set_contrasts.xlsx` |
| [`03_Module_Contrasts`](03_Module_Contrasts/README.md) | limma on module eigengenes and fry on modules as protein sets, same design; kME against protein t | `module_tests.rds`, `module_contrasts.xlsx` |
| [`04_Variance`](04_Variance/README.md) | variance shares (person, time, group, blood), DiffVar between groups, gap statistic, PCs against covariates | `variance.xlsx` |
| [`05_Module_Preservation`](05_Module_Preservation/README.md) | WGCNA modules built inside each group, preservation of each in the other | `preservation.rds`, `module_preservation.xlsx` |
| [`06_Classification`](06_Classification/README.md) | single-feature AUC (Wilcoxon, DeLong interval) for proteins, set scores and modules in five tasks | `classification.rds`, `classification.xlsx` |
| [`07_Group_Trait_Slopes`](07_Group_Trait_Slopes/README.md) | Spearman within each group; limma HR slope minus LR slope per feature, window and trait | `group_trait_slopes.rds`, `group_trait_slopes.xlsx` |

`helpers.R` holds `person_values()` (one column per person: a level or a change) and the trait
windows. Steps 06 and 07 source it, as do `03_Continuous` steps 01 and 03.

## Inputs

`01_Processing` outputs (`proteins.rds`, `proteins_corrected.rds`, `gene_sets.rds`,
`singscore.rds`, `modules.rds`) and `00_Setup/Input/phenotype.csv`.

## Run

```sh
Rscript run_all.R 02_HRvLR
```

Steps 02 and 03 read `01_Protein_Contrasts/c_data/fit.rds`, so 01 runs first; the others read only
`01_Processing`.

## Read with care

No protein, module or fry set test separates HR from LR at FDR < 0.05 in an HR − LR contrast or
an interaction, and one single-feature AUC of 17,470 does (a GO Slim set). fgsea, which treats
proteins as independent, calls hundreds of sets in the same contrasts. Nominal protein counts in
the HR − LR contrasts run from 0.8 to 1.4 times chance. HR and LR are the top and bottom 8 of 25 by
the composite hypertrophy score, so differences in growth between them are built in.
