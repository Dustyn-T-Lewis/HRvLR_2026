# 01_Trait_Associations

Which proteins, gene-set scores and module eigengenes follow muscle size and strength with all 16
people pooled and no group label? Four analyses per feature type: a sample model on the T1 and T2
biopsies (within- and between-person trait terms, person blocked), a change score per window and
trait (limma slope, Spearman rho, bicor), a two-predictor fit (growth + strength gain) and their
interaction.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | blood-corrected protein values and the biopsy table |
| `01_Processing/05_Set_Scores/c_data/singscore.rds` | singscore per set and biopsy |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets and their databases |
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes |
| `00_Setup/Input/phenotype.csv` | per-person size and strength traits at T1 and T2, composite score, volume load |
| `02_HRvLR/helpers.R` | `trait_windows()`, `window_values()` |

## Run

```sh
quarto render 03_Continuous/01_Trait_Associations/a_script/trait_associations.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/associations.rds` | `counts`, `model`, `change_score`, `two_predictor`, `interaction`; steps 02 and 03 read the protein change scores |
| `c_data/trait_associations.xlsx` | the same five tables |
| `b_reports/trait_associations.html` | the rendered walkthrough |

## Read with care

- Six windows × 42 window-trait screens. Change score, proteins: 3,682 of 86,478 tests at
  p < 0.05 (4.3%; chance 5%), 3 at FDR < 0.05 (EHD1 with T1 fibre area, SEC31A and APOBEC2 with
  T2 muscle area). Sets: 2,930 of 59,724 (4.9%), 144 at FDR < 0.05, 118 of them T1 sets against
  T1 leg press. Modules: 18 of 546, none.
- Sample model: 1,225 of 24,708 protein tests at p < 0.05 (5.0%), 4 at FDR < 0.05.
- HR and LR are the top and bottom 8 of 25 by growth, so correlations across these 16 run stronger
  than they would across all 25. Pooling avoids a split, not confounding.
- Proteins limpa never observed at a time point carry one estimate for everyone there; their rho
  and bicor are empty in that window.
