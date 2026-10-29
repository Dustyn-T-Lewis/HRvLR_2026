# 07_Group_Trait_Slopes

How much did each group change in fibre area, muscle area and strength (geometric-mean % with
95% CI per group, Welch's t-test for HR − LR, on log ratios)? And does the link between a feature
and a muscle trait differ between HR and LR? Per feature type
(proteins, set scores, eigengenes), window and trait: Spearman's rho within each group (exact p,
descriptive), and limma on feature ~ group × trait, whose interaction is HR's slope minus LR's with
its 95% interval. The pooled, label-free version is `03_Continuous/01_Trait_Associations`.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/blood_correction.rds` | blood-corrected protein values and the biopsy table |
| `01_Processing/05_Set_Scores/c_data/set_scores.rds` | singscore per set and biopsy |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets and their databases |
| `01_Processing/06_Module_Building/c_data/module_building.rds` | eigengenes |
| `00_Setup/Input/phenotype.csv` | per-person size and strength traits at T1 and T2, composite score, volume load |
| `00_Setup/Shared_functions/person_traits.R` | `trait_windows()`, `window_values()` |

## Run

```sh
quarto render 02_HRvLR/07_Group_Trait_Slopes/a_script/group_trait_slopes.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/group_trait_slopes.rds` | `trait_effects`, `counts`, `within_group`, `group_difference`, for the figures |
| `c_data/group_trait_slopes.xlsx` | the same tables |
| `b_reports/group_trait_slopes.html` | the rendered walkthrough |

## Read with care

- Trait change: HR grew in fibre area (+16%) while LR's shrank (−15%); muscle area +11% against
  +3%. Muscle area and mixed fibre area chose the groups, so those gaps are built in. Strength rose
  alike: leg press +34% against +31% (p 0.85), leg extension +17% against +13% (p 0.43).
- Slope difference: 4,900 of 136,149 tests at p < 0.05 (3.6%; chance 5%), 50 at FDR < 0.05, 44 of
  them gene sets in the training window against leg-press gain. Proteins: 3,207 of 80,301 (4.0%),
  3 at FDR < 0.05.
- Within groups: 12,927 of 292,126 Spearman tests at p < 0.05 (4.4%), 15 at FDR < 0.05. Each rests
  on 8 people, and a group split on growth has a narrowed range of growth inside it.
- The composite hypertrophy score defined the groups, so it is left out of the slope difference.
- Proteins limpa never observed at a time point carry one estimate for everyone there; they have no
  rank correlation in that window and are not counted as tested. The slope test counts them, with a
  slope difference of 0 and p = 1 (46 tests).
