# 08_Module_Detail

What sits underneath and around each module (`01_Processing/06_Module_Building`):

- hubs: each module's 15 highest-kME proteins with their training, acute and interaction results;
- coupling: Spearman correlation between two modules' eigengene changes across people, within HR
  and within LR, and the HR − LR difference (Fisher z);
- rewiring: each member's kME within HR biopsies and within LR biopsies;
- per person: each person's eigengene change per module, with their growth score;
- PLS: sPLS-DA (HR against LR) and sPLS regression (trait gains) on the eigengenes, 4-fold
  cross-validation repeated 5 times (folds balanced for HR and LR) with the panel size tuned inside
  each training set, 1,000 permutations, and the loadings of a model on all 16 people.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/06_Module_Building/c_data/module_building.rds` | eigengenes, membership with kME, biopsy table |
| `01_Processing/03_Blood_Correction/c_data/blood_correction.rds` | values for the within-group kME |
| `02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx` | `protein_by_contrast` for the hubs |
| `00_Setup/Input/phenotype.csv` | group, growth score and trait gains |
| `00_Setup/Shared_functions/person_traits.R` | `person_values()`, `trait_windows()` |

## Run

```sh
quarto render 02_HRvLR/08_Module_Detail/a_script/module_detail.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/module_detail.rds` | `hubs`, `coupling`, `rewiring`, `per_person`, `pls`, `pls_loadings`, for S7 (`04_Figures/F5/a_script/supp/S7_composite.R`) |
| `c_data/module_detail.xlsx` | the same tables |
| `b_reports/module_detail.html` | the rendered walkthrough |

## Read with care

- Everything here is exploratory: 8 people per group, nothing corrected across the step.
- Coupling uses Fisher's z on Spearman correlations from 8 people each; treat its p as a guide.
- Within-group kME rests on 24 biopsies from 8 people; a shift needs to be large to mean much.
- The permutation null reruns the whole procedure on shuffled outcomes, so tuning cannot inflate
  the scores. Balanced folds avoid the leave-one-out bias, where a model with no signal predicts
  against the held-out person's group and the AUC falls below 0.5.
