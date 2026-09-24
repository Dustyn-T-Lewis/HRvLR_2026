# 02 · Differential Expression

Fits the normalised DAList and tests each protein. Three sub-stages, so a design problem surfaces
before anything is fitted.

```
01_Preprocess/02_Normalization/c_data/DAList_normalized.rds
  01_Design        design formula, nine contrasts, one diagnostic  -> design.rds
  02_Differential  fit, results, classification, hit matrices      -> fit.rds
  03_Phenotype     every protein against the ten phenotypes        -> phenotype.rds
```

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd       --output-dir ../b_reports
```

## The nine comparisons

| Name | Role | Compares |
|---|---|---|
| `Training_Interaction` | primary | HR training change minus LR training change |
| `Acute_Interaction` | secondary | HR acute change minus LR acute change |
| `Baseline_HRvLR` | floor | the arms before training |
| `Training_HR`, `Training_LR` | descriptive | T2 − T1 within an arm |
| `Acute_HR`, `Acute_LR` | descriptive | T3 − T2 within an arm |
| `Trained_HRvLR`, `Acute_HRvLR` | descriptive | the arms at T2, at T3 |

`Training_Interaction` equals `Trained_HRvLR` minus `Baseline_HRvLR` exactly, so these are not
nine independent questions.

## Notes

No contrast has a BH hit. The lowest adjusted p across all nine is 0.071 (Acute_HR). The
primary contrast could detect a median effect of 1.04 log2 units at 80% power, so the empty
interaction cannot exclude effects smaller than about twofold.

The floor is not a negative control. The HR/LR label was cut from the training outcome, so a
real baseline difference would be predictive. It has no hits, and its null-calibration numbers are
the between-arm null reference.

Subject is a random effect. HR and LR are different people, so a fixed subject term would
absorb every between-arm contrast. The within-subject correlation is 0.189.
