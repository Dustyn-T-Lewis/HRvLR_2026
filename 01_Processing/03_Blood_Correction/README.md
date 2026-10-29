# 03_Blood_Correction

Removes each protein's blood slope from its observed values, for methods that cannot carry the
blood index as a covariate (modules, set scores, correlations, PCA). Three ways of removing the
slope are compared first; the one kept is the only one that lowers blood correlation for every
protein group.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Quantification/c_data/quantification.rds` | limpa protein values, standard errors, observed precursor counts, blood index |

## Run

```sh
quarto render 01_Processing/03_Blood_Correction/a_script/blood_correction.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/blood_correction.rds` | the protein `EList` with `E` blood-corrected on observed values |
| `c_data/blood_correction.xlsx` | median \|rho\| with blood by detection band for the uncorrected values and the three variants |
| `b_reports/blood_correction.html` | the rendered walkthrough |

## Read with care

- Removing the slope from every value (`removeBatchEffect`) or with precision weights raises the
  blood correlation of proteins observed in 20 or fewer biopsies (median |rho| 0.14 to 0.41 and
  0.47): their slope rests mostly on values estimated from the detection curve. The kept variant
  fits and removes the slope on observed values only.
- The comparisons in `02_HRvLR` do not read this file; they carry blood as a covariate.
