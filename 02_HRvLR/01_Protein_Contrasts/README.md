# 01_Protein_Contrasts

Tests every protein in nine HR-against-LR contrasts and two trajectory F-tests with limpa's
`dpcDE()`, and describes the whole proteome per time point (PCA, PERMANOVA).

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Quantification/c_data/proteins.rds` | limpa protein values and standard errors, biopsy table with the blood index |
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | blood-corrected values, for the PCA and PERMANOVA |

## Run

```sh
quarto render 02_HRvLR/01_Protein_Contrasts/a_script/protein_contrasts.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/fit.rds` | the `dpcDE()` fit, the moderated contrasts (`eb`), design and contrast matrix |
| `c_data/protein_contrasts.xlsx` | per protein and contrast results, contrast counts, trajectory F-tests, sample weights, fry carry-over (all four sets per pairing, untested ones blank), PCA per time point and per change (training, acute), PERMANOVA for each |
| `b_reports/protein_contrasts.html` | the rendered walkthrough |

## Read with care

- About 103 proteins reach p < 0.05 in any one contrast by chance. Within-group contrasts give
  155 to 204; HR − LR gives 86 at T1, 141 at T2 and 124 at T3; the interactions give 68 (training)
  and 119 (acute); the two F-tests give 98 and 84. No HR − LR or interaction test reaches
  FDR < 0.05.
- Sample weights act only in the `dpcDE()` fits. They do not differ between HR and LR (Wilcoxon
  test in the report).
