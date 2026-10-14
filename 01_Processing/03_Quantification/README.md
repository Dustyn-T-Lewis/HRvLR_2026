# 03_Quantification

Gives every protein a value and a standard error in every run with limpa, which models missing
values as faint rather than filling them. Also writes a copy with each protein's blood slope
removed, for methods that take no covariate.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Normalisation/c_data/DAList_normalized.rds` | normalised log2 abundance with `NA` where undetected |

## Run

```sh
quarto render 01_Processing/03_Quantification/a_script/quantification.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/quantification.html` | the walkthrough: detection curve, before and after, blood check |
| `c_data/proteins.rds` | limpa's fit with standard errors; read by the protein-level contrasts |
| `c_data/DAList_limpa.rds` | the complete matrix |
| `c_data/DAList_corrected.rds` | the complete matrix with blood slopes removed; read by set, module, classification and trait steps |
| `c_data/quantification.xlsx` | curve, per-protein summaries, blood index, run correlations; read by `05_Figures/Processing/S01_QC` |

The blood removal raises, rather than lowers, the blood correlation of the 85 proteins measured in
20 or fewer runs. The walkthrough shows the numbers; the fix is pending.
