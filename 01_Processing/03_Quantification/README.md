# 03_Quantification

How does every protein get a value in every sample without treating estimates as measurements?

| | |
|---|---|
| Reads | `01_Processing/02_Normalisation/c_data/DAList_normalized.rds` |
| Writes | `c_data/DAList_corrected.rds`, `c_data/DAList_limpa.rds`, `c_data/proteins.rds`, `c_data/quantification.xlsx`, `b_reports/quantification_report.pdf` |
| Run | `Rscript 01_Processing/03_Quantification/a_script/quantification.R` |

## Key result

- 12.3% of values are missing; the detection curve has intercept -10.45 and slope 0.67 (logit scale).
- 1,900 proteins get a value in all 45 biopsies; 931 were measured in every one, and the median protein in 44.

## Method

- `limpa::dpc()` fits one detection-probability curve: the chance a value is detected, as a logistic function of its log2 intensity.
- `limpa::dpcQuantByRow()` gives each protein a value and a standard error in every sample. Observed values keep a small standard error; estimated ones get a wide one, so later fits weight them less.
- `DAList_corrected.rds` removes each protein's blood slope (`limma::removeBatchEffect`, design `~ 0 + timepoint`) for the methods that cannot take a covariate.

## Read with care

- `prop_observed` in `proteins` says how much of a row was measured; a low share rests on the curve more than on data.

## Sheets

`quantification.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `curve` | DPC intercept and slope on the logit scale, and the share of values missing. |
| `proteins` | Per protein: values observed, mean log2 abundance where observed and over the complete matrix, and the median standard error where observed and where missing. |
| `blood` | Per sample: blood index (HBB and HBA1 against the sample median) and its centred value. |
| `sample_correlation` | Pearson correlation between samples, each protein centred, on the blood-corrected matrix. |
