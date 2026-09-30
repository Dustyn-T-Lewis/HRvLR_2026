# 02_Contrasts

Fits the nine contrasts with proteoDA and limma on the unimputed matrix.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `01_Preprocess/02_Normalization/c_data/DAList_unadjusted.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf`, and proteoDA's per-contrast HTML reports in `b_reports/` (ignored) |
| Run | `Rscript 02_Differential_Expression/02_Contrasts/a_script/02_contrasts.R` |

## Method

`fit_limma_model()` fits every protein on the samples where it was seen, subject blocked, then
`eBayes(robust = TRUE)`; `extract_DA_results()` applies BH within each contrast. A protein with an
empty cell in a contrast is not tested there. `treat()` retests at a 1.15-fold floor. The fry
concordance check takes one arm's responders and tests them in the other arm's contrast on the
imputed matrix, since fry takes no missing value. `fit.rds` holds the full proteoDA result; the
pathway and module contrasts read it.

The `unadjusted` sheet refits the nine contrasts on the matrix before blood adjustment, without the
blood term: the sensitivity check for that adjustment.

ANCOVA sits beside the interactions. Per protein and window, `lm(later ~ earlier + arm)` over
subjects with both biopsies and at least three per arm: T2 on T1 for training, T3 on T2 for the
acute bout. The arm term is HR minus LR at equal starting value, so an arm that starts low and
catches up does not read as a different response. BH within window.

## Figures

Two overview pages: every protein across the 45 biopsies, centred per protein, with arm, time point,
fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the biopsy
traits are blank at T3, where none were measured); then log2 fold change in the nine contrasts
beside a status map. Then the p-value histogram, a volcano per contrast, and the ANCOVA volcano for both windows.

## Workbook

| Sheet | Holds |
|---|---|
| `DEP_matrix` | logFC, p and BH FDR per protein and contrast, with samples observed. |
| `contrast_summary` | Tested, p < 0.05 and BH counts per contrast, with its role. |
| `null_calibration` | propTrueNull, median |t| and share of p below 0.2 per contrast. |
| `treat` | Proteins at BH < 0.05 under treat() with a 1.15-fold floor, per contrast. |
| `detectable_effect` | log2 effect the primary contrast detects at 80% power, median and 90th percentile. |
| `fry_concordance` | fry: one arm's responders tested in the other arm's contrast. |
| `unadjusted` | logFC and p per protein and contrast before blood adjustment, for sensitivity. |
| `ancova` | ANCOVA per protein and window: arm effect at equal starting value, t, p and BH FDR. |
