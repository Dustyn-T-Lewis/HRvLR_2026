# 02_Contrasts

Fits the nine contrasts with limpa and limma.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Preprocess/03_Quantify/c_data/proteins.rds`, `01_Preprocess/03_Quantify/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 02_Proteins/02_Contrasts/a_script/02_contrasts.R` |

## Method

`fit_limpa()` in `R/helpers.R`: limpa's `dpcDE()` turns each value's standard error into a
precision weight (vooma), fits the design with subject blocked and blood as a covariate, then
`contrasts.fit()` and `eBayes(robust = TRUE)`; BH within each contrast. limpa gives every protein a
value in every sample, so every contrast tests all 1,900 proteins, and uncertain values count less.
`treat()` retests at a 1.15-fold floor. The fry concordance check takes one arm's responders and
tests them in the other arm's contrast on limpa's complete matrix. `fit.rds` holds the fit, its
design, limpa's matrix and object; the pathway and module contrasts read it.

The `unadjusted` sheet refits the nine contrasts without the blood term: the sensitivity check for
the blood covariate.

ANCOVA sits beside the interactions. Per protein and window, `lm(later ~ earlier + arm)` on the
blood-corrected values where observed, over subjects with both biopsies and at least three per arm:
T2 on T1 for training, T3 on T2 for the acute bout. The arm term is HR minus LR at equal starting
value, so an arm that starts low and catches up does not read as a different response. Residual
variances are moderated across proteins (`limma::squeezeVar`), as limma moderates its own fits. BH
within window.

Leave one subject out reruns the tests that feed the candidate tables once per subject dropped,
through `leave_one_out()` in `R/helpers.R`. A hit that stays at p < 0.05 in every rerun does not
rest on one person.

## Figures

Two overview pages: every protein across the 45 biopsies, centred per protein, with arm, time point,
fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the biopsy
traits are blank at T3, where none were measured); then log2 fold change in the nine contrasts
beside a status map. Then the p-value histogram, a volcano per contrast, and the ANCOVA volcano for both windows.

## Workbook

| Sheet | Holds |
|---|---|
| `DEP_matrix` | logFC, p, BH FDR and Storey q per protein and contrast, with samples observed. |
| `contrast_summary` | Tested, p < 0.05 and BH counts per contrast, with its role. |
| `null_calibration` | propTrueNull, median |t| and share of p below 0.2 per contrast. |
| `treat` | Proteins at BH < 0.05 under treat() with a 1.15-fold floor, per contrast. |
| `detectable_effect` | log2 effect the primary contrast detects at 80% power, median and 90th percentile. |
| `fry_concordance` | fry: one arm's responders tested in the other arm's contrast. |
| `unadjusted` | logFC and p per protein and contrast before blood adjustment, for sensitivity. |
| `ancova` | ANCOVA per protein and window: arm effect at equal starting value, t, p, BH FDR and Storey q. |
| `leave_one_out` | Per protein and HR-against-LR test: reruns with one subject left out, reruns still at p < 0.05, the largest p and the subject whose removal gave it. |
