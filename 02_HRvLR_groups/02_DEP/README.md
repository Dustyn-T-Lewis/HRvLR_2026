# 02_Contrasts

Fits the nine contrasts with limpa and limma.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Trunk/03_Quantify/c_data/proteins.rds`, `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/02_DEP/a_script/02_dep.R` |

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
limpa's complete blood-corrected values, over subjects with both biopsies and at least three per arm:
T2 on T1 for training, T3 on T2 for the acute bout. The arm term is HR minus LR at equal starting
value, so an arm that starts low and catches up does not read as a different response. Residual
variances are moderated across proteins (`limma::squeezeVar`), as limma moderates its own fits. BH
within window.

The whole proteome is tested once, beside the per-protein fits. PCA runs on the complete
blood-corrected matrix with proteins scaled, and PERMANOVA (`vegan::adonis2`, 9,999 permutations)
on its Euclidean distances. Time is tested within subject, so permutations never move a sample to
another person; arm is tested on one mean profile per subject, and at each time point. Time
explains 6.3% of variance (p = 0.0009), most of it after the bout (acute 5.6%, p = 0.0005;
training 4.2%, p = 0.052). Arm explains none beyond chance (subject means p = 0.42; T3
p = 0.21). The same test on the six size and strength traits at T1 and T2 (time within subject)
gives R² 3.6%, p = 0.037; arm is not tested there, because HR and LR were split on how much
those traits changed.

## Figures

Two overview pages: every protein across the 45 biopsies, centred per protein, with arm, time point,
fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the biopsy
traits are blank at T3, where none were measured); then log2 fold change in the nine contrasts
beside a status map. Then the p-value histogram, a volcano per contrast, the ANCOVA volcano for both windows, and the PCA.

## Workbook

| Sheet | Holds |
|---|---|
| `DEP_matrix` | logFC, p, BH FDR and Storey q per protein and contrast, with samples observed. |
| `contrast_summary` | Tested, p < 0.05 and BH counts per contrast, with its role. |
| `null_calibration` | propTrueNull, median |t| and share of p below 0.2 per contrast. |
| `treat` | Proteins at BH < 0.05 under treat() with a 1.15-fold floor, per contrast. |
| `detectable_effect` | log2 effect the primary contrast detects at 80% power, median and 90th percentile. |
| `fry_concordance` | fry: one arm's responders in the other arm's same contrast, and each arm's training sets on both arms' acute contrasts; with mixed p, mean log2 FC in both contrasts and the share of members keeping sign. A training set on its own arm's acute contrast shares T2 and reads as a reversal by construction. |
| `set_trajectory` | Each arm's p < 0.05 training sets: mean change from each subject's own T1, with SE, per arm and time point (blood-corrected complete matrix). |
| `unadjusted` | logFC and p per protein and contrast before blood adjustment, for sensitivity. |
| `ancova` | ANCOVA per protein and window: arm effect at equal starting value, t, p, BH FDR and Storey q. |
| `pca` | PC1 and PC2 per sample, from the complete blood-corrected matrix with proteins scaled. |
| `pca_variance` | Share of variance on the first four components. |
| `permanova` | PERMANOVA: proteome time within subject per window and arm, arm on subject means and per time point; phenotype time within subject. |
