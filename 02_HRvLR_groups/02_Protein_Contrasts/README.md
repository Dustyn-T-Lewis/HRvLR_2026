# 02_Protein_Contrasts

Which proteins change between time points, and which differ between the groups?

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/proteins.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/protein_contrasts.xlsx`, `b_reports/protein_contrasts_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/02_Protein_Contrasts/a_script/protein_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_protein_contrasts.R`, shared with the cluster groupings |

## Key result

- 1,900 proteins tested in each of the nine contrasts; 0 at FDR < 0.05 across all nine (lowest FDR 0.051).
- Proteins at p < 0.05 run from 41 (Training_Interaction) to 189 (Acute_HR); chance alone gives about 95 per contrast.
- PERMANOVA: time explains 6.3% of the distance between biopsies (p < 0.001); group, on one mean profile per person, 6.6% (p = 0.424).
- The Training_Interaction effect detectable at 80% power is 0.95 log2 units for the median protein.

## Method

- `fit_limpa_contrasts()` (`00_Setup/Shared_functions/analysis_reporting_helpers.R`): `limpa::dpcDE()` with person as a random block and limpa's standard errors as weights, then `contrasts.fit()` and `eBayes(robust = TRUE)`.
- Each logFC with its 95% interval (`topTable(confint = TRUE)`), BH within contrast, Storey's q.
- Checks: the same fit without the blood term (`unadjusted`), null calibration, the effect detectable at 80% power.
- `fry` tests one group's p < 0.05 proteins in the other group's contrast (`fry_concordance`).
- PCA and PERMANOVA (`vegan::adonis2`, 9,999 permutations, time permuted within person).

## Read with care

- Seven or eight people per group: read logFC with its interval, and nominal p as a lead.
- A protein that changes in one group and not the other is not, on that basis, a group difference; the interactions test that.

## Sheets

`protein_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `protein_by_contrast` | One row per contrast and protein: logFC with its 95% CI, t, p, BH and Storey q, and the samples observed. |
| `contrast_summary` | Tested, p < 0.05 and BH counts per contrast, with its role. |
| `null_calibration` | propTrueNull, median \|t\| and share of p below 0.2 per contrast. |
| `detectable_effect` | log2 effect the primary contrast detects at 80% power, median and 90th percentile. |
| `fry_concordance` | fry: one arm's responders tested in the other arm's contrast, and training sets on each arm's acute contrast; with mixed p, mean log2 FC in both contrasts and sign kept. |
| `set_trajectory` | Training p < 0.05 sets: mean change from each subject's T1 (and SE) per arm and time point. |
| `unadjusted` | logFC and p per protein and contrast before blood adjustment, for sensitivity. |
| `pca` | PC1 and PC2 per sample, from the complete blood-corrected matrix with proteins scaled. |
| `pca_variance` | Share of variance on the first four components. |
| `permanova` | PERMANOVA (9,999 permutations): proteome time per window and arm, arm per time point; phenotype time. |
