# 03_Protein_Contrasts

Which proteins change between time points, and which differ between the groups?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`02_HRvLR_groups/02_Protein_Contrasts`](../../../02_HRvLR_groups/02_Protein_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_protein_contrasts.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/03_Quantification/c_data/proteins.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `<input>/c_data/fit.rds`, `<input>/c_data/protein_contrasts.xlsx`, `<input>/b_reports/protein_contrasts_report.pdf`, `c_data/protein_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/03_Protein_Contrasts/a_script/protein_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_protein_contrasts.R` |

## Key result

- `T1`
  - 1,900 proteins tested in each of the nine contrasts; 130 at FDR < 0.05 across all nine (lowest FDR 0.000).
  - Proteins at p < 0.05 run from 73 (Acute_Interaction) to 408 (Baseline_C1vC2); chance alone gives about 95 per contrast.
  - PERMANOVA: time explains 6.7% of the distance between biopsies (p < 0.001); group, on one mean profile per person, 16.2% (p < 0.001).
  - The Training_Interaction effect detectable at 80% power is 0.99 log2 units for the median protein.
- `T2`
  - 1,900 proteins tested in each of the nine contrasts; 4 at FDR < 0.05 across all nine (lowest FDR 0.015).
  - Proteins at p < 0.05 run from 52 (Acute_C1vC2) to 268 (Training_C2); chance alone gives about 95 per contrast.
  - PERMANOVA: time explains 6.6% of the distance between biopsies (p = 0.001); group, on one mean profile per person, 6.4% (p = 0.640).
  - The Training_Interaction effect detectable at 80% power is 0.93 log2 units for the median protein.
- `T3`
  - 1,900 proteins tested in each of the nine contrasts; 271 at FDR < 0.05 across all nine (lowest FDR 0.000).
  - Proteins at p < 0.05 run from 47 (Training_Interaction) to 493 (Acute_C1vC2); chance alone gives about 95 per contrast.
  - PERMANOVA: time explains 6.6% of the distance between biopsies (p = 0.001); group, on one mean profile per person, 12.4% (p = 0.009).
  - The Training_Interaction effect detectable at 80% power is 1.03 log2 units for the median protein.
- `training`
  - 1,900 proteins tested in each of the nine contrasts; 162 at FDR < 0.05 across all nine (lowest FDR 0.000).
  - Proteins at p < 0.05 run from 46 (Acute_C1vC2) to 407 (Training_C1); chance alone gives about 95 per contrast.
  - PERMANOVA: time explains 7.0% of the distance between biopsies (p < 0.001); group, on one mean profile per person, 6.4% (p = 0.772).
  - The Training_Interaction effect detectable at 80% power is 0.96 log2 units for the median protein.
- `acute`
  - 1,900 proteins tested in each of the nine contrasts; 180 at FDR < 0.05 across all nine (lowest FDR 0.000).
  - Proteins at p < 0.05 run from 71 (Baseline_C1vC2) to 408 (Acute_C1); chance alone gives about 95 per contrast.
  - PERMANOVA: time explains 6.6% of the distance between biopsies (p = 0.001); group, on one mean profile per person, 5.9% (p = 0.788).
  - The Training_Interaction effect detectable at 80% power is 1.03 log2 units for the median protein.

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
