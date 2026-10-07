# S01_QC

Did the samples hold up, and does the analysis behave?

It reads a workbook from `02_HRvLR_groups`, so `run_all.R` runs it with the figures at the end.

## The idea

A supplementary figure in ten panels. The top row asks whether the data are sound: how many
proteins survived filtering, how complete each run is, whether normalisation lined the samples up,
how wide the dynamic range is, and whether missing values follow intensity as limpa assumes. The
bottom row asks whether the statistics behave: whether biopsies from one person resemble each
other, how much blood each biopsy carries and what the blood covariate does, whether p-values are
calibrated, and what effect size the primary contrast can detect.

Two of the checks need a premise:

- Calibration: if no protein changed, p-values would spread evenly between 0 and 1 and the
  median \|t\| would sit near 0.67. A pile of small p-values on a flat floor is signal; a histogram
  bent everywhere suggests a mis-specified model.
- Power: for each protein the moderated standard error from limma fixes the smallest log2 fold
  change the training interaction could detect at 80% power, two-sided 0.05. Comparing the observed
  effects with that size says whether a null result means "no effect" or "too small to see here".

The script tests nothing itself. Its own summaries (ordering, within- and between-person
correlation, histograms) describe the saved outputs it reads.

## What the code does

| Panel | Shows | Reads |
|---|---|---|
| a | Proteins left after each filtering step | `protein_filtering.xlsx`, `filter_log` |
| b | Share of proteins missing in each of the 48 runs, the three dropped runs marked | `protein_filtering.xlsx`, `outlier_diagnostics`; kept runs are those in `quantification.xlsx`, `blood` |
| c | Each sample's log2 intensity quantiles (median, 25-75%, 5-95%) before and after cyclic loess | `normalisation.xlsx`, `intensity_quantiles` |
| d | Proteins ranked by mean abundance over limpa's complete matrix, with muscle landmarks named (MYH7, MYH2, ACTA1, TTN, DES, CKM, MB, ATP5F1A, NDUFA4, PSMA1) | `quantification.xlsx`, `proteins` |
| e | Share of values missing against mean observed abundance, with limpa's detection curve drawn as 1 minus detection probability | `quantification.xlsx`, `proteins`, `curve` |
| f | Correlation between the 45 samples (each protein centred, blood-corrected), clustered by average linkage; subtitle gives the median within and between people | `quantification.xlsx`, `sample_correlation`, `blood` |
| g | Blood index per biopsy by time point and arm, with the arm-by-T3 estimate and permutation p | `quantification.xlsx`, `blood`; `protein_filtering.xlsx`, `blood_by_arm` |
| h | Each protein's log2 FC with and without the blood covariate, nine contrasts pooled, with their correlation | `protein_contrasts.xlsx`, `protein_by_contrast`, `unadjusted` |
| i | p-value histogram per contrast against a flat null, with median \|t\| | `protein_contrasts.xlsx`, `protein_by_contrast`, `null_calibration` |
| j | Observed \|log2 FC\| in the training interaction against the detectable effect for the median protein and the 90th percentile | `protein_contrasts.xlsx`, `protein_by_contrast`, `detectable_effect` |

Before writing, `stopifnot()` checks that every plotted count traces to its table: 48 runs, three
excluded, 45 kept, nine contrasts, two power quantiles.

## How to read the outputs

- Panel b's three dropped runs are the consensus outliers from `01_Protein_Filtering`.
- Panel e is the visual case for limpa: if missingness falls as abundance rises, a missing value
  says "faint", which is what the detection curve models.
- Panel f is drawn on a capped scale (r from -0.6 to 0.6); the diagonal is blank.
- Panel h shows how far the blood covariate moves effects. Points on the diagonal are untouched.
- Panel j answers the power question for one contrast, `Training_Interaction`, only.

## Files

| | |
|---|---|
| Reads | `01_Trunk/01_Protein_Filtering/c_data/protein_filtering.xlsx` (`filter_log`, `outlier_diagnostics`, `blood_by_arm`), `01_Trunk/02_Normalisation/c_data/normalisation.xlsx` (`intensity_quantiles`), `01_Trunk/03_Quantification/c_data/quantification.xlsx` (`curve`, `proteins`, `blood`, `sample_correlation`), `02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`protein_by_contrast`, `unadjusted`, `null_calibration`, `detectable_effect`) |
| Writes | `b_reports/s01_qc.pdf` (300 x 190 mm, cairo PDF), `c_data/s01_qc.xlsx` |
| Run | `Rscript 01_Trunk/Figures/S01_QC/a_script/s01_qc.R`, or `run_all.R` with the figures at the end (it reads the HR/LR protein contrasts) |


```
S01_QC/
├── a_script/s01_qc.R
├── b_reports/s01_qc.pdf
└── c_data/s01_qc.xlsx
```

## Workbook: `s01_qc.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | Panel | One row is | Key columns |
|---|---|---|---|
| `filter_log` | a | a filtering step | `step`, `n_removed`, `n_after` |
| `detection` | b | a run | `pct_missing`, `kept` |
| `intensity_quantiles` | c | a sample before or after | `stage`, `q05` to `q95` |
| `dynamic_range` | d | a protein | `rank`, `mean_complete` |
| `missingness` | e | a protein | `mean_observed`, `prop_observed` |
| `correlation_summary` | f | same person or not | `same`, `r` (median) |
| `blood` | g | a sample | `blood_rel` |
| `blood_adjustment` | h | the pooled comparison | `r`, `n_points` |
| `calibration` | i | a contrast | `prop_true_null`, `median_abs_t`, `p_below_0.2` |
| `power` | j | a quantile | `log2_fc` |
