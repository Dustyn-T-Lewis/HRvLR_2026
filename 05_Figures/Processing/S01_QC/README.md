# S01_QC

Did the samples hold up, and does the analysis behave? Ten panels: the top row checks the data (filtering, completeness, normalisation, dynamic range, the detection curve), the bottom row the statistics (within-person correlation, blood, p-value calibration, power). The script tests nothing; it reads `03_Group_comparisons/HRvLR/c_data/02_protein_contrasts.xlsx` as well, so `run_all.R` runs it last.

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

| | |
|---|---|
| Reads | `01_Processing/01_Protein_Filtering/c_data/protein_filtering.xlsx` (`filter_log`, `outlier_diagnostics`, `blood_by_arm`), `01_Processing/02_Normalisation/c_data/normalisation.xlsx` (`intensity_quantiles`), `01_Processing/03_Quantification/c_data/quantification.xlsx` (`curve`, `proteins`, `blood`, `sample_correlation`), `03_Group_comparisons/HRvLR/c_data/02_protein_contrasts.xlsx` (`protein_by_contrast`, `unadjusted`, `null_calibration`, `detectable_effect`) |
| Writes | `b_reports/s01_qc.pdf` (300 x 190 mm, cairo PDF), `c_data/s01_qc.xlsx` |
| Run | `Rscript 05_Figures/Processing/S01_QC/a_script/s01_qc.R`, or `run_all.R` with the figures at the end (it reads the HR/LR protein contrasts) |

## Reading the panels

- Calibration (i): with no change, p-values spread evenly and the median |t| sits near 0.67; a pile of small p on a flat floor is signal.
- Power (j): the smallest `Training_Interaction` log2 FC each protein could show at 80% power, against the observed effects.
- Panel e is the case for limpa: missingness falls as abundance rises, so a missing value says "faint".
- The script stops unless every plotted count traces to its table: 48 runs, three excluded, 45 kept, nine contrasts, two power quantiles.

## Sheets

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
