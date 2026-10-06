# S04_QC

Supplement to Figure 1: whether the samples held up, and whether the analysis behaves.

| | |
|---|---|
| Reads | `01_Trunk/01_Filter/c_data/01_filter.xlsx` (`filter_log`, `outlier_diagnostics`, `blood_by_arm`), `01_Trunk/02_Normalize/c_data/02_normalize.xlsx` (`intensity_quantiles`), `01_Trunk/03_Quantify/c_data/03_quantify.xlsx` (`curve`, `proteins`, `blood`, `sample_correlation`), `02_HRvLR_groups/02_DEP/c_data/02_contrasts.xlsx` (`DEP_matrix`, `unadjusted`, `null_calibration`, `detectable_effect`) |
| Writes | `b_reports/s04_qc.pdf`, `b_reports/s04_qc.png` (300 x 190 mm, 600 dpi), `c_data/s04_qc.xlsx` |
| Run | `Rscript 01_Trunk/Figures/S01_QC/a_script/s01_qc.R` |

The script reads saved outputs and tests nothing; its own summaries (ordering, within- and
between-person correlation, histograms) describe them.

## Panels

| Panel | Shows | Reading |
|---|---|---|
| a | Proteins left after each filtering step | 2,400 to 1,900: contaminant classes, then detection |
| b | Share of proteins missing in each of the 48 runs | kept runs miss 12-39%; the three excluded runs miss 44% (S29_T1), 57% (S28_T2) and 68% (S28_T3) |
| c | Each sample's log2 intensity quantiles before and after cyclic loess | medians aligned after normalisation |
| d | Proteins ranked by mean abundance, with muscle landmarks named | 1,900 proteins over 19.1 log2 units; myoglobin, creatine kinase, actin and myosins at the top |
| e | Share of values missing against mean observed abundance, with limpa's detection curve | 12.3% missing, concentrated in faint proteins (curve slope 0.67), so values are missing because they are low, which limpa models |
| f | Correlation between the 45 samples, each protein centred, clustered | median r 0.17 between two samples from one person and -0.03 between different people |
| g | Blood index per biopsy | LR - HR at T3 0.84 log2, permutation p = 0.23 |
| h | Each protein's log2 FC with and without the blood covariate, nine contrasts pooled | r = 0.895; the covariate trims effects without reordering them |
| i | p-value histograms per contrast, against a flat null | median \|t\| 0.57 to 0.88 (0.67 under the null); Acute HR shows the clearest excess near zero |
| j | Observed \|log2 FC\| in the training interaction against the effect it detects at 80% power | 99% of observed effects fall below the median detectable 0.95 log2 |

Panels i and j explain the null single-protein result: the tests are calibrated, and with seven or
eight people per arm the primary contrast detects a 0.95 log2 change for the median protein, a
size almost no observed effect reaches.

## Workbook

| Sheet | Holds |
|---|---|
| `filter_log` | Panel a: proteins removed and remaining at each filtering step. |
| `detection` | Panel b: share of proteins missing per run, and whether the run was kept. |
| `intensity_quantiles` | Panel c: per-sample log2 intensity quantiles before and after cyclic loess. |
| `dynamic_range` | Panel d: proteins ranked by mean log2 abundance over the complete matrix. |
| `missingness` | Panel e: per protein, mean observed abundance and share of values observed. |
| `correlation_summary` | Panel f: median sample correlation within and between people. |
| `blood` | Panel g: blood index per sample. |
| `blood_adjustment` | Panel h: correlation of log2 FC with and without the blood covariate. |
| `calibration` | Panel i: per contrast, propTrueNull, median \|t\| and share of p below 0.2. |
| `power` | Panel j: log2 FC the training interaction detects at 80% power. |
