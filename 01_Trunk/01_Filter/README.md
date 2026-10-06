# 01_Filter

Removes contaminant proteins by identity, drops outlier samples, then filters on detection.

| | |
|---|---|
| Reads | `00_Input/HRvLR_raw.xlsx`, `metadata.csv`, `blood_contaminants.csv`, `HPA_annotations_full.tsv`, `RBC_proteome_reference.tsv` |
| Writes | `c_data/DAList_filtered.rds`, `c_data/01_filter.xlsx`, `b_reports/01_filter_figures.pdf` |
| Run | `Rscript 01_Trunk/01_Filter/a_script/01_filter.R` |

## Method

Contaminants go by identity: the curated list in `00_Input/blood_contaminants.csv`, plus HPA plasma,
immunoglobulin and erythrocyte tags. A muscle rescue (myonuclei at least 20 nCPM and plasma
concentration below 1e9 pg/L) overturns an HPA tag. `blood_cor` is reported beside every call and
gates nothing. Each sample's `blood_rel` (mean log2 of HBB and HBA1 minus the sample's median
log2 intensity) goes into the metadata; every later model carries it as a covariate. The two adult chains are seen
in all 48 runs. An earlier index also averaged HBG2, which was missing in the 20 low-blood runs, so
averaging what was present lifted exactly those runs. A sample goes when three of four outlier methods flag it. The detection filter runs
last, so no discarded sample counts: a protein needs 5 detections in at least one group cell.

## Figures

One page: the filtering cascade, each sample's contamination share by panel before removal, and the
outlier flags.

## Workbook

| Sheet | Holds |
|---|---|
| `filter_log` | Proteins left after each step. |
| `contaminants_removed` | Every protein removed as a contaminant, with its verdict and reason. |
| `rescued` | Blood-tagged proteins kept by the muscle rescue. |
| `protein_calls` | One verdict per protein, with the annotation behind it. |
| `contamination_index` | Share of each sample's signal per contamination panel, before removal. |
| `blood_index` | Mean log2 haemoglobin per sample, raw and relative to the sample's median log2 intensity. |
| `outlier_diagnostics` | The four outlier flags per sample. |
| `blood_by_arm` | Arm by T3 blood term on two scales, with a permutation p. |
| `blood_cor_null` | 99.9th percentile of the blood_cor permutation null, by observation count. |
