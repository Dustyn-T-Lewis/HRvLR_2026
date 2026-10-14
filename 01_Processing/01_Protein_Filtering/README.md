# 01_Protein_Filtering

Which proteins and which runs are fit to analyse?

|  |  |
|------------------------------------|------------------------------------|
| Reads | `00_Setup/Input/HRvLR_raw.xlsx`, `00_Setup/Input/metadata.csv`, `01_Processing/01_Protein_Filtering/references/blood_contaminants.csv`, `01_Processing/01_Protein_Filtering/references/HPA_annotations_full.tsv`, `01_Processing/01_Protein_Filtering/references/RBC_proteome_reference.tsv` |
| Writes | `c_data/DAList_filtered.rds`, `c_data/protein_filtering.xlsx`, `b_reports/protein_filtering_report.pdf` |
| Run | `Rscript 01_Processing/01_Protein_Filtering/a_script/protein_filtering.R` |

## Key result

-   2,400 proteins in the raw report; 222 removed as contaminants, 278 by the detection filter; 1,900 kept.
-   3 of 48 runs fail three or more of the four checks and are dropped: S29_T1, S28_T2, S28_T3.

## Method

-   Contaminants go by identity, never by behaviour: the 30 accessions written into the script (keratins, globins, immunoglobulins, complement), the 95 in `01_Processing/01_Protein_Filtering/references/blood_contaminants.csv`, and proteins the Human Protein Atlas tags as blood unless their myonuclear RNA marks them as muscle.
-   A run goes when three of four checks flag it: missingness, Mahalanobis distance on the first three principal components, median intensity (Hampel, 3 MAD) and median correlation with other runs.
-   A protein stays when it is detected in at least 5 samples of at least one of the six group-by-time cells (`proteoDA::filter_proteins_by_group()`).
-   The blood index is the mean log2 intensity of HBB and HBA1 minus the sample's median log2 intensity (`blood_rel`). Later steps use it as a covariate or regress it out.

## Read with care

-   The detection filter is the one step in 01_Processing that reads the HR/LR label. It counts detections, never abundance.
-   `blood_cor` and `in_rbc` are reported and decide nothing; `blood_cor_null` shows why no single correlation cut fits every protein.

## Sheets

`protein_filtering.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|------------------------------------|------------------------------------|
| `filter_log` | Proteins left after each step. |
| `contaminants_removed` | Every protein removed as a contaminant, with its verdict and reason. |
| `rescued` | Blood-tagged proteins kept by the muscle rescue. |
| `protein_calls` | One verdict per protein, with the annotation behind it. |
| `contamination_index` | Share of each sample's signal per contamination panel, before removal. |
| `blood_index` | Mean log2 haemoglobin per sample, raw and relative to the sample's median log2 intensity. |
| `outlier_diagnostics` | The four outlier flags per sample. |
| `blood_by_arm` | Arm by T3 blood term on two scales, with a permutation p. |
| `blood_cor_null` | 99.9th percentile of the blood_cor permutation null, by observation count. |