# 01_Protein_Filtering

Which proteins and which runs are fit to analyse?

## The idea

The mass spectrometer reports every protein it identifies in the sample, and a needle biopsy is
not pure muscle. Skin keratins come in on the needle and from dust; blood comes in with the tissue
and brings haemoglobin, plasma proteins, immunoglobulins and complement. These proteins say how
much skin or blood a biopsy caught, not what the muscle did, so they go before anything is
measured.

They go by identity, never by behaviour. A protein that rises and falls with haemoglobin across
samples might be blood, or might be a muscle protein that happens to track the blood load; a
correlation cut cannot tell the two apart. So a protein goes when a curated list names it, or when
the Human Protein Atlas (HPA) tags it as plasma, immunoglobulin or red-cell and nothing marks it
as muscle. The "muscle rescue" keeps an HPA-tagged protein whose RNA is clearly present in
muscle nuclei (myonuclei) and whose blood concentration is not high.

Some runs fail outright: too little material, too few proteins identified. Four independent checks
look at each run, and a run goes only when three of the four flag it. The detection filter then
drops proteins seen too rarely to model, and it runs last so a discarded run cannot count as a
detection.

Blood that stays in the tissue is measured rather than removed. The blood index is each sample's
mean log2 intensity of the two adult haemoglobin chains (HBB, HBA1), which every run detects.
`blood_rel` is that index minus the sample's median log2 intensity, so the amount of protein loaded
does not enter it. It rides along in the sample sheet; later steps use it as a model covariate or
regress it out.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Read | `readxl::read_excel()`, `readr::read_csv()` | The protein report and the sample sheet; the script stops if their run names differ. |
| Deduplicate | `dplyr::slice_max()` | One row per accession: the row with the highest mean intensity. |
| Blood index | `log2()`, `colMeans()` | Mean log2 HBB and HBA1 per sample, taken before the globins are removed. `blood_rel` goes into the sample sheet. `blood_cor`, each protein's Spearman correlation with the index (proteins seen in 15 or more samples), is reported and gates nothing. |
| Contamination share | `colSums()` | Each sample's share of summed intensity in five marker panels (keratin, blood, leukocyte, adipose, myofibre), before removal. Descriptive. |
| Contaminant calls | `dplyr::case_when()` | Removed: the 30 accessions written into the script (18 keratins and other skin proteins, most from cRAP's dust and contact section; 7 globins; 3 immunoglobulins; 2 complement proteins), matched on accession only, and the 95 in `00_Input/blood_contaminants.csv`; or an HPA tag (secreted to blood, immunoglobulin gene class, erythrocyte RNA at 5,000 nCPM or more) without rescue. Rescue: myonuclei RNA at 20 nCPM or more and HPA blood concentration below 1e9 pg/L or not given. Curated entries cannot be rescued. Membership of the red-cell reference is reported as `in_rbc` and removes nothing. |
| Build | `proteoDA::DAList()`, `proteoDA::zero_to_missing()` | The kept proteins as a DAList (intensities, annotation, sample sheet); zeros become missing. |
| Outlier runs | `prcomp()`, `mahalanobis()`, `mad()`, `cor()` | Four flags: missingness above the Tukey fence (Q3 + 1.5 IQR) on its share or on its spread within the subject; Mahalanobis distance on the first three principal components beyond the chi-square 0.99 quantile; median log2 intensity more than 3 MAD from the median (Hampel); median correlation with other runs more than 3 MAD below the median. Three of four removes the run (`proteoDA::filter_samples()`). |
| Detection | `proteoDA::filter_proteins_by_group()` | A protein stays when it is detected in at least 5 samples of at least one of the six arm-by-time cells (`group`). |
| Blood by arm | `lme4::lmer()`, `permutation_p()` in `R/helpers.R` | `arm * timepoint + (1 \| subject)` on the blood index and on the blood panel share; the arm-by-T3 term, with a p from 999 shuffles of arm across subjects. A check, not a correction. |
| `blood_cor` null | `cor(method = "spearman")` | 300 shuffles of the blood index; the 99.9th percentile of \|rho\| by how many samples saw the protein (15-25, 26-40, 41-48). Records why a single correlation cut was dropped. |

## How to read the outputs

- `filter_log` counts proteins after each step. `n_samples` drops from 48 to 45 at the outlier
  step; the detection count uses the 45.
- `verdict` gives the rule that decided each protein; `reason` gives the curated reason where one
  exists. A protein marked `keep: rescued (muscle-expressed)` carried an HPA blood tag and was kept.
- `in_rbc` and `blood_cor` inform and decide nothing.
- The detection filter is the one step in the trunk that reads the HR/LR label, through the six
  `group` cells. It counts detections, never abundance.
- `blood_by_arm` asks whether one arm carries more blood at T3. A blood difference between arms
  would not cancel in an interaction contrast, which is why every later model carries the index.
- Keeping the three failed runs was tried on 2026-10-06 and reversed: they flattened limpa's
  detection curve in `03_Quantification`, and limma's sample weights did not single them out.

## Files

| | |
|---|---|
| Reads | `00_Input/HRvLR_raw.xlsx`, `00_Input/metadata.csv`, `00_Input/blood_contaminants.csv`, `00_Input/HPA_annotations_full.tsv`, `00_Input/RBC_proteome_reference.tsv` |
| Writes | `c_data/DAList_filtered.rds` (unnormalised intensities, annotation, sample sheet with `blood_rel`; read by `02_Normalisation`), `c_data/protein_filtering.xlsx`, `b_reports/protein_filtering_report.pdf` |
| Run | `Rscript 01_Trunk/01_Protein_Filtering/a_script/protein_filtering.R` |

```
01_Protein_Filtering/
├── a_script/protein_filtering.R
├── b_reports/protein_filtering_report.pdf
└── c_data/
    ├── protein_filtering.xlsx
    └── DAList_filtered.rds
```

## Workbook: `protein_filtering.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `filter_log` | a filtering step | `step`, `n_removed`, `n_after`, `pct_of_raw`, `n_samples` |
| `contaminants_removed` | a removed protein | `verdict`, `reason`, `secretome`, `blood_conc`, `ery`, `myo` |
| `rescued` | a blood-tagged protein kept by the rescue | `myo`, `blood_conc` |
| `protein_calls` | a protein in the report | `verdict`, `reason`, `contaminant`, `in_rbc`, `blood_cor` |
| `contamination_index` | a sample and marker panel | `panel`, `pct_signal` |
| `blood_index` | a sample | `blood_index`, `blood_rel` |
| `outlier_diagnostics` | a run | `pct_missing`, `delta_missing`, `miss_flag`, `pca_flag`, `mad_flag`, `cor_flag`, `n_flags`, `consensus_outlier` |
| `blood_by_arm` | a scale (log2 index, panel share) | `lr_minus_hr_at_t3`, `t`, `permutation_p` |
| `blood_cor_null` | an observation-count bin | `observations`, `proteins`, `null_999` |

## Report: `protein_filtering_report.pdf`

One page, three panels:

1. The filtering cascade: proteins left after each step, with the number each step removed.
2. Each sample's share of signal in the keratin, blood, leukocyte and adipose panels before
   removal, by time point, HR filled and LR open.
3. The outlier flags for every run flagged at least once.
