# 01_Filtering

Removes contaminant proteins by identity, drops outlier samples, then filters on detection.

| | |
|---|---|
| Reads | `00_Input/HRvLR_raw.xlsx`, `metadata.csv`, `blood_contaminants.csv`, `HPA_annotations_full.tsv`, `RBC_proteome_reference.tsv` |
| Writes | `c_data/DAList_filtered.rds`, `c_data/01_filter.xlsx`, `b_reports/01_filter_figures.pdf` |
| Run | `Rscript 01_Preprocess/01_Filtering/a_script/01_filter.R` |
| Cost | about 35 s |

The figure PDF is one page: the filtering cascade, each sample's contamination share by panel
before removal, and the outlier flags for every flagged sample. `blood_cor_null` holds the
permutation null that retired the old `blood_cor` cut.
