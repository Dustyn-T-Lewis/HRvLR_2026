# 01_Filtering

Removes contaminant proteins by identity, drops outlier samples, then filters on detection.

| | |
|---|---|
| Reads | `00_Input/HRvLR_raw.xlsx`, `HRvLR_meta.csv`, `blood_contaminants.csv`, `HPA_annotations_full.tsv`, `RBC_proteome_reference.tsv` |
| Writes | `c_data/DAList_filtered.rds`, `c_data/01_filter.xlsx`, `b_reports/01_filter_figures.pdf` |
| Run | `quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports` |
| Cost | about 35 s |
