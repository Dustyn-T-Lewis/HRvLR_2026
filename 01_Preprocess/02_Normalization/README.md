# 02_Normalization

Compares proteoDA's normalisation methods and applies cyclic loess.

| | |
|---|---|
| Reads | `01_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/02_normalize.xlsx`, `b_reports/02_normalize_figures.pdf`, and proteoDA's `norm_comparison.pdf`, `qc_pre.pdf`, `qc_post.pdf` in `b_reports/` |
| Run | `quarto render 01_Preprocess/02_Normalization/a_script/02_normalize.qmd --output-dir ../b_reports` |
| Cost | about 25 s |
