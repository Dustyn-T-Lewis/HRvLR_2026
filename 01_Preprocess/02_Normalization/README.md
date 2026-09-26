# 02_Normalization

Compares proteoDA's normalisation methods and applies cyclic loess.

| | |
|---|---|
| Reads | `01_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/02_normalize.xlsx`, `b_reports/02_normalize_figures.pdf`, and proteoDA's `norm_comparison.pdf`, `qc_pre.pdf`, `qc_post.pdf` in `b_reports/` |
| Run | `Rscript 01_Preprocess/02_Normalization/a_script/02_normalize.R` |
| Cost | about 25 s |

The figure PDF is one page: PCA of the samples after normalisation, gaps median-filled for display
only. proteoDA writes its own warnings to the log (MD plots use the first two groups; partial NA
coefficients for 34 proteins); both concern its QC plots, not the normalised matrix.
