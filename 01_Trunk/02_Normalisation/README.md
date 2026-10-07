# 02_Normalize

Compares proteoDA's normalisation methods and applies cyclic loess.

| | |
|---|---|
| Reads | `01_Filter/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/02_normalize.xlsx`, `b_reports/02_normalize_figures.pdf`, and proteoDA's `norm_comparison.pdf`, `qc_pre.pdf`, `qc_post.pdf` in `b_reports/` |
| Run | `Rscript 01_Trunk/02_Normalize/a_script/02_normalize.R` |

## Method

`normalize_data("cycloess")`, which calls `limma::normalizeCyclicLoess(method = "fast")`. With
`adaptive.span = TRUE` the span comes from `chooseLowessSpan(nrow)`, so it moves whenever filtering
changes the protein count. proteoDA writes its method comparison and QC reports to `b_reports/`, and
its log warnings concern those plots, not the normalised matrix.

Run medians vary little (SD 0.23 log2) and cyclic loess scores best of proteoDA's comparison, so
the choice barely moves results. Missing values stay missing; `03_Quantify` models them.

## Figures

One page: PCA of the samples after normalisation, gaps median-filled for display only.

## Workbook

| Sheet | Holds |
|---|---|
| `normalized` | Cyclic-loess log2 abundance, one row per protein. |
| `components` | Variance explained by the first four principal components. |
| `pca_scores` | Each sample's PC1 and PC2 score. |
| `eta_squared` | Share of each protein's variance explained by the six group cells. |
| `intensity_quantiles` | Each sample's log2 intensity quantiles (5, 25, 50, 75, 95%) before and after cyclic loess. |
