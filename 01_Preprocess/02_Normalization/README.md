# 02_Normalization

Compares proteoDA's normalisation methods, applies cyclic loess, then removes each protein's blood
slope.

| | |
|---|---|
| Reads | `01_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/DAList_unadjusted.rds`, `c_data/02_normalize.xlsx`, `b_reports/02_normalize_figures.pdf`, and proteoDA's `norm_comparison.pdf`, `qc_pre.pdf`, `qc_post.pdf` in `b_reports/` |
| Run | `Rscript 01_Preprocess/02_Normalization/a_script/02_normalize.R` |

## Method

`normalize_data("cycloess")`, which calls `limma::normalizeCyclicLoess(method = "fast")`. With
`adaptive.span = TRUE` the span comes from `chooseLowessSpan(nrow)`, so it moves whenever filtering
changes the protein count. proteoDA writes its method comparison and QC reports to `b_reports/`, and
its log warnings concern those plots, not the normalised matrix.

LR's T3 biopsies carry more blood than any other cell (`01_filter.xlsx`, `blood_by_arm`). Each
protein is fitted on the six group means plus the centred `blood_rel`, subject blocked
(`duplicateCorrelation()`), and its blood slope times the sample's blood value is subtracted. The
group means enter the fit, so the slope comes from spread within cells and no contrast leaks into
it; a protein with no estimable slope stays as it was. `DAList_normalized.rds` holds the adjusted
matrix that every later step reads; `DAList_unadjusted.rds` keeps the matrix before adjustment for
the sensitivity sheet in `02_Differential_Expression/02_Contrasts`.

## Figures

One page: PCA of the samples after normalisation, gaps median-filled for display only.

## Workbook

| Sheet | Holds |
|---|---|
| `normalized` | Cyclic-loess log2 abundance with the blood slope removed, one row per protein. |
| `components` | Variance explained by the first four principal components. |
| `pca_scores` | Each sample's PC1 and PC2 score. |
| `eta_squared` | Share of each protein's variance explained by the six group cells. |
| `blood_slope` | log2 change per unit of the centred blood index, with its moderated p, per protein. |
