# 02_Normalisation

How are the 45 samples made comparable before any protein is compared across them?

| | |
|---|---|
| Reads | `01_Processing/01_Protein_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/normalisation.xlsx`, `b_reports/norm_comparison.pdf`, `b_reports/normalisation_report.pdf`, `b_reports/qc_post.pdf`, `b_reports/qc_pre.pdf` |
| Run | `Rscript 01_Processing/02_Normalisation/a_script/normalisation.R` |

## Key result

- 1,900 proteins normalised; the six group-by-time cells explain a median 12% of a protein's variance (eta squared).
- PC1 and PC2 carry 11.6% and 10.2% of the variance.

## Method

- log2, then cyclic loess (`proteoDA::normalize_data(norm_method = "cycloess")`). Missing values stay missing; `03_Quantification` models them.
- proteoDA's standard reports compare every normalisation method (`norm_comparison.pdf`) and show sample QC before (`qc_pre.pdf`) and after (`qc_post.pdf`).
- `eta_squared` is the share of each protein's variance the six group-by-time cells explain (`limma::lmFit(~ group)`); descriptive, read by nothing downstream.

## Read with care

- The PCA fills gaps with the sample median for display only; the saved matrix keeps them missing.

## Sheets

`normalisation.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `normalized` | Cyclic-loess log2 abundance, one row per protein. |
| `components` | Variance explained by the first four principal components. |
| `pca_scores` | Each sample's PC1 and PC2 score. |
| `eta_squared` | Share of each protein's variance explained by the six group cells. |
| `intensity_quantiles` | Each sample's log2 intensity quantiles (5, 25, 50, 75, 95%) before and after cyclic loess. |
