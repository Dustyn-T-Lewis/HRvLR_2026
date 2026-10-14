# 02_Normalisation

Makes the 45 runs comparable with cyclic loess on log2 intensities. Missing values stay missing;
`03_Quantification` models them.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/01_Protein_Filtering/c_data/DAList_filtered.rds` | filtered raw intensities, 1,900 proteins by 45 runs |

## Run

```sh
quarto render 01_Processing/02_Normalisation/a_script/normalisation.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/normalisation.html` | the walkthrough: before and after, eta squared, PCA |
| `b_reports/norm_comparison.pdf` | proteoDA's comparison of every normalisation method |
| `b_reports/qc_pre.pdf`, `b_reports/qc_post.pdf` | proteoDA's run-level QC before and after |
| `c_data/DAList_normalized.rds` | normalised log2 abundance; read by `03_Quantification` and `04_Gene_Sets` |
| `c_data/normalisation.xlsx` | normalised matrix, PCA, eta squared, intensity quantiles; read by `05_Figures/Processing/S01_QC` |
