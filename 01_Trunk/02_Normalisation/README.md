# 02_Normalisation

How do we make the 45 samples comparable before any protein is compared across them?

## The idea

Two runs of the same sample never give identical intensities. One biopsy yields more protein than
another, the instrument drifts over a batch, and a run can respond differently to faint and bright
peptides. These shifts touch every protein in a sample at once, so they look like biology without
being biology.

Normalisation removes them. It rests on one assumption: most proteins do not change between
samples, so whatever moves the whole distribution of a sample is technical. After it, a difference
left in one protein is a difference in that protein.

Cyclic loess (limma's `normalizeCyclicLoess`, fast method) goes further than shifting each sample's
median. On the log2 scale it plots, for every protein, the gap between a sample and the average
of all samples against the protein's overall intensity, fits a smooth curve (loess) through that
cloud, and subtracts the curve. A sample that reads high only for bright proteins is corrected
for bright proteins only. The span of the smoother comes from `chooseLowessSpan()`, so it moves
when the number of proteins changes.

Missing values stay missing. A missing value in this experiment usually means the protein was too
faint to detect, and `03_Quantification` models that. Filling it here with a guess would hide the
information.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Compare methods | `proteoDA::write_norm_report()` | Every normalisation method proteoDA offers, drawn side by side with `group` as the grouping, to `norm_comparison.pdf`. |
| QC before | `proteoDA::write_qc_report()` | Sample-level QC plots coloured by `group`, to `qc_pre.pdf`. |
| Normalise | `proteoDA::normalize_data(norm_method = "cycloess")` | log2, then cyclic loess. The before-values for the QC table are logged the same way. |
| QC after | `proteoDA::write_qc_report()` | The same plots after normalisation, to `qc_post.pdf`. |
| Group share | `limma::lmFit(~ group)` | eta squared: the share of each protein's variance the six arm-by-time cells explain. Descriptive; proteins seen in fewer than 4 samples get none. |
| PCA | `prcomp(scale. = TRUE)` | Gaps filled with the sample median for this display only; the saved matrix keeps them missing. |
| Quantiles | `quantile()` | Each sample's 5th, 25th, 50th, 75th and 95th percentile of log2 intensity, before and after. |

## How to read the outputs

- After normalisation the sample medians and spreads in `intensity_quantiles` should line up. A
  sample that stays offset after cyclic loess differs in more than loading.
- The PCA fills gaps with medians, which pulls faint proteins toward the middle. Read it as a
  picture of sample structure, not as data.
- `eta_squared` uses the HR/LR cells and is descriptive. Nothing downstream reads it.
- proteoDA's three PDFs are its standard reports. `norm_comparison.pdf` and `qc_pre.pdf` show the
  data before cyclic loess, `qc_post.pdf` after; none of them feeds a later step.

## Files

| | |
|---|---|
| Reads | `01_Trunk/01_Protein_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds` (log2 cyclic-loess DAList, gaps kept; read by `03_Quantification` and `04_Gene_Sets`), `c_data/normalisation.xlsx`, `b_reports/normalisation_report.pdf`, and proteoDA's `b_reports/norm_comparison.pdf`, `b_reports/qc_pre.pdf`, `b_reports/qc_post.pdf` |
| Run | `Rscript 01_Trunk/02_Normalisation/a_script/normalisation.R` |

```
02_Normalisation/
├── a_script/normalisation.R
├── b_reports/
│   ├── normalisation_report.pdf
│   ├── norm_comparison.pdf
│   ├── qc_pre.pdf
│   └── qc_post.pdf
└── c_data/
    ├── normalisation.xlsx
    └── DAList_normalized.rds
```

## Workbook: `normalisation.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `normalized` | a protein | `uniprot_id`, `gene`, then log2 abundance per sample |
| `components` | a principal component (PC1 to PC4) | `variance_explained` |
| `pca_scores` | a sample | `PC1`, `PC2`, `subject`, `arm`, `timepoint` |
| `eta_squared` | a protein | `eta2` |
| `intensity_quantiles` | a sample before or after | `stage`, `q05`, `q25`, `q50`, `q75`, `q95` |

## Report: `normalisation_report.pdf`

One page: PCA of the samples after normalisation, coloured by time point, shaped by arm.
