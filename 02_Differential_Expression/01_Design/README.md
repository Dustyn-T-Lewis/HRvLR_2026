# 02_Differential_Expression / 01_Design

Attaches the design formula and the nine contrasts to the DAList, and measures the correlation
the random effect carries. Separate from the fit so a design problem surfaces in seconds.

| | |
|---|---|
| Script | `a_script/01_design.qmd` |
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/design.rds`, `c_data/01_design.xlsx` |

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports
```

## The design

```r
add_design(dal, "~ 0 + group + (1 | subject)")
```

45 samples by 6 cell means, one per responder-and-timepoint combination, 39 residual degrees of
freedom. proteoDA reads the formula's terms from metadata column names, so the notebook relabels
the sample sheet first (`group`, `subject`, `responder`, `time`).

Subject is random because HR and LR are different people: a fixed subject term would absorb the
between-arm contrasts. `(1 | subject)` makes proteoDA estimate one consensus within-subject
correlation and pass it to `lmFit()` as a block.

## Two assertions before anything is fitted

The design must be full rank, or the contrasts come back as silent NAs. Every column name must
survive `make.names()`, because `add_contrasts()` parses the contrast strings as R code.

## What design.rds holds

The DAList with design and contrasts attached, the nine contrast strings, the role table and the
name of the floor contrast. `02_Differential` fits from it directly.

## The correlation diagnostic

`duplicateCorrelation()` on the normalised matrix, blocked by subject: 0.189. `fit_limma_model()`
estimates the same number and discards it, so the workbook's `correlation_strata` sheet is where
it is recorded.
