# 02_Differential_Expression / 02_Differential

Fits the model, writes the results, and classifies every protein on eight tasks.

| | |
|---|---|
| Script | `a_script/02_differential.qmd` |
| Reads | `01_Design/c_data/design.rds`; for the appendices, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` and, when present, V1's `05_results.xlsx` |
| Writes | `c_data/fit.rds`, `c_data/differential.rds`, `c_data/02_differential.xlsx`, `b_reports/02_differential_figures.pdf` |

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

## The sequence

```r
fit <- fit_limma_model(dal)
res <- extract_DA_results(fit, pval_thresh = 0.05, lfc_thresh = 0, adj_method = "BH")
write_limma_plots(res, grouping_column = "group", ...)
```

Inside `fit_limma_model()`: `duplicateCorrelation()` by subject, `lmFit()` with that block,
`contrasts.fit()`, `eBayes(robust = TRUE)`. The matrix is unimputed, so each protein is fitted on
the samples where it was seen. `write_limma_plots()` knits its own R Markdown and cannot run inside
a Quarto render, so it runs in a fresh session through `callr`. Its per-contrast HTML reports are
not tracked.

## What comes out

`02_differential.xlsx`:

| Sheet | Holds |
|---|---|
| `DEP_matrix` | one row per protein: `logFC`, `P.Value`, `adj.P.Val` per contrast, and samples observed |
| `contrast_summary` | tested, nominal, pi-score, BH counts per contrast |
| `null_calibration` | `propTrueNull`, median \|t\|, share of p < 0.2 |
| `treat` | BH < 0.05 against a 1.15-fold floor |
| `pi_ranking` | top ten by pi-score in the four within-arm contrasts |
| `protein_auc` | AUC, Wilcoxon p and BH per protein and task |
| `chance_expectation` | nominal hits against chance per task, with the set and module levels' columns and labels |
| `v1_equivalence` | maximum difference from V1's committed results |
| `fry_concordance` | one arm's training signature tested in the other |

`fit.rds` is the full proteoDA result, including the limma fit. `differential.rds` holds every
sheet above except `DEP_matrix` and `v1_equivalence`, with provenance.

## Multiple testing

BH within each contrast, never pooled. A protein with an empty cell in a contrast is not tested
there, so each contrast has its own denominator (1,877 to 1,892).

## Reading a null contrast

| Contrast | prop. true null | median \|t\| | p < 0.2 |
|---|---:|---:|---:|
| Training_Interaction *(primary)* | 1.00 | 0.64 | 0.16 |
| Acute_Interaction *(secondary)* | 0.96 | 0.73 | 0.24 |
| Baseline_HRvLR *(floor)* | 1.00 | 0.61 | 0.16 |
| Training_HR | 0.90 | 0.79 | 0.25 |
| Training_LR | 1.00 | 0.64 | 0.18 |
| Acute_HR | 0.85 | 0.85 | 0.30 |
| Acute_LR | 0.98 | 0.69 | 0.24 |
| Trained_HRvLR | 1.00 | 0.63 | 0.17 |
| Acute_HRvLR | 0.94 | 0.76 | 0.24 |

Under a calibrated null these sit at 1, 0.674 and 0.200. Five contrasts sit above the null on all
three: Training_HR, Acute_HR, Acute_LR, Acute_HRvLR and Acute_Interaction. None reaches a BH hit.
The rest, the floor and the primary contrast among them, sit at or below it. `treat()` at 1.15-fold returns
nothing anywhere.

## Classification

Eight tasks: training (T1→T2) and acute (T2→T3) within each arm, paired; and HR against LR on the
T1 level, the T2 level, the training change and the acute change. AUC from `pROC` with
`direction = "<"`, p from the Wilcoxon test. No task has a BH hit.

## Hit matrices

The figure bundle holds the p-value histogram, then every protein nominal in any contrast and
every protein nominal in any task, as dot matrices 75 rows to a page.

## Appendices

The V1 equivalence appendix reproduces V1's nine contrasts to at most 5.5e-12. It is skipped
when the sibling V1 tree is absent.

The fry concordance appendix tests the proteins that responded in one arm as a set in the other
arm's contrast. Two of sixteen set-ranking tests are concordant at FDR < 0.05.
