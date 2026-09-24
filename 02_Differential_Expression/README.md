# 02 · Differential Expression

Fits the normalised DAList, tests each protein on nine contrasts and eight classification tasks,
and correlates each protein with ten phenotypes.

| Step | Runs | Writes |
|---|---|---|
| [`01_Design`](01_Design/README.md) | design formula, nine contrasts, within-subject correlation | `design.rds` |
| [`02_Differential`](02_Differential/README.md) | proteoDA fit, results, null checks, classification, hit matrices | `fit.rds` |
| [`03_Phenotype`](03_Phenotype/README.md) | every protein against the ten phenotypes in three windows | `phenotype.rds` |

## Six cell means, subject random

```r
add_design(dal, "~ 0 + group + (1 | subject)")
```

One cell mean per arm and timepoint: 45 samples, 39 residual degrees of freedom. Subject is a
random effect because HR and LR are different people; a fixed subject term would absorb every
between-arm contrast. The within-subject correlation is 0.189. The design must be full rank and
every column name must survive `make.names()`, because `add_contrasts()` parses the contrast
strings as R code; the notebook asserts both.

| Contrast | Role | Definition |
|---|---|---|
| `Training_HR` | descriptive | HR_T2 − HR_T1 |
| `Training_LR` | descriptive | LR_T2 − LR_T1 |
| `Acute_HR` | descriptive | HR_T3 − HR_T2 |
| `Acute_LR` | descriptive | LR_T3 − LR_T2 |
| `Baseline_HRvLR` | floor | HR_T1 − LR_T1 |
| `Trained_HRvLR` | descriptive | HR_T2 − LR_T2 |
| `Acute_HRvLR` | descriptive | HR_T3 − LR_T3 |
| `Training_Interaction` | primary | (HR_T2 − HR_T1) − (LR_T2 − LR_T1) |
| `Acute_Interaction` | secondary | (HR_T3 − HR_T2) − (LR_T3 − LR_T2) |

`Training_Interaction` equals `Trained_HRvLR` minus `Baseline_HRvLR` exactly. The floor is not a
negative control: the HR/LR label was cut from the training outcome, so a real baseline difference
would be predictive.

## No contrast has a BH hit

`fit_limma_model()` runs `duplicateCorrelation()` by subject, `lmFit()` with that block,
`contrasts.fit()` and `eBayes(robust = TRUE)` on the unimputed matrix. BH runs within each
contrast. A protein with an empty cell in a contrast is not tested there, so the denominators run
from 1,877 to 1,892.

The lowest FDR across all nine is 0.071 (Acute_HR). The primary contrast
could detect a median effect of 1.04 log2 units at 80% power.

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

A calibrated null sits at 1, 0.674 and 0.200. Five contrasts depart from it on all three:
Training_HR, Acute_HR, Acute_LR, Acute_HRvLR and Acute_Interaction. `treat()` at 1.15-fold returns
nothing in any contrast.

The eight classification tasks are training (T1 to T2) and the acute bout (T2 to T3) within each
arm, paired, and HR against LR on the T1 level, the T2 level, the training change and the acute
change. AUC is the rank-sum statistic over the product of the group sizes; p is the Wilcoxon test,
signed-rank when paired. No task has a BH hit.

The `fry` appendix tests the proteins that responded in one arm as a set in the other arm's
contrast, on the imputed matrix. Two of sixteen tests reach FDR < 0.05.

## Two protein-outcome pairs survive BH

Each subject contributes one value per protein in three windows: training change (T2 − T1, the
window the phenotypes share), T1 level and acute change (T3 − T2). Spearman with the t
approximation, BH within window and outcome. A protein is tested only when two thirds of the
window's subjects observed it. Above nine observations `cor.test`'s "exact" Spearman p is an
Edgeworth series that returns 0 in the far tail; it put RPS4X at FDR 0.

| Window | Outcome | Nominal / chance | BH < 0.05 |
|---|---|---:|---:|
| training | d_1rm_ext | 1.96 | 0 |
| training | d_mcsa | 0.95 | 1 (RPS4X, rho 0.91, n 14) |
| baseline | d_1rm_ext | 2.36 | 0 |
| baseline | d_1rm_legpress | 1.25 | 1 (ACSL3, rho −0.87, n 15) |
| baseline | volume_load | 1.46 | 0 |
| acute | comp_hypertrophy | 1.86 | 0 |
| acute | d_mcsa | 1.52 | 0 |

The other 23 cells sit between 0.44 and 1.18. The same correlations inside each arm rest on 4 to
8 subjects; the workbook keeps those with p < 0.05.
