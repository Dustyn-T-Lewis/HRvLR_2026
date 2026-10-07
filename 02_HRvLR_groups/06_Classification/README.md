# 06_Classification

Can a single protein, gene-set score or module eigengene tell one group from the other?

## The idea

The contrast stages ask whether a group's average moved. This stage asks a different question
about people: if you knew only one feature's value for a person, how well could you say which
group they are in? Three substages ask it at three levels and share one function,
`classify_step()` in `R/classify.R`.

| Substage | Feature | Value per biopsy |
|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | a protein | log2 abundance, blood-corrected complete matrix |
| [`02_Sets`](02_Sets/README.md) | a gene set | singscore |
| [`03_Modules`](03_Modules/README.md) | a module | eigengene |

Every feature is scored on five tasks. Each compares the two groups (HR and LR in this arc); none
pairs a person with themselves, since within-group change is the subject of the contrast stages.

| Task | Compares |
|---|---|
| `Baseline_HRvLR` | the groups' values at T1 |
| `Trained_HRvLR` | the groups' values at T2 |
| `Acute_HRvLR` | the groups' values at T3 |
| `Training_change_HRvLR` | each person's T2 - T1 change, people with both biopsies |
| `Acute_change_HRvLR` | each person's T3 - T2 change, people with both biopsies |

The AUC (area under the ROC curve) measures separation. Pick one person from each group at random:
the AUC is the chance that the first group's person has the higher value, ties counting half. 0.5
is no separation, 1 means every member of the first group sits above every member of the second, 0
means every one sits below. It equals the Mann-Whitney count W (the number of cross-group pairs
in which the first group is higher) divided by n1 x n2, the number of such pairs. The ROC curve
behind the name slides a cut-off along the values and plots, at each cut, the share of the first
group above it (sensitivity) against the share of the second group above it (1 - specificity);
the area under that curve is the same number.

Nothing is fitted. The Wilcoxon rank-sum test gives the p for AUC = 0.5, and BH FDR and Storey's q
adjust it over the features in a task. Because the test uses ranks and takes no covariate, the
values it reads already have the blood slope removed.

An interval is the range of AUC the data allow. Two are given:

- DeLong's 95% interval (`pROC::ci.auc(method = "delong")`) on every feature. It comes from a
  formula for the AUC's standard error, which is approximate with seven or eight people per group.
  At perfect separation (AUC 0 or 1) its variance is zero, so no interval is given.
- A stratified bootstrap 95% interval (`pROC::ci.auc(method = "bootstrap")`, 2,000 resamples,
  seeded) on the features drawn in the report. It resamples people within each group, keeping the
  group sizes, recomputes the AUC each time and takes the middle 95%. It is slow, so only drawn
  features get it.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Tasks | `by_subject()`, `change()` inside `classify_step()` | One column per person at a time point, or the change over people sampled at both; split into first and second group. |
| AUC and p | `classify_row()` in `R/helpers.R`: `wilcox.test()`, `pROC::roc(direction = "<")`, `pROC::ci.auc(method = "delong")` | A feature needs three people per side. `direction = "<"` pins the first group as the higher side; pROC would otherwise flip any AUC below 0.5. |
| Adjustment | `p.adjust("BH")`, `storey_q()` | Within task (and within collection for sets). |
| Counts | `dplyr::summarise()` | Per task: n per group, features tested, at p < 0.05 and at BH < 0.05, and `min_p`, the smallest p the rank-sum test can reach, 2 / choose(n1 + n2, n1). |
| Hits | `paginate()` in `R/helpers.R` | Features at p < 0.05 per task, ordered by p, then distance of AUC from 0.5; at most 24 per task, 12 to a page. |
| Bootstrap | `pROC::ci.auc(method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)` | `set.seed(1)`; drawn features only. |
| Report | `ComplexHeatmap`, `pROC::ggroc()` | Overview heatmaps, then ROC curves. |

## How to read the outputs

- Read the share of features at p < 0.05 against the 5% that chance gives with no group
  difference at all. It summarises the screen; it does not measure how well the groups separate.
  Features are correlated, so the share wanders further from 5% than independent tests would.
- With seven or eight people per group the rank-sum p is discrete and has a floor (`min_p`). With
  ties, `wilcox.test()` switches to its normal approximation.
- The best AUC among hundreds of features is high even when nothing separates the groups. A
  drawn feature is one of many; read its interval, not its point AUC.
- A perfectly separating feature gives a bootstrap interval of [1, 1] (or [0, 0]): every resample
  separates perfectly too. That says nothing about how a new sample of people would split.
- `Baseline_HRvLR` carries the same caution as the floor contrast: HR and LR come from an outcome
  of training, so a T1 difference shows how far the split separates people before training.

## Files

Each substage writes the same three things:

```
0x_Level/
├── a_script/<level>_classification.R
├── b_reports/<level>_classification_report.pdf
└── c_data/<level>_classification.xlsx
```

Each runs on its own with `Rscript`; none reads another's output. `Figures/F3_People` reads the
`auc` and `counts` sheets of all three.

## Workbook: `<level>_classification.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a task (and collection, for sets) | `n_positive`, `n_negative` (people per group), `min_p`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `auc` | a feature in a task | `auc`, `auc_lower`, `auc_upper` (DeLong), `p`, `fdr`, `qvalue`, `boot_lower`, `boot_upper` (drawn features only), `n_positive`, `n_negative` |
| `plot_index` | an overview page or a drawn feature | `task`, `feature`, `pdf_page`, `position` (1 to 12) |

## Report: `<level>_classification_report.pdf`

1. Overview: AUC for every feature in the five tasks (0 to 1, white at 0.5) beside a p/FDR status
   map.
2. ROC curves for the features at p < 0.05 in each task, twelve to a page, at most 24 per task.
   Each header gives the AUC with its bootstrap interval, p, BH q and the feature's n per group;
   the subtitle gives the task's tested, p < 0.05 and BH < 0.05 counts.
