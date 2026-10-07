# 06_Classification/03_Modules

Can one module eigengene tell one group from the other?

## The idea

The method, the five tasks, the AUC and its two intervals are explained in
[`06_Classification`](../README.md). This substage runs them on every module eigengene.

A module is a group of proteins that rise and fall together across biopsies (WGCNA, in
`01_Trunk/06_Module_Building`). Its eigengene is the first principal component of its members:
one number per biopsy that summarises the module, signed to follow the members' average. The
eigengenes are scored on blood-corrected abundance, so they carry differences between people;
grey, the proteins in no module, has none.

There are few modules, so BH has few tests to adjust over and Storey's q equals BH
(`storey_q()` returns BH below 100 tests).

## What the code does

| Step | Function | Detail |
|---|---|---|
| Values | `readRDS()` | `eigengenes` and `meta` from `modules.rds`. The script stops if the columns are out of order. |
| Labels | `set_names(rownames(values))` | Each module is named by its colour. |
| Row order | `module_summary$module` | Module size, largest first. |
| Scoring | `classify_step()` in `R/classify.R` | As described in the parent README. |

## How to read the outputs

- With few modules, the 5% chance share is a fraction of one module per task. Read each module's
  AUC and interval directly rather than the count.
- An eigengene summarises a module's dominant pattern. A module can separate the groups through a
  few members without its eigengene doing so; `01_Proteins` would show those members.

## Files

| | |
|---|---|
| Reads | `01_Trunk/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/module_classification.xlsx` (read by `Figures/F3_People`), `b_reports/module_classification_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/03_Modules/a_script/module_classification.R` |

```
03_Modules/
├── a_script/module_classification.R
├── b_reports/module_classification_report.pdf
└── c_data/module_classification.xlsx
```

## Workbook: `module_classification.xlsx`

Sheets as in the parent README.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a task | `n_positive`, `n_negative`, `min_p`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `auc` | a module in a task | `auc`, `auc_lower`, `auc_upper`, `p`, `fdr`, `qvalue`, `boot_lower`, `boot_upper` |
| `plot_index` | an overview page or a drawn module | `task`, `feature`, `pdf_page`, `position` |

## Report: `module_classification_report.pdf`

1. One overview page: AUC for every module in the five tasks beside a p/FDR status map, rows by
   size.
2. ROC curves for the modules at p < 0.05 in each task.
