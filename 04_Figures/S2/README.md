# S2: single-feature classification

How well one protein, set score or module eigengene tells HR from LR, from
`02_HRvLR/06_Classification`. No model is fitted and no AUC is measured on held-out people.

| Panel | Shows |
|---|---|
| a `a_auc_top` | Per task (T1, T2, T3 levels; T2 − T1 and T3 − T2 changes) and feature type, the three features with the smallest Wilcoxon p (ties broken by AUC furthest from 0.5): AUC with DeLong's 95% interval, none at AUC 0 or 1. Dashed line: 0.5, chance. Shade: which group is higher. Shape: filled at FDR < 0.05, open at p < 0.05 only, a cross otherwise. Under each cell: features at p < 0.05 against the number chance gives. |

Set names come from `pathway_label()` (`style/pathway_grid.R`); modules are named by their leading
set with the colour after it (`style/module_labels.R`).

## Inputs

| File | Used for |
|---|---|
| `02_HRvLR/06_Classification/c_data/classification.rds` | `auc` and `counts` |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | module names |

## Run

```sh
Rscript 04_Figures/S2/a_script/S2_composite.R
Rscript 04_Figures/S2/a_script/panels/a_auc_top.R
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/a_auc_top.pdf` | the panel alone |
| `b_reports/S2_composite.pdf` | the figure, 180 mm wide |
| `c_data/S2.xlsx` | the numbers drawn, column dictionary, input md5 |
