# 06_Preserve

Builds modules inside each arm and tests their preservation in the other arm.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/06_preserve.xlsx`, `b_reports/06_preserve_figures.pdf` |
| Run | `Rscript 04_Network/06_Preserve/a_script/06_preserve.R` |

## Method

Each arm gets its own network with the full-cohort settings. `modulePreservation()` tests each arm's
modules in the other arm, 200 permutations. Zsummary above 10 is strong preservation, 2 to 10
moderate, below 2 none (Langfelder et al. 2011). Each arm module is labelled with its
best-overlapping full-cohort module.

## Figures

One page: Zsummary and median rank against module size, per direction.

## Workbook

| Sheet | Holds |
|---|---|
| `preservation` | Zsummary and median rank of each arm's modules tested in the other arm. |
| `arm_networks` | Samples, soft power and module count of each arm's network. |
| `arm_membership` | Every protein's module in each arm's network. |
