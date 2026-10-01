# 06_Preserve

Builds modules inside each arm and tests their preservation in the other arm.

| | |
|---|---|
| Reads | `01_Preprocess/03_Quantification/c_data/DAList_corrected.rds`, `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/06_preserve.xlsx`, `b_reports/06_preserve_figures.pdf` |
| Run | `Rscript 04_Network/06_Preserve/a_script/06_preserve.R` |

## Method

Each arm gets its own network with the full-cohort settings. `modulePreservation()` tests each arm's
modules in the other arm, 200 permutations. Zsummary above 10 is strong preservation, 2 to 10 weak
to moderate, below 2 none (Langfelder et al. 2011). Each arm module is labelled with its
best-overlapping full-cohort module.

Each arm network rests on 20 to 24 samples, where even one population split in two scores
modestly. The calibration splits the 16 subjects into random halves of 8, arm ignored, 50 times,
and runs the same test (50 permutations each): the Zsummary two halves of one population reach at
this size. Arm modules inside that range are as preserved as the cohort is with itself.

## Figures

Two pages: Zsummary and median rank against module size, per direction; then the arm modules'
Zsummary beside the random-halves calibration.

## Workbook

| Sheet | Holds |
|---|---|
| `preservation` | Zsummary and median rank of each arm's modules tested in the other arm. |
| `calibration` | Zsummary of each module across random halves of all subjects, the reference at this size. |
| `arm_networks` | Samples, soft power and module count of each arm's network. |
| `arm_membership` | Every protein's module in each arm's network. |
