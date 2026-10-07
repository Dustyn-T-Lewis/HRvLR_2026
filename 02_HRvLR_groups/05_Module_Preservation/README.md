# 05_Module_Preservation

Builds modules inside each arm and tests their preservation in the other arm.

| | |
|---|---|
| Reads | `01_Trunk/03_Quantification/c_data/DAList_corrected.rds`, `01_Trunk/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/module_preservation.xlsx`, `b_reports/06_preserve_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/05_Module_Preservation/a_script/02_preserve.R` |

## Method

Each arm gets its own network with the full-cohort settings. `modulePreservation()` tests each arm's
modules in the other arm, 200 permutations. Zsummary above 10 is strong preservation, 2 to 10 weak
to moderate, below 2 none (Langfelder et al. 2011). Each arm module is labelled with its
best-overlapping full-cohort module.

Each arm network rests on 20 to 24 samples, so the test has little power to show that a module is
missing from the other arm. A Zsummary above 2 says a module is found again; a Zsummary below 2 is
weak evidence that it is absent; neither shows the two networks are the same.

## Figures

One page: Zsummary and median rank against module size, for each direction.

## Workbook

| Sheet | Holds |
|---|---|
| `preservation` | Zsummary and median rank of each arm's modules tested in the other arm. |
| `arm_networks` | Samples, soft power and module count of each arm's network. |
| `arm_membership` | Every protein's module in each arm's network. |
