# 05_Module_Preservation

Does each group's co-expression structure reappear in the other group?

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/06_Module_Building/c_data/modules.rds` |
| Writes | `c_data/module_preservation.xlsx`, `b_reports/module_preservation_report.pdf` |
| Run | `Rscript 03_Group_comparisons/HRvLR/07_Module_Preservation/a_script/module_preservation.R` |

## Key result

- HR modules tested in LR: 5 strongly preserved, 13 weak to moderate, 1 not preserved.
- LR modules tested in HR: 4 strongly preserved, 10 weak to moderate, 1 not preserved.

## Method

- Modules built inside each group with the 01_Processing settings (`WGCNA::blockwiseModules()`).
- `WGCNA::modulePreservation()` in both directions, 200 permutations; Zsummary below 2 is none, 2 to 10 weak to moderate, above 10 strong.
- Each group module is named by the full-cohort module it overlaps most.

## Read with care

- Preservation compares correlation structure, not abundance: a module can be preserved while its level differs between the groups.

## Sheets

`module_preservation.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `preservation` | Zsummary and median rank of each arm's modules tested in the other arm. |
| `arm_networks` | Samples, soft power and module count of each arm's network. |
| `arm_membership` | Every protein's module in each arm's network. |
