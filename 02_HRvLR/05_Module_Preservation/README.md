# 05_Module_Preservation

Builds WGCNA modules inside HR and inside LR with the network settings of `06_Module_Building`,
then tests each group's modules in the other group's data with `modulePreservation()` (Zsummary and
median rank, 200 permutations). The soft power is the WGCNA FAQ value for a signed
network of 24 biopsies, 16; the full network in `06_Module_Building` uses 10, chosen by
leave-one-person-out stability.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | limpa's complete matrix with blood slopes removed |
| `01_Processing/06_Module_Building/c_data/modules.rds` | full-cohort modules, used to label each group module |

## Run

```sh
quarto render 02_HRvLR/05_Module_Preservation/a_script/module_preservation.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/preservation.rds` | `preservation` (Zsummary, median rank, best full-cohort match), `group_networks`, `group_membership` |
| `c_data/module_preservation.xlsx` | the same tables |
| `b_reports/module_preservation.html` | the rendered walkthrough |

## Read with care

- HR's network has 13 modules (273 of 2,058 proteins unassigned); LR's has 16 (508 unassigned).
  One protein has no spread inside a group and is dropped from both. Each network rests on 24
  biopsies from 8 people.
- Zsummary runs from 1.3 to 9.8. Three of 29 group modules fall below 2 (one HR, two LR); none
  passes 10.
- Every group module maps best to a full-cohort module, not to the unassigned `grey` proteins,
  holding 22% to 77% of its proteins in that module (`full_share`).
- With 8 people per group the test has little power to show a module is missing.
