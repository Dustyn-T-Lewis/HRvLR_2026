# S01_Volcano

Supplement: which proteins move in each within-group and interaction contrast, and which fgsea sets do they belong to? Each volcano plots logFC against -log10 p, colours proteins at BH FDR < 0.05, and rings up to eight non-redundant fgsea sets at FDR < 0.05, each an arc coloured by NES with ticks to its leading-edge proteins.

| Page | Contrast | Labels |
|---|---|---|
| 1 | `Training_Interaction` | up to five proteins at FDR < 0.05, smallest FDR first |
| 2 | `Acute_Interaction` | as above |
| 3 to 6 | `Training_HR`, `Training_LR`, `Acute_HR`, `Acute_LR` | as above |
| 7, 8 | `Training_HR`, `Training_LR` again | the five proteins with the smallest Π |

| | |
|---|---|
| Reads | `03_Group_comparisons/HRvLR/03_Set_Contrasts/c_data/set_tests.rds` (`protein_results`, `set_tests`) |
| Writes | `b_reports/s01_volcano.pdf`, `c_data/s01_volcano.xlsx` |
| Run | `Rscript 05_Figures/HRvLR/S01_Volcano/a_script/s01_volcano.R` |

## Read with care

- A protein's colour (protein FDR) and a set's arc (set FDR) are separate tests; a protein under a ringed set need not be significant itself.
- Pages 7 and 8 use Π (p^|logFC|) to choose labels only; a protein named there is not a hit.
- One group's training volcano beside the other's does not compare the groups; the interaction pages do.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `ringed_sets` | a set ringed on one volcano | `contrast`, `rank_by`, `set_id`, `collection`, `pathway`, `nes`, `fdr` |
| `labelled_proteins` | a protein named on one volcano | `contrast`, `rank_by` (fdr or pi), `label` |
