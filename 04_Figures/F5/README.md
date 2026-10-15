# F5: modules

Co-expression modules from `01_Processing/06_Module_Building` (WGCNA, soft power 12). Colours are
labels only, so every module is named by its top over-represented set (FDR < 0.05) with the colour
after it (`style/module_labels.R`). Rows run largest module first in a and b.

| Panel | Shows |
|---|---|
| a_identity | Each module's size, its name, and STRING enrichment of its members (observed over expected edges at score ≥ 700; filled at FDR < 0.05). |
| b_eigengene_matrix | limma on each module's eigengene in the nine contrasts, in the F2 grammar: colour logFC, size -log10 p, filled at p < 0.05, black ring at FDR < 0.05. A grey square marks fry p < 0.05 for the module's proteins as a set. The subtitle sets the p < 0.05 count against chance (4 of 117 tests, 5.9 expected). |
| c_trajectories | Each person's eigengene at T1, T2 and T3, HR and LR, with the group means, for the module with the smallest eigengene p in a training, acute or HR - LR contrast (interactions left out). Exploratory: the script picks the module by p, and the subtitle gives its contrast, p and FDR (0.27). |
| d_preservation | `modulePreservation()` Zsummary for each group's own modules tested in the other group; dashed lines at 2 and 10, shading below 2. Each row names the group module and the full-cohort module holding most of its proteins. Point size is module size. |

## Inputs

| File | Used for |
|---|---|
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes, biopsy table, module sizes |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | `labels` (top set), `string` |
| `02_HRvLR/03_Module_Contrasts/c_data/module_tests.rds` | `eigengene_tests`, `fry_tests` |
| `02_HRvLR/05_Module_Preservation/c_data/preservation.rds` | `preservation`, `group_networks` |

## Run

```sh
Rscript 04_Figures/F5/a_script/F5_composite.R                 # figure and workbook
Rscript 04_Figures/F5/a_script/panels/b_eigengene_matrix.R    # one panel alone
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel |
| `b_reports/F5_composite.pdf` | the figure, 180 × 165 mm |
| `c_data/F5.xlsx` | one sheet per panel with the numbers drawn |
