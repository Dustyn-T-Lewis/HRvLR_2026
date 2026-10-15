# 03_Module_Contrasts

Tests every co-expression module in the nine HR-against-LR contrasts two ways: limma on the
module's eigengene, and fry on the module's proteins as a set. A third table describes, within each
module, whether members that follow the module closely (high kME) are the ones that changed most.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes (module × biopsy), membership with kME |
| `02_HRvLR/01_Protein_Contrasts/c_data/fit.rds` | the protein fit (`fit$EList`, correlation), moderated contrasts, design, contrast matrix |

## Run

```sh
quarto render 02_HRvLR/03_Module_Contrasts/a_script/module_contrasts.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/module_tests.rds` | `eigengene_tests`, `fry_tests`, `membership_contrast`, for the figures |
| `c_data/module_contrasts.xlsx` | the same three tables |
| `b_reports/module_contrasts.html` | the rendered walkthrough |

## Read with care

- 13 modules × 9 contrasts, about 5.9 tests per method at p < 0.05 by chance. Eigengene tests: 4
  at p < 0.05, none at FDR < 0.05. fry: 12 at p < 0.05 and 4 at FDR < 0.05, all within-group
  changes (tan up with training in LR; blue, green and pink down after the bout in LR). No HR − LR
  or interaction test reaches FDR < 0.05 by either method.
- The eigengene fit uses the protein design and person block, with the within-person correlation
  estimated on the eigengenes; limpa's sample weights are not applied to eigengenes.
- The kME-against-t correlations and their p are descriptions: members were chosen for moving
  together.
- Module count and colours come from `06_Module_Building` and change if its soft power changes.
