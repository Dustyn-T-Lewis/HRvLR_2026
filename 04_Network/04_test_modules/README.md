# 04_test_modules

The modules on the nine contrasts, by eigengene and by fry, and each module's membership against
protein-level significance.

| | |
|---|---|
| Reads | `modules.rds`, `design.rds`, `fit.rds`, `phenotype.rds`, `DAList_imputed.rds` |
| Writes | `module_tests.rds`, `04_test_modules.xlsx`, 5 figures |

The eigengene test is `lmFit()` with the protein design and subject block, correlation re-estimated on the
eigengenes (0.095), `eBayes(robust = TRUE)`, BH within contrast. 5 of 108 tests are nominal
against 5.4 expected, none at BH < 0.05.

fry tests each module as a protein set on the imputed matrix, with the same design and block and
a correlation estimated on that matrix (0.176). fry's FDR runs over the twelve modules within each contrast. 7 of
108 are nominal; red rising in Acute_LR is the one at FDR < 0.05 (0.045).

Membership against significance (WGCNA's MM against GS) is, within each module, the Spearman
correlation between a member's kME and its moderated t per contrast, and its rho with each phenotype over training. On
the primary contrast blue (rho 0.32) and green (0.32) run positive and pink negative (−0.32).
Members share a module, so these correlations are descriptive, not tests.

Figures: eigengene contrasts, fry contrasts, membership against each contrast and each phenotype,
and kME against moderated t for Training_Interaction in every module.
