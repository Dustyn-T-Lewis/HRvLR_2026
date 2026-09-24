# 01_run_fgsea_and_fry

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `gene_sets.rds`, `fit.rds`, `design.rds`, `DAList_imputed.rds` |
| Writes | `set_tests.rds`, `set_tests.csv`, `01_run_fgsea_and_fry.xlsx`, dot plots, hit matrices |

`topTable()` rebuilds the nine contrasts from the saved fit, keeping `02_Differential`'s BH.
fgsea ranks proteins by moderated t, seeded; a protein untested in a contrast has no t and leaves
that ranking. `collapsePathways` re-tests each significant set against a stronger set's leading
edge; survivors carry `main = TRUE`. fry reads the imputed matrix with the fit's design, the
subject block and a within-subject correlation estimated on that matrix (0.176).

The floor contrast prints first. `set_summary` gives, per contrast, the sets each method tested
(fry all 1,378, fgsea those with at least 15 measured proteins in that ranking) and the counts called.

## Outputs

`set_tests` has one row per set, contrast and method. `NES`, `leadingEdge` and `main` are
fgsea-only. `leadingEdge` holds gene symbols; it is a list column in the RDS and `;`-joined in the
CSV.

`b_reports/` has one folder per collection plus `all_db/`, one dot plot per drawn contrast (the
ten strongest collapse survivors), and `hits/03_set_hits.pdf`: every set nominal under fry in
at least one contrast, as a dot matrix filled by fgsea NES, 75 sets to a page.
