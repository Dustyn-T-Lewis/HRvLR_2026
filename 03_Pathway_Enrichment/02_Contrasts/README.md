# 02_Contrasts

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `00_Gene_Sets/c_data/gene_sets.rds`, `01_Scores/c_data/singscore.rds`, `00_Input/phenotype.csv`, `02_Differential_Expression/02_Contrasts/c_data/fit.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/set_tests.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/02_Contrasts/a_script/02_contrasts.R` |

## Method

fgsea ranks proteins by moderated t from the protein fit and asks whether a set sits at one end
(seeded). fry asks whether the set moved at all under the design, subject block and within-subject
correlation, on the imputed matrix. Each method's BH runs within contrast over the five collections
pooled. `collapsePathways()` marks non-redundant fgsea hits in `main` and deletes nothing.
`set_tests.rds` holds `set_tests` (one row per set, contrast and method) and `protein_results`.

## Figures

Two overview pages per collection: singscore across the 45 biopsies, centred per set, with arm, time
point, type I share, fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the
biopsy traits are blank at T3, where none were measured); then fgsea NES in the nine contrasts
beside a status map. Then a dot plot per contrast of every fgsea collapse survivor at FDR < 0.05,
placed by NES, and one per contrast of every set at fry p < 0.05, placed by signed -log10 p, forty
sets to a page.

## Workbook

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast, by method, and fgsea survivors of collapsePathways. |
| `significant` | Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on the imputed matrix. |
