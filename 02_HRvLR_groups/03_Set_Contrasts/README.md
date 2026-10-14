# 03_Set_Contrasts

Which gene sets shift in each of the nine contrasts?

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds`, `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/set_contrasts.xlsx`, `c_data/set_tests.rds`, `b_reports/set_contrasts_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/03_Set_Contrasts/a_script/set_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_set_contrasts.R`, shared with the cluster groupings |

## Key result

- 1,378 gene sets tested per contrast.
- fgsea calls 0 (Trained_HRvLR) to 299 (Acute_HR) sets at FDR < 0.05; the baseline contrast (the floor) calls 70.
- fry calls 0 sets at FDR < 0.05 across the nine contrasts.

## Method

- fgsea (`fgsea::fgsea()`, seeded) ranks each set against the other proteins on the protein moderated t: does the set sit high or low in the ranking?
- fry (`limma::fry()`) tests whether the set's proteins move consistently across people, on limpa's complete matrix with person blocked.
- Each method's BH runs within contrast; `fgsea::collapsePathways()` marks the non-redundant fgsea calls (`main`).

## Read with care

- Read the baseline (floor) count first: it is the number of calls the split alone produces.
- fgsea treats proteins as exchangeable, which they are not, so its calls run far ahead of fry's. An fgsea call without a fry call is a ranking result, not evidence that people moved consistently.

## Sheets

`set_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast by method, and fgsea collapse survivors. |
| `significant` | Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on limpa's matrix. |
