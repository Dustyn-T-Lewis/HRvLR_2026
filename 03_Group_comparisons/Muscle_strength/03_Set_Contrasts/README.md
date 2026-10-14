# 04_Set_Contrasts

Which gene sets shift in each of the nine contrasts?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`03_Group_comparisons/HRvLR/03_Set_Contrasts`](../../../03_Group_comparisons/HRvLR/03_Set_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_set_contrasts.R`.

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `03_Group_comparisons/Muscle_strength/02_Protein_Contrasts/c_data/fit.rds`, `03_Group_comparisons/Muscle_strength/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/set_contrasts.xlsx`, `c_data/set_tests.rds`, `b_reports/set_contrasts_report.pdf` |
| Run | `Rscript 03_Group_comparisons/Muscle_strength/03_Set_Contrasts/a_script/set_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_set_contrasts.R` |

## Key result

- 1,378 gene sets tested per contrast.
- fgsea calls 57 (Training_Interaction) to 285 (Baseline_HivLo) sets at FDR < 0.05; the baseline contrast (the floor) calls 285.
- fry calls 0 sets at FDR < 0.05 across the nine contrasts.

## Sheets

`set_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast by method, and fgsea collapse survivors. |
| `significant` | Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on limpa's matrix. |
