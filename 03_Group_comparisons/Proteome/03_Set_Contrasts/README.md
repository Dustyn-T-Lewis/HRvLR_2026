# 04_Set_Contrasts

Which gene sets shift in each of the nine contrasts?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`03_Group_comparisons/HRvLR/03_Set_Contrasts`](../../../03_Group_comparisons/HRvLR/03_Set_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_set_contrasts.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `03_Group_comparisons/Proteome/02_Protein_Contrasts/<input>/c_data/fit.rds`, `03_Group_comparisons/Proteome/01_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `<input>/c_data/set_contrasts.xlsx`, `<input>/c_data/set_tests.rds`, `<input>/b_reports/set_contrasts_report.pdf`, `c_data/set_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Group_comparisons/Proteome/03_Set_Contrasts/a_script/set_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_set_contrasts.R` |

## Key result

- `T1`
  - 1,378 gene sets tested per contrast.
  - fgsea calls 22 (Acute_Interaction) to 234 (Training_C1) sets at FDR < 0.05; the baseline contrast (the floor) calls 189.
  - fry calls 207 sets at FDR < 0.05 across the nine contrasts.
- `T2`
  - 1,378 gene sets tested per contrast.
  - fgsea calls 62 (Acute_C2) to 396 (Acute_C1) sets at FDR < 0.05; the baseline contrast (the floor) calls 260.
  - fry calls 18 sets at FDR < 0.05 across the nine contrasts.
- `T3`
  - 1,378 gene sets tested per contrast.
  - fgsea calls 41 (Training_Interaction) to 471 (Acute_C1vC2) sets at FDR < 0.05; the baseline contrast (the floor) calls 220.
  - fry calls 1,216 sets at FDR < 0.05 across the nine contrasts.
- `training`
  - 1,378 gene sets tested per contrast.
  - fgsea calls 0 (Trained_C1vC2) to 209 (Acute_C2) sets at FDR < 0.05; the baseline contrast (the floor) calls 9.
  - fry calls 0 sets at FDR < 0.05 across the nine contrasts.
- `acute`
  - 1,378 gene sets tested per contrast.
  - fgsea calls 20 (Trained_C1vC2) to 343 (Acute_C1) sets at FDR < 0.05; the baseline contrast (the floor) calls 148.
  - fry calls 1,168 sets at FDR < 0.05 across the nine contrasts.

## Sheets

`set_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast by method, and fgsea collapse survivors. |
| `significant` | Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on limpa's matrix. |
