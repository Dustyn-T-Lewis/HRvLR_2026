# 05_Module_Contrasts

Which co-expression modules shift in each of the nine contrasts?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`03_Group_comparisons/HRvLR/04_Module_Contrasts`](../../../03_Group_comparisons/HRvLR/04_Module_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_module_contrasts.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `03_Group_comparisons/Proteome/01_Design_Contrasts/<input>/c_data/design.rds`, `03_Group_comparisons/Proteome/02_Protein_Contrasts/<input>/c_data/fit.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `<input>/c_data/module_contrasts.xlsx`, `<input>/b_reports/module_contrasts_report.pdf`, `c_data/module_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Group_comparisons/Proteome/04_Module_Contrasts/a_script/module_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_module_contrasts.R` |

## Key result

- `T1`
  - 15 modules by 9 contrasts: 26 eigengene tests at p < 0.05 and 10 at FDR < 0.05.
  - fry: 31 module tests at p < 0.05 and 12 at FDR < 0.05.
- `T2`
  - 15 modules by 9 contrasts: 26 eigengene tests at p < 0.05 and 10 at FDR < 0.05.
  - fry: 28 module tests at p < 0.05 and 12 at FDR < 0.05.
- `T3`
  - 15 modules by 9 contrasts: 24 eigengene tests at p < 0.05 and 21 at FDR < 0.05.
  - fry: 26 module tests at p < 0.05 and 19 at FDR < 0.05.
- `training`
  - 15 modules by 9 contrasts: 23 eigengene tests at p < 0.05 and 13 at FDR < 0.05.
  - fry: 28 module tests at p < 0.05 and 14 at FDR < 0.05.
- `acute`
  - 15 modules by 9 contrasts: 35 eigengene tests at p < 0.05 and 23 at FDR < 0.05.
  - fry: 36 module tests at p < 0.05 and 23 at FDR < 0.05.

## Sheets

`module_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `eigengene_tests` | limma on module eigengenes: logFC, moderated t, p, BH and Storey q per contrast. |
| `fry_tests` | fry on each module as a protein set: direction, p and FDR per contrast. |
| `membership_contrast` | Spearman rho (Fisher-z 95% CI) between member kME and protein moderated t, per module and contrast. |
| `correlation` | Within-subject correlation on the eigengenes and on limpa's matrix. |
