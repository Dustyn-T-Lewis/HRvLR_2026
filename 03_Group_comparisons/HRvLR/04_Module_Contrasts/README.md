# 04_Module_Contrasts

Which co-expression modules shift in each of the nine contrasts?

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `03_Group_comparisons/HRvLR/01_Design_Contrasts/c_data/design.rds`, `03_Group_comparisons/HRvLR/02_Protein_Contrasts/c_data/fit.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/module_contrasts.xlsx`, `b_reports/module_contrasts_report.pdf` |
| Run | `Rscript 03_Group_comparisons/HRvLR/04_Module_Contrasts/a_script/module_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_module_contrasts.R`, shared with the cluster groupings |

## Key result

- 15 modules by 9 contrasts: 9 eigengene tests at p < 0.05 and 0 at FDR < 0.05.
- fry: 8 module tests at p < 0.05 and 1 at FDR < 0.05.

## Method

- Eigengenes: `lmFit()` with the design from `design.rds`, person blocked (`duplicateCorrelation()` on the eigengenes), `eBayes(robust = TRUE)`, logFC with its 95% interval.
- fry: each module as a protein set on limpa's complete matrix.
- `membership_contrast`: within each module, Spearman between a member's kME and its protein t, with a Fisher-z interval.

## Read with care

- Few modules means BH has little to adjust over; a nominal module result is still one of several tests per contrast.

## Sheets

`module_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `eigengene_tests` | limma on module eigengenes: logFC, moderated t, p, BH and Storey q per contrast. |
| `fry_tests` | fry on each module as a protein set: direction, p and FDR per contrast. |
| `membership_contrast` | Spearman rho (Fisher-z 95% CI) between member kME and protein moderated t, per module and contrast. |
| `correlation` | Within-subject correlation on the eigengenes and on limpa's matrix. |
