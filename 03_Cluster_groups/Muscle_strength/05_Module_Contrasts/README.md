# 05_Module_Contrasts

Which co-expression modules shift in each of the nine contrasts?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`02_HRvLR_groups/04_Module_Contrasts`](../../../02_HRvLR_groups/04_Module_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_module_contrasts.R`.

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/c_data/fit.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/module_contrasts.xlsx`, `b_reports/module_contrasts_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/05_Module_Contrasts/a_script/module_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_module_contrasts.R` |

## Key result

- 15 modules by 9 contrasts: 10 eigengene tests at p < 0.05 and 0 at FDR < 0.05.
- fry: 12 module tests at p < 0.05 and 0 at FDR < 0.05.

## Sheets

`module_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `eigengene_tests` | limma on module eigengenes: logFC, moderated t, p, BH and Storey q per contrast. |
| `fry_tests` | fry on each module as a protein set: direction, p and FDR per contrast. |
| `membership_contrast` | Spearman rho (Fisher-z 95% CI) between member kME and protein moderated t, per module and contrast. |
| `correlation` | Within-subject correlation on the eigengenes and on limpa's matrix. |
