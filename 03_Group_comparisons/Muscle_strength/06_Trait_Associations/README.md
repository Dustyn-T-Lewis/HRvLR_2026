# 07_Trait_Associations

Does a protein's, set's or module's link to a trait differ between the groups?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`03_Group_comparisons/HRvLR/06_Trait_Associations`](../../../03_Group_comparisons/HRvLR/06_Trait_Associations/README.md)), run through the same function in `00_Setup/Shared_analyses/run_group_trait_associations.R`.

| | |
|---|---|
| Reads | `03_Group_comparisons/Muscle_strength/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/trait_associations.xlsx`, `b_reports/trait_associations_report.pdf` |
| Run | `Rscript 03_Group_comparisons/Muscle_strength/06_Trait_Associations/a_script/trait_associations.R` |
| Code | `00_Setup/Shared_analyses/run_group_trait_associations.R` |

## Key result

- Group difference (limma, feature ~ group * trait): 118,548 tests, 5,484 at p < 0.05, 15 at FDR < 0.05.
- Within each group (Spearman): 276,612 tests, 10,978 at p < 0.05, 15 at FDR < 0.05.

## Sheets

`trait_associations.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per test, level, window, trait and group. |
| `within_group` | Spearman rho with its Fisher-z 95% CI and exact p in each group, every feature and trait. |
| `group_difference` | limma, feature ~ group * trait: Hi slope minus Lo slope with its 95% CI, t, p, BH and Storey q. |
