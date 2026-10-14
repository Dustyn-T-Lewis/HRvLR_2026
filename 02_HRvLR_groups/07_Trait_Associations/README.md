# 07_Trait_Associations

Does a protein's, set's or module's link to a trait differ between the groups?

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/trait_associations.xlsx`, `b_reports/trait_associations_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/07_Trait_Associations/a_script/trait_associations.R` |
| Code | `00_Setup/Shared_analyses/run_group_trait_associations.R`, shared with the cluster groupings |

## Key result

- Group difference (limma, feature ~ group * trait): 128,427 tests, 5,098 at p < 0.05, 40 at FDR < 0.05.
- Within each group (Spearman): 276,612 tests, 11,758 at p < 0.05, 16 at FDR < 0.05.

## Method

- Within each group: Spearman rho with a Fisher-z 95% interval (`spearman_by_row()`).
- Between groups: limma on one value per person, `feature ~ group * trait`; the interaction is one group's slope minus the other's, with its 95% interval (`fit_subject_model()`).
- The traits that defined the groups stay out of the group-difference test. Each group needs five people with the trait measured.

## Read with care

- Opposite within-group rho with no group-difference result describes small samples, not a finding.

## Sheets

`trait_associations.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per test, level, window, trait and group. |
| `within_group` | Spearman rho with its Fisher-z 95% CI and exact p in each group, every feature and trait. |
| `group_difference` | limma, feature ~ group * trait: HR slope minus LR slope with its 95% CI, t, p, BH and Storey q. |
