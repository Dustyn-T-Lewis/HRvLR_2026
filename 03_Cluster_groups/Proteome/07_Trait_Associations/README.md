# 07_Trait_Associations

Does a protein's, set's or module's link to a trait differ between the groups?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md)), run through the same function in `00_Setup/Shared_analyses/run_group_trait_associations.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `<input>/c_data/trait_associations.xlsx`, `<input>/b_reports/trait_associations_report.pdf`, `c_data/trait_associations_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/07_Trait_Associations/a_script/trait_associations.R` |
| Code | `00_Setup/Shared_analyses/run_group_trait_associations.R` |

## Key result

- `T1`
  - Group difference (limma, feature ~ group * trait): 138,306 tests, 6,358 at p < 0.05, 89 at FDR < 0.05.
  - Within each group (Spearman): 276,612 tests, 8,681 at p < 0.05, 6 at FDR < 0.05.
- `T2`
  - Group difference (limma, feature ~ group * trait): 138,306 tests, 7,494 at p < 0.05, 6 at FDR < 0.05.
  - Within each group (Spearman): 276,612 tests, 12,502 at p < 0.05, 23 at FDR < 0.05.
- `T3`
  - Group difference (limma, feature ~ group * trait): 121,841 tests, 11,485 at p < 0.05, 821 at FDR < 0.05.
  - Within each group (Spearman): 276,612 tests, 9,038 at p < 0.05, 61 at FDR < 0.05.
- `training`
  - Group difference (limma, feature ~ group * trait): 138,306 tests, 5,425 at p < 0.05, 14 at FDR < 0.05.
  - Within each group (Spearman): 276,612 tests, 11,363 at p < 0.05, 27 at FDR < 0.05.
- `acute`
  - Group difference (limma, feature ~ group * trait): 0 tests, 0 at p < 0.05, 0 at FDR < 0.05.
  - Within each group (Spearman): 276,612 tests, 5,260 at p < 0.05, 21 at FDR < 0.05.

## Sheets

`trait_associations.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per test, level, window, trait and group. |
| `within_group` | Spearman rho with its Fisher-z 95% CI and exact p in each group, every feature and trait. |
| `group_difference` | limma, feature ~ group * trait: C1 slope minus C2 slope with its 95% CI, t, p, BH and Storey q. |
