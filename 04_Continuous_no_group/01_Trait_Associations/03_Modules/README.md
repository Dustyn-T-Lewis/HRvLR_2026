# 03_Modules

Which module eigengenes track each size and strength trait, with every person pooled?

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/module_trait_associations.xlsx`, `b_reports/module_trait_associations_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Trait_Associations/03_Modules/a_script/module_trait_associations.R` |
| Code | `R/associate.R` |

## Key result

- Sample model: 180 tests, 6 at p < 0.05, 2 at FDR < 0.05.
- Change score: 630 tests, 21 at p < 0.05, 2 at FDR < 0.05.
- Two predictor: 180 tests, 6 at p < 0.05, 1 at FDR < 0.05.
- Interaction: 90 tests, 4 at p < 0.05, 1 at FDR < 0.05.

## Method

- Sample model: limma on every T1 and T2 biopsy with person blocked; `within` (a biopsy's trait against that person's mean) and `between` (the person's mean) slopes with 95% intervals.
- Change score: one value per person and window; limma slope with its interval, Spearman rho with a Fisher-z interval, and bicor.
- Two predictor (`~ growth + strength`) and interaction (`~ growth * strength`) per change window and 1RM. All in `associate_step()` (`R/associate.R`); BH within each screen.

## Read with care

- Sixteen people: read every slope with its interval.
- `within` is the stronger evidence for a link to change, because each person is their own control.

## Sheets

`module_trait_associations.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and trait. |
| `model` | Sample model: slope per trait unit with its 95% CI, t, p, BH and Storey q per module, trait and term. |
| `change_score` | Change score, every subject pooled: limma slope with its 95% CI, t, p, BH, Storey q; Spearman rho with its Fisher-z 95% CI and p_rho; bicor and its p; per module, window and trait. |
| `two_predictor` | Two predictor: growth and strength slopes with 95% CIs, t, p, BH and Storey q per module, window and test. |
| `interaction` | Interaction, feature ~ growth x strength gain: interaction term with its 95% CI, t, p, BH and Storey q per module, window and test. |
| `plot_index` | PDF page of each overview, and page and panel position of every p < 0.05 row. |
