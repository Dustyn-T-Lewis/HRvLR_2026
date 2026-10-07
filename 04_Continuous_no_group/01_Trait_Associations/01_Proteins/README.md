# 01_Proteins

Which proteins track each size and strength trait, with every person pooled?

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `01_Processing/03_Quantification/c_data/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/associate.rds`, `c_data/protein_trait_associations.xlsx`, `b_reports/protein_trait_associations_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Trait_Associations/01_Proteins/a_script/protein_trait_associations.R` |
| Code | `R/associate.R` |

## Key result

- Sample model: 22,800 tests, 1,175 at p < 0.05, 1 at FDR < 0.05.
- Change score: 79,800 tests, 3,551 at p < 0.05, 4 at FDR < 0.05.
- Two predictor: 22,800 tests, 1,003 at p < 0.05, 2 at FDR < 0.05.
- Interaction: 11,400 tests, 576 at p < 0.05, 20 at FDR < 0.05.

## Method

- Sample model: limma on every T1 and T2 biopsy with person blocked; `within` (a biopsy's trait against that person's mean) and `between` (the person's mean) slopes with 95% intervals.
- Change score: one value per person and window; limma slope with its interval, Spearman rho with a Fisher-z interval, and bicor.
- Two predictor (`~ growth + strength`) and interaction (`~ growth * strength`) per change window and 1RM. All in `associate_step()` (`R/associate.R`); BH within each screen.

## Read with care

- Sixteen people: read every slope with its interval.
- `within` is the stronger evidence for a link to change, because each person is their own control.

## Sheets

`protein_trait_associations.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and trait. |
| `model` | Sample model: slope per trait unit with its 95% CI, t, p, BH and Storey q per protein, trait and term. |
| `change_score` | Change score, every subject pooled: limma slope with its 95% CI, t, p, BH, Storey q; Spearman rho with its Fisher-z 95% CI and p_rho; bicor and its p; per protein, window and trait. |
| `two_predictor` | Two predictor: growth and strength slopes with 95% CIs, t, p, BH and Storey q per protein, window and test. |
| `interaction` | Interaction, feature ~ growth x strength gain: interaction term with its 95% CI, t, p, BH and Storey q per protein, window and test. |
| `plot_index` | PDF page of each overview, and page and panel position of every p < 0.05 row. |
