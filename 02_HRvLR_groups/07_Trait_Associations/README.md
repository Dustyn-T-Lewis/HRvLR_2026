# 07_Trait_Associations

Does a protein's, set's or module's link to a trait differ between the two groups?

## The idea

`04_Continuous_no_group/01_Trait_Associations` links every feature to the size and strength traits
with all people pooled and no group label. This stage asks only what the grouping adds: is the
link the same in both groups?

Two analyses answer it, and only one of them is a test of a difference.

Within group, descriptive. In each group separately, the Spearman correlation between each
feature and each trait, one value per person. Spearman's rho is the ordinary correlation of the
ranks: +1 when the person with the most of the feature also has the most of the trait, and so on
down; -1 when the order is reversed; 0 when the orders are unrelated. Ranks make it robust to one
extreme person. Each group holds six to eight people with both values, so rho is unstable, and
the groups were split on an outcome of training, which narrows the range of that outcome inside
each group and pulls rho toward zero. The p is exact (`cor.test(exact = TRUE)`), and the 95%
interval comes from Fisher's z transform (`DescTools::CorCI()`), approximate at this n.

Group difference, the test. limma fits `feature ~ group * trait` across both groups, one value
per person. That model draws one line per group: an intercept and a slope of feature on trait for
each. The interaction term `group:trait` is the difference between the two slopes. LR is the
reference level, so in this arc it is HR's slope minus LR's (`slope_difference`). Its 95% interval
and p say whether the link differs. A link significant in one group and not the
other is not, on that basis, a difference; the slopes can be close while one p crosses 0.05 and
the other does not.

The trait the groups were split on (`comp_hypertrophy` here) stays out of the group difference:
testing it would re-measure the split. A trait enters the group difference only when each group
has at least five people with it measured.

Windows pair a feature value with trait values, from `phenotype_tables()` in `R/helpers.R`:

| Window | Feature value | Traits |
|---|---|---|
| `training` | change T1 to T2 | % change in each trait, `comp_hypertrophy`, `volume_load_total_kg` |
| `baseline` | level at T1 | the same change outcomes: does where a person starts predict how they change? |
| `acute` | change T2 to T3 | the same change outcomes |
| `baseline_level` | level at T1 | trait levels at T1 |
| `post_level` | level at T2 | trait levels at T2 |
| `post_bout_level` | level at T3 | trait levels at T2 (no trait was measured at T3) |

The six traits are fibre cross-sectional area (fCSA mixed, type I, type II), muscle
cross-sectional area (mCSA), leg press 1RM and leg extension 1RM. Every analysis runs at three
levels: proteins (blood-corrected complete matrix), gene-set singscores (qualifying sets) and
module eigengenes.

FDR (BH) and Storey's q run within level, window, trait and group, and within collection for sets.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Pairings | `phenotype_tables()`, `subject_values()`, `subject_change()` in `R/helpers.R` | One value per person per window; change windows use people with both biopsies. |
| Within group | `spearman_by_row(exact = TRUE)` in `R/helpers.R`: `cor.test(method = "spearman")`, `DescTools::CorCI()` | In each group separately. Exact p for nine people or fewer; with ties `cor.test` falls back to the normal approximation. No interval below four people. |
| Group difference | `subject_fit()` in `R/helpers.R`: `model.matrix(~ group * outcome)`, `limma::lmFit()`, `eBayes(robust = TRUE)`, `topTable(confint = TRUE)` | LR is the reference level, so the `groupHR:outcome` coefficient is HR's slope minus LR's. People missing the trait are dropped. |
| Gate | `enough()` | Five people per group with the trait measured. |
| Adjustment | `p.adjust("BH")`, `storey_q()` | Within level, window, trait, group and collection. |
| Counts | `dplyr::summarise()` | Tested, p < 0.05 and BH < 0.05 per test, level, window, trait and group. |

## How to read the outputs

- `slope_difference` is in feature units per trait unit: log2 abundance, singscore or eigengene per
  percentage point for a % change trait, per unit of the trait for a level. Read it with its
  interval.
- Compare the share at p < 0.05 in `counts` with the 5% chance gives. Features are correlated, so
  one real signal can lift many rows together.
- The composite that defines HR and LR has fibre-area change among its ingredients, so the groups
  differ by construction in the range of fCSA % change. A slope difference on those traits
  compares slopes fitted over different ranges.
- Within-group rho of opposite sign in the two groups, with no group-difference result, is a
  description of small samples, not a finding.
- One value per person means no subject block is needed, and the sample size is the number of
  people, not biopsies.

## Files

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds` (proteins, metadata), `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (qualifying sets, collections), `01_Processing/06_Module_Building/c_data/modules.rds` (eigengenes), `00_Input/phenotype.csv` (traits and group) |
| Writes | `c_data/trait_associations.xlsx` (read by `Figures/F3_People`), `b_reports/trait_associations_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/07_Trait_Associations/a_script/trait_associations.R` |

```
07_Trait_Associations/
├── a_script/trait_associations.R
├── b_reports/trait_associations_report.pdf
└── c_data/trait_associations.xlsx
```

## Workbook: `trait_associations.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a test, level, window, trait (and group) | `test`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `within_group` | a feature, trait and window in one group | `level`, `collection`, `label`, `group`, `n_subjects`, `rho`, `rho_lower`, `rho_upper`, `p`, `fdr`, `qvalue` |
| `group_difference` | a feature, trait and window | `level`, `collection`, `label`, `n_subjects`, `slope_difference`, `slope_difference_lower`, `slope_difference_upper`, `t`, `p`, `fdr`, `qvalue` |

## Report: `trait_associations_report.pdf`

1. One page: the share of features at p < 0.05 for each test (within group in each group, and
   the group difference), by trait, in a grid of level against window. The dashed line marks 5%.
