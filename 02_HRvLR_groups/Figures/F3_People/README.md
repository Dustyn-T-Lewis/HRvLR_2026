# F3_People

Parked draft: the script runs and its outputs are rebuilt, but the figure is not final.

Figure 3: can one feature tell the groups apart, and does a protein's link to strength gain differ
between them?

## The idea

The figure turns from group averages to people. Panel a asks whether any single protein, set or
module separates HR from LR (`06_Classification`). Panels b to d take the protein with the
strongest group-by-trait difference (`07_Trait_Associations`) and look at it three ways.

| Panel | Draws | Rule that picks the rows |
|---|---|---|
| a | For each level (proteins, sets, modules) and each of the five classification tasks, the best single feature's AUC with its 95% interval: pROC's stratified bootstrap where the classification stage drew the feature, DeLong's otherwise (solid and dotted lines). Beside it, the feature's name and the count of features at p < 0.05 against the count chance gives (5% of those tested). | smallest Wilcoxon p per task and level, then AUC furthest from 0.5 |
| b | Every protein's group-difference t (HR slope minus LR slope, from `feature ~ group * trait`) against its t for the continuous growth-by-strength interaction from the pooled arc. The subtitle gives the Spearman correlation of the two across proteins, a description. Up to eight group-difference hits at FDR < 0.05 are named. | the lead protein's window and trait |
| c | The same proteins' within-group Spearman rho with the trait, HR on x and LR on y, with group-difference hits at FDR < 0.05 marked and named; the subtitle counts how many of them change sign between groups. | as b |
| d | The lead protein's training change against the trait, one point per person, a least-squares line per group, each group's rho in the legend and the group-difference t, p and FDR in the subtitle. | the lead: smallest group-difference p among proteins |

AUC is the chance that a random member of the first group has a higher value than a random member
of the second (0.5 no separation, 1 or 0 complete). The group difference is the interaction term
of `feature ~ group * trait`: one group's slope minus the other's. The continuous interaction in
`04_Continuous_no_group` is `feature ~ growth * strength gain` with no groups. HR and LR are a
median split of growth (the composite hypertrophy score), so the group-by-trait test is a
two-category version of the continuous one; panel b shows how far they agree.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Lead | `slice_min(p)` on `group_difference`, proteins only | Smallest p across every window and trait. |
| Panel a | `panel_auc_top()` in `R/figure_panels.R` | `auc` and `counts` sheets of the three classification workbooks; tasks under Between groups, Training and Acute strips. |
| Panel b | `panel_group_vs_continuous()` | Joins `group_difference` to the continuous `interaction` sheet by protein, window and strength test; `stats::cor(method = "spearman")`. |
| Panel c | `panel_within_rho()` | `within_group`, widened to one column per group. |
| Panel d | `panel_group_slope()`; `subject_change(..., "T1", "T2")` | Change from `DAList_corrected.rds`; trait values from `phenotype_tables()`. |
| Assembly | `patchwork`, `save_figure()` | 300 x 175 mm, `cairo_pdf`. |
| Checks | `stopifnot()` | Panel a has one row per level and task and its tested counts equal the sources; panel b has one row per protein; panel d's per-group Spearman, recomputed from the plotted points, equals the `within_group` sheet. |

## How to read the outputs

- The script completes only when the lead protein's window is `training` and its trait is
  `pct_leg_press_1rm` or `pct_leg_ext_1rm`. The continuous interaction covers only those two
  strength traits, so any other trait leaves panel b empty; panel d always draws the T2 - T1
  change, so any other window fails the rho check. Both checks stop the script rather than draw a
  wrong figure.
- Panel a's best feature per task is the extreme of hundreds of tests. A high AUC there is what a
  screen produces by chance; the count against chance beside it is the fairer summary.
- The lead protein is chosen for the smallest p of many, so its group difference is inflated by
  the choice. Panels b and c place it among all proteins for that window and trait.
- Within-group rho rests on six to eight people and is descriptive. A group-difference hit need
  not change sign between groups; when it does, the opposite signs are part of the same result,
  not separate evidence.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/06_Classification/01_Proteins/c_data/protein_classification.xlsx`, `02_HRvLR_groups/06_Classification/02_Sets/c_data/set_classification.xlsx`, `02_HRvLR_groups/06_Classification/03_Modules/c_data/module_classification.xlsx` (each `auc`, `counts`), `02_HRvLR_groups/07_Trait_Associations/c_data/trait_associations.xlsx` (`group_difference`, `within_group`), `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/protein_trait_associations.xlsx` (`interaction`), `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `b_reports/f3_people.pdf`, `c_data/f3_people.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F3_People/a_script/f3_people.R` |

```
F3_People/
├── a_script/f3_people.R
├── b_reports/f3_people.pdf
└── c_data/f3_people.xlsx
```

## Workbook: `f3_people.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `panel_sources` | a panel | `panel`, `shows`, `file`, `sheet`, `rule` |
| `panel_a` | a level and task | `feature`, `label`, `auc`, `auc_lower`, `auc_upper`, `boot_lower`, `boot_upper`, `p`, `fdr`, `n_positive`, `n_negative`, `n_nominal`, `n_tested` |
| `panel_b` | a protein | `t_group`, `p_group`, `fdr_group`, `t_continuous`, `p_continuous`, `fdr_continuous`, `hit` |
| `panel_c` | a protein | `rho_HR`, `rho_LR`, `hit` |
| `panel_d` | a person | `subject`, `value`, `trait`, `arm` |

Module preservation between the groups is in `05_Module_Preservation` and is not drawn here.
