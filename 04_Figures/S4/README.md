# S4: HR slope minus LR slope

Does a protein, set score or eigengene follow a trait differently in HR and LR? From
`02_HRvLR/07_Group_Trait_Slopes`: limma on feature ~ group x trait (the interaction is HR's slope
minus LR's), and Spearman's rho within each group of 8. Panel d, everyone pooled, is the protein
profile groups of `03_Continuous/03_Profile_Clusters`.

| Panel | Shows | Reads |
|---|---|---|
| a `a_slope_counts` | share of features whose slope difference sits at p < 0.05, per feature type, window and trait, white at the 5% chance gives; the count at FDR < 0.05 in bold below. The composite score defined the groups and is not tested | `group_trait_slopes.rds` (`counts`, checked against `group_difference`) |
| b `b_within_rho` | the screen holding the smallest slope-difference p: every feature's rho in HR against its rho in LR, coloured by the slope-difference call; the six strongest FDR hits named | `group_trait_slopes.rds` (`within_group`, `group_difference`) |
| c `c_top_slope` | that smallest-p feature person by person, with a least-squares line and Spearman rho per group | `group_trait_slopes.rds`, the feature's values (`proteins_corrected.rds`, `singscore.rds` or `modules.rds`), `phenotype.csv` through `window_values()` and `trait_windows()` |
| d `d_profile_clusters` | the 12 largest profile groups: their members' mean limma t in each of the 42 window-trait screens, named by the top over-represented set where it passes FDR < 0.05 | `03_Continuous/03_Profile_Clusters/c_data/profile_clusters.rds` |

Panels b and c pick their feature by rule (smallest p across every feature type, window and trait),
so they follow the data on a rerun.

## Run

```sh
Rscript 04_Figures/S4/a_script/S4_composite.R
Rscript 04_Figures/S4/a_script/panels/c_top_slope.R   # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<letter>_<name>.pdf` | each panel alone |
| `b_reports/S4_composite.pdf` | the figure, 250 x 230 mm landscape |
| `c_data/S4.xlsx` | one sheet per panel with the numbers drawn, the files read with md5, packages |

## Read with care

HR and LR were split on growth, so the range of growth inside each group is narrow, and each
within-group rho rests on 8 people. Panel a is the check on panels b and c: a screen whose share
sits near 5% has no more slope differences than chance gives.

Panel c's feature (cardiac muscle tissue morphogenesis, training change against leg press gain)
does not rest on one person: dropping any one of the 8 HR people leaves HR's rho between -0.93 and
-0.96, and the HR slope (-384) stays between -281 (without HR_S03, whose leg press rose 108%) and
-515 (without HR_S27).

- Profile groups have a mean silhouette near zero (see `03_Profile_Clusters`); panel d is a way to
  read the protein profiles, not evidence of structure.
