# F4: continuous traits, everyone pooled

The 16 people as one group, each trait kept as a number. No panel uses the HR/LR label in a test;
panel a shows it only so the reader can see where the groups sit.

| Panel | Shows | Reads |
|---|---|---|
| a `a_trait_changes` | each person's T1 to T2 change in fibre area, muscle area and 1RM (% change), the composite hypertrophy score and total volume load; HR dark circles, LR light triangles, bar at the group median | `00_Setup/Input/phenotype.csv` through `trait_windows()` in `00_Setup/Shared_functions/person_traits.R` |
| b `b_links_vs_chance` | percent of proteins, set scores and eigengenes at p < 0.05 per window and trait (change score), white at the 5% chance gives; the count at FDR < 0.05 in bold below | `03_Continuous/01_Trait_Associations/c_data/trait_associations.rds` (`counts`, checked against `change_score`) |
| c `c_trait_pathways` | per window, one set per tree cluster (Jaccard height 0.85): the most specific, as in F3 (the smallest set whose best p is within 10-fold of the cluster's best), up to ten by rank, against every trait of the window: fill fgsea NES, size −log10 FDR. Change windows on the left, level windows on the right | `03_Continuous/02_Trait_Pathways/c_data/trait_pathways.rds` |

Window columns: T1, T2 and T3 level pair a time point's values with trait levels; "T1, gain"
pairs T1 values with the % changes; T2 - T1 and T3 - T2 pair the proteome change with the % changes.

## Read with care


- There are few modules, so one module is a large share of those tested; panel b's module column
  moves in large steps.
- Most of panel b sits at or below chance. Two screens stand out: T1 levels against leg press 1RM
  (17% of set scores, 118 at FDR < 0.05) and T2 levels against muscle area (19%, 26). Proteins
  reach FDR < 0.05 in two cells, three proteins in all.

## Supplementary figure S4: HR slope minus LR slope

Does a protein, set score or eigengene follow a trait differently in HR and LR? From
`02_HRvLR/07_Group_Trait_Slopes`: limma on feature ~ group x trait (the interaction is HR's slope
minus LR's), and Spearman's rho within each group of 8. Panel d, everyone pooled, is the protein
profile groups of `03_Continuous/03_Profile_Clusters`.

| Panel | Shows | Reads |
|---|---|---|
| a `a_slope_counts` | share of features whose slope difference sits at p < 0.05, per feature type, window and trait, white at the 5% chance gives; the count at FDR < 0.05 in bold below. The composite score defined the groups and is not tested | `group_trait_slopes.rds` (`counts`, checked against `group_difference`) |
| b `b_within_rho` | the screen holding the smallest slope-difference p: every feature's rho in HR against its rho in LR, coloured by the slope-difference call; the six strongest FDR hits named | `group_trait_slopes.rds` (`within_group`, `group_difference`) |
| c `c_top_slope` | that smallest-p feature person by person, with a least-squares line and Spearman rho per group | `group_trait_slopes.rds`, the feature's values (`blood_correction.rds`, `set_scores.rds` or `module_building.rds`), `phenotype.csv` through `window_values()` and `trait_windows()` |
| d `d_profile_clusters` | the 12 largest profile groups: their members' mean limma t in each of the 42 window-trait screens, named by the top over-represented set where it passes FDR < 0.05 | `03_Continuous/03_Profile_Clusters/c_data/profile_clusters.rds` |

Panels b and c pick their feature by rule (smallest p across every feature type, window and trait),
so they follow the data on a rerun.

### Read with care


HR and LR were split on growth, so the range of growth inside each group is narrow, and each
within-group rho rests on 8 people. Panel a is the check on panels b and c: a screen whose share
sits near 5% has no more slope differences than chance gives.

Panel c's feature (cardiac muscle tissue morphogenesis, training change against leg press gain)
does not rest on one person: dropping any one of the 8 HR people leaves HR's rho between -0.93 and
-0.96, and the HR slope (-384) stays between -281 (without HR_S03, whose leg press rose 108%) and
-515 (without HR_S27).

- Profile groups have a mean silhouette near zero (see `03_Profile_Clusters`); panel d is a way to
  read the protein profiles, not evidence of structure.

## Supplementary figure S6: elastic net on traits, everyone pooled

From `03_Continuous/04_Elastic_Net`. The lead drawn in b to d is picked by rule
(`style/elastic_net.R`): the fit with the smallest permutation p, ties broken by the higher Q²; its
pairing and trait are then shown for every feature set.

| Panel | Shows |
|---|---|
| a_fit_overview | every fit (6 pairings × traits × 3 feature sets): held-out Q², outlined at permutation p < 0.05 (solid) and < 0.10 (dashed) |
| b_selected | for the lead, each feature kept in at least 40% of 200 subsamples of 12 people, by feature set; colour: whether it rises or falls with the trait; lines at 0.4 (recurring) and 0.6 (stable) |
| c_predictions | for the lead, each person's held-out prediction against the observed trait, HR and LR marked, with Q² and p per feature set |
| d_hr_lr | for the lead's features, Spearman rho with the trait within HR and within LR, and the HR − LR slope p (`02_HRvLR/07_Group_Trait_Slopes`) |

### Inputs


| File | Used for |
|---|---|
| `03_Continuous/04_Elastic_Net/c_data/elastic_net.rds` | `summary`, `selection`, `predictions` |
| `02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds` | `within_group`, `group_difference` (d) |

### Read with care


- 9 of 81 fits are leads at permutation p < 0.10 and 2 at p < 0.05, about what chance gives (8
  and 4). The lead stands out because three feature sets agree on it, not because it is rare.
- The lead: the acute change (the bout, in trained muscle) against the muscle-area gain. Modules
  Q² 0.24 (p 0.008), proteins Q² 0.14 (p 0.019), Hallmark sets Q² 0.06 (p 0.061). Proteins kept
  in at least 60% of subsamples: OLA1, AKR1B1, JSRP1, CSNK2B, STOML2, RPL37A, TTC37; kept in 40
  to 60%: PSMA3, GSTM2, STBD1, EIF4E. All rise more after the bout in people whose muscle grew
  more. Module: yellow (91%), whose top over-represented set does not pass FDR < 0.05 and which
  does not re-form reliably when one person is left out (S5).
- Panel d: the link runs the same way within each group (rho 0.05 to 0.78 in HR, 0.39 to 0.98 in
  LR for the proteins); no HR − LR slope reaches p < 0.05 (the closest, STOML2, 0.06). It is a
  link with growth, not a difference between the groups.
- Q² of 0.14 to 0.24 is the share of held-out variation explained; it is not a predictor for new
  people.

## Run

```sh
Rscript run_all.R 04_Figures/F4                       # F4 and S4, S6
Rscript 04_Figures/F4/a_script/F4_composite.R        # F4 alone
```

## Outputs

| File | Holds |
|---|---|
| `a_script/`, `b_reports/` | F4's panel scripts and PDFs (`a_*.R`, `a_*.pdf`, ...) and its composite (`F4_composite.pdf`: the figure, 250 x 170 mm; the profile groups are in S4 d) |
| `a_script/supp/`, `b_reports/supp/` | S4: `S4_*.R`, `S4_*.pdf` and `S4_composite.pdf` (the figure, 250 x 230 mm landscape) |
| `a_script/supp/`, `b_reports/supp/` | S6: `S6_*.R`, `S6_*.pdf` and `S6_composite.pdf` (the figure, 250 × 260 mm) |
| `c_data/F4.xlsx`, `c_data/S4.xlsx`, `c_data/S6.xlsx` | one sheet per panel, the numbers drawn |
