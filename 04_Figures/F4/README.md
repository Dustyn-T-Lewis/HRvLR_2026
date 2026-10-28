# F4: continuous traits, everyone pooled

The 16 people as one group, each trait kept as a number. No panel uses the HR/LR label in a test;
panel a shows it only so the reader can see where the groups sit.

| Panel | Shows | Reads |
|---|---|---|
| a `a_trait_changes` | each person's T1 to T2 change in fibre area, muscle area and 1RM (% change), the composite hypertrophy score and total volume load; HR dark circles, LR light triangles, bar at the group median | `00_Setup/Input/phenotype.csv` through `trait_windows()` in `02_HRvLR/helpers.R` |
| b `b_links_vs_chance` | percent of proteins, set scores and eigengenes at p < 0.05 per window and trait (change score), white at the 5% chance gives; the count at FDR < 0.05 in bold below | `03_Continuous/01_Trait_Associations/c_data/associations.rds` (`counts`, checked against `change_score`) |
| c `c_trait_pathways` | per window, one set per tree cluster (Jaccard height 0.85): the most specific, as in F3 (the smallest set whose best p is within 10-fold of the cluster's best), up to ten by rank, against every trait of the window: fill fgsea NES, size −log10 FDR. Change windows on the left, level windows on the right | `03_Continuous/02_Trait_Pathways/c_data/trait_pathways.rds` |

Window columns: T1, T2 and T3 level pair a time point's values with trait levels; "T1, gain"
pairs T1 values with the % changes; T2 - T1 and T3 - T2 pair the proteome change with the % changes.

## Run

```sh
Rscript 04_Figures/F4/a_script/F4_composite.R
Rscript 04_Figures/F4/a_script/panels/b_links_vs_chance.R   # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<letter>_<name>.pdf` | each panel alone |
| `b_reports/F4_composite.pdf` | the figure, 250 x 170 mm; the profile groups are in S4 d |
| `c_data/F4.xlsx` | one sheet per panel with the numbers drawn, the files read with md5, packages |

## Read with care

- There are few modules, so one module is a large share of those tested; panel b's module column
  moves in large steps.
- Most of panel b sits at or below chance. Two screens stand out: T1 levels against leg press 1RM
  (17% of set scores, 118 at FDR < 0.05) and T2 levels against muscle area (19%, 26). Proteins
  reach FDR < 0.05 in two cells, three proteins in all.
