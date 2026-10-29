# S6: elastic net on traits, everyone pooled

From `03_Continuous/04_Elastic_Net`. The lead drawn in b to d is picked by rule
(`style/elastic_net.R`): the fit with the smallest permutation p, ties broken by the higher Q²; its
pairing and trait are then shown for every feature set.

| Panel | Shows |
|---|---|
| a_fit_overview | every fit (6 pairings × traits × 3 feature sets): held-out Q², outlined at permutation p < 0.05 (solid) and < 0.10 (dashed) |
| b_selected | for the lead, each feature kept in at least 40% of 200 subsamples of 12 people, by feature set; colour: whether it rises or falls with the trait; lines at 0.4 (recurring) and 0.6 (stable) |
| c_predictions | for the lead, each person's held-out prediction against the observed trait, HR and LR marked, with Q² and p per feature set |
| d_hr_lr | for the lead's features, Spearman rho with the trait within HR and within LR, and the HR − LR slope p (`02_HRvLR/07_Group_Trait_Slopes`) |

## Inputs

| File | Used for |
|---|---|
| `03_Continuous/04_Elastic_Net/c_data/elastic_net.rds` | `summary`, `selection`, `predictions` |
| `02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds` | `within_group`, `group_difference` (d) |

## Run

```sh
Rscript run_all.R 04_Figures/S6
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel alone |
| `b_reports/S6_composite.pdf` | the figure, 250 × 260 mm |
| `c_data/S6.xlsx` | one sheet per panel with the numbers drawn |

## Read with care

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
