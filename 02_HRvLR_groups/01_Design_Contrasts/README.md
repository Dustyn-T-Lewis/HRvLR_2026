# 01_Design_Contrasts

Which comparisons does every later model make, and how is that model built?

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/proteins.rds` |
| Writes | `c_data/design_contrasts.xlsx`, `c_data/design.rds` |
| Run | `Rscript 02_HRvLR_groups/01_Design_Contrasts/a_script/design_contrasts.R` |
| Code | `R/design_contrasts.R`, shared with the cluster groupings |

## Key result

- 45 biopsies from 16 people: HR 8 people (21 biopsies), LR 8 people (24 biopsies).
- Within-subject correlation from limpa's fit: 0.17.

## Method

- `model.matrix(~ 0 + group + blood)`: six cell means (group by time point) and the blood index.
- Nine contrasts (`limma::makeContrasts()`): each group's training (T2 - T1) and acute (T3 - T2) change, the groups compared at T1, T2 and T3, and the two interactions (one group's change minus the other's).
- Roles: `Training_Interaction` primary, `Acute_Interaction` secondary, the baseline contrast the floor (the groups come from an outcome, so a T1 difference shows how far the split separates people before training).
- `design.rds` also names the two groups (`arms`), who is in each (`group_of`) and the traits that defined them; every later step of the arc reads it.

## Sheets

`design_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `contrasts` | The nine contrasts, their definitions and roles. |
| `correlation_strata` | Within-subject correlation from limpa's fit, blocked by subject. |
