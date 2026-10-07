# 02_Design_Contrasts

Which comparisons does every later model make, and how is that model built?

Groups: Hi and Lo, the strength clusters.

The method is the HR/LR step's ([`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md)), run through the same function in `R/design_contrasts.R`.

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/proteins.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/groups.rds` |
| Writes | `c_data/design_contrasts.xlsx`, `c_data/design.rds` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/02_Design_Contrasts/a_script/design_contrasts.R` |
| Code | `R/design_contrasts.R` |

## Key result

- 42 biopsies from 15 people: Hi 8 people (22 biopsies), Lo 7 people (20 biopsies).
- Within-subject correlation from limpa's fit: 0.17.

## Sheets

`design_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `contrasts` | The nine contrasts, their definitions and roles. |
| `correlation_strata` | Within-subject correlation from limpa's fit, blocked by subject. |
