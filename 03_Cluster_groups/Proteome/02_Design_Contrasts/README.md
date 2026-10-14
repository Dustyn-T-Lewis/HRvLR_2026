# 02_Design_Contrasts

Which comparisons does every later model make, and how is that model built?

Groups: C1 and C2, the proteome clusters at each input.

The method is the HR/LR step's ([`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md)), run through the same function in `00_Setup/Shared_analyses/run_design_contrasts.R`. It runs once per clustering input, each into its own `<input>/` folder; `stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/proteins.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/groups/<input>.rds` |
| Writes | `<input>/c_data/design_contrasts.xlsx`, `<input>/c_data/design.rds`, `c_data/design_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/02_Design_Contrasts/a_script/design_contrasts.R` |
| Code | `00_Setup/Shared_analyses/run_design_contrasts.R` |

## Key result

- `T1`
  - 43 biopsies from 15 people: C1 9 people (25 biopsies), C2 6 people (18 biopsies).
  - Within-subject correlation from limpa's fit: 0.14.
- `T2`
  - 44 biopsies from 15 people: C1 7 people (21 biopsies), C2 8 people (23 biopsies).
  - Within-subject correlation from limpa's fit: 0.19.
- `T3`
  - 44 biopsies from 15 people: C1 5 people (15 biopsies), C2 10 people (29 biopsies).
  - Within-subject correlation from limpa's fit: 0.17.
- `training`
  - 42 biopsies from 14 people: C1 7 people (21 biopsies), C2 7 people (21 biopsies).
  - Within-subject correlation from limpa's fit: 0.21.
- `acute`
  - 44 biopsies from 15 people: C1 4 people (12 biopsies), C2 11 people (32 biopsies).
  - Within-subject correlation from limpa's fit: 0.20.

## Sheets

`design_contrasts.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `contrasts` | The nine contrasts, their definitions and roles. |
| `correlation_strata` | Within-subject correlation from limpa's fit, blocked by subject. |
