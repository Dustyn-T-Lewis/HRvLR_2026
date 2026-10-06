# 02_Design

A copy of `02_HRvLR_groups/01_Design` run on the muscle size clusters: groups Hi and Lo from `01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

Attaches the design and the nine contrasts to the DAList and records the within-subject correlation.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/proteins.rds` |
| Writes | `c_data/design.rds`, `c_data/01_design.xlsx` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/02_Design/a_script/02_design.R` |

## Method

`model.matrix(~ 0 + group + blood)`: one cell mean per arm and time point, and the centred blood
index as a covariate. Subject is blocked rather than fixed, since Hi and Lo are different people
and a fixed subject term would absorb every between-arm contrast; limpa's `dpcDE()` estimates the
within-subject correlation. Contrasts come from `limma::makeContrasts()`. The script stops if the
design is rank deficient or a column name fails `make.names()`. `design.rds` keeps the fields
proteoDA's object used to carry (data, annotation, metadata, design and contrast matrices), so the
later steps read it unchanged.

| Contrast | Role | Definition |
|---|---|---|
| `Training_Hi` | descriptive | Hi_T2 - Hi_T1 |
| `Training_Lo` | descriptive | Lo_T2 - Lo_T1 |
| `Acute_Hi` | descriptive | Hi_T3 - Hi_T2 |
| `Acute_Lo` | descriptive | Lo_T3 - Lo_T2 |
| `Baseline_HivLo` | floor | Hi_T1 - Lo_T1 |
| `Trained_HivLo` | descriptive | Hi_T2 - Lo_T2 |
| `Acute_HivLo` | descriptive | Hi_T3 - Lo_T3 |
| `Training_Interaction` | primary | (Hi_T2 - Hi_T1) - (Lo_T2 - Lo_T1) |
| `Acute_Interaction` | secondary | (Hi_T3 - Hi_T2) - (Lo_T3 - Lo_T2) |

## Figures

None.

## Workbook

| Sheet | Holds |
|---|---|
| `contrasts` | The nine contrasts, their definitions and roles. |
| `correlation_strata` | Within-subject correlation from limpa's fit, blocked by subject. |
