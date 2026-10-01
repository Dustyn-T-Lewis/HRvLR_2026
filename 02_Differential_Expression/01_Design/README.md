# 01_Design

Attaches the design and the nine contrasts to the DAList and records the within-subject correlation.

| | |
|---|---|
| Reads | `01_Preprocess/03_Quantification/c_data/proteins.rds` |
| Writes | `c_data/design.rds`, `c_data/01_design.xlsx` |
| Run | `Rscript 02_Differential_Expression/01_Design/a_script/01_design.R` |

## Method

`model.matrix(~ 0 + group + blood)`: one cell mean per arm and time point, and the centred blood
index as a covariate. Subject is blocked rather than fixed, since HR and LR are different people
and a fixed subject term would absorb every between-arm contrast; limpa's `dpcDE()` estimates the
within-subject correlation. Contrasts come from `limma::makeContrasts()`. The script stops if the
design is rank deficient or a column name fails `make.names()`. `design.rds` keeps the fields
proteoDA's object used to carry (data, annotation, metadata, design and contrast matrices), so the
later steps read it unchanged.

| Contrast | Role | Definition |
|---|---|---|
| `Training_HR` | descriptive | HR_T2 - HR_T1 |
| `Training_LR` | descriptive | LR_T2 - LR_T1 |
| `Acute_HR` | descriptive | HR_T3 - HR_T2 |
| `Acute_LR` | descriptive | LR_T3 - LR_T2 |
| `Baseline_HRvLR` | floor | HR_T1 - LR_T1 |
| `Trained_HRvLR` | descriptive | HR_T2 - LR_T2 |
| `Acute_HRvLR` | descriptive | HR_T3 - LR_T3 |
| `Training_Interaction` | primary | (HR_T2 - HR_T1) - (LR_T2 - LR_T1) |
| `Acute_Interaction` | secondary | (HR_T3 - HR_T2) - (LR_T3 - LR_T2) |

## Figures

None.

## Workbook

| Sheet | Holds |
|---|---|
| `contrasts` | The nine contrasts, their definitions and roles. |
| `correlation_strata` | Within-subject correlation from limpa's fit, blocked by subject. |
