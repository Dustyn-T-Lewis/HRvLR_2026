# 02_Groups

Decides who is compared, before any comparison runs. Four groupings, each written as one
`groups.rds` that every step in `03` reads the same way:

| Grouping | Groups | Defined by |
|---|---|---|
| `HRvLR` | HR, LR | the study: top and bottom 8 of 25 on the composite hypertrophy score |
| `Muscle_size` | Hi, Lo | PAM on % change in fibre and muscle area; reproduces HR/LR |
| `Muscle_strength` | Hi, Lo | PAM on % change in leg-press and leg-extension 1RM |
| `Proteome` | C1, C2 per input | PAM on module eigengenes at T1, T2, T3, over training and after the bout |

The gap statistic picks K = 1 for every clustering: the data hold no natural groups. K = 2 is
imposed on purpose, so every clustered comparison downstream is exploratory.

## Inputs

| File | Holds |
|---|---|
| `00_Setup/Input/phenotype.csv` | per-person size, strength, composite score and volume |
| `00_Setup/Input/metadata.csv` | run to person and label (HR/LR only) |
| `01_Processing/06_Module_Building/c_data/modules.rds` | module eigengenes (Proteome only) |

## Run

```sh
Rscript run_all.R 02_Groups
```

or one grouping, from the project root:

```sh
quarto render 02_Groups/a_script/study_groups.qmd -P grouping:HRvLR --output-dir ../HRvLR/b_reports
quarto render 02_Groups/a_script/clustering.qmd -P grouping:Proteome --output-dir ../Proteome/b_reports
```

## Outputs

| File | Holds |
|---|---|
| `HRvLR/b_reports/study_groups.html` | composite score check and labels |
| `HRvLR/c_data/groups.rds`, `study_groups.xlsx` | HR/LR membership; the recomputed composite |
| `<grouping>/b_reports/clustering.html` | gap statistic, silhouette, bootstrap Jaccard, people in cluster space |
| `<grouping>/c_data/groups.rds` | membership (`Proteome/c_data/groups/<input>.rds`, one per input) |
| `<grouping>/c_data/clustering.xlsx` | K choice and reason, diagnostics, stability, assignments; read by the F1 figures and `08_Cluster_Phenotype` |
