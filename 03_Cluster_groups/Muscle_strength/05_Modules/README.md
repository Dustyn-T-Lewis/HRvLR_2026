# 05_Modules

Which modules differ? A copy of `02_HRvLR_groups/04_Modules/01_Contrasts` run on the muscle strength clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/modules.rds`, `03_Cluster_groups/Muscle_strength/02_Design/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/03_DEP/c_data/fit.rds`, `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/05_Modules/a_script/05_modules.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/04_Modules/01_Contrasts`](../../../02_HRvLR_groups/04_Modules/01_Contrasts/README.md), with HR and LR read as Hi and Lo.
