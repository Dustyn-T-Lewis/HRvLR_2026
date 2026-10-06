# 03_DEP

Which proteins differ? A copy of `02_HRvLR_groups/02_DEP` run on the muscle strength clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/02_Design/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/proteins.rds`, `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/02_contrasts.xlsx`, `c_data/fit.rds`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/03_DEP/a_script/03_dep.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/02_DEP`](../../../02_HRvLR_groups/02_DEP/README.md), with HR and LR read as Hi and Lo.
