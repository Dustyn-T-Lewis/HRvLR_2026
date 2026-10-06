# 06_Classify/03_Modules

Does one module separate the groups? A copy of `02_HRvLR_groups/05_Classify/03_Modules` run on the muscle strength clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/modules.rds` |
| Writes | `c_data/03_modules.xlsx`, `b_reports/03_modules_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classify/03_Modules/a_script/classify_modules.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/05_Classify/03_Modules`](../../../../02_HRvLR_groups/05_Classify/03_Modules/README.md), with HR and LR read as Hi and Lo.
