# 06_Classify/01_Proteins

Does one protein separate the groups? A copy of `02_HRvLR_groups/05_Classify/01_Proteins` run on the muscle strength clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_corrected.rds`, `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_limpa.rds` |
| Writes | `c_data/01_proteins.xlsx`, `b_reports/01_proteins_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classify/01_Proteins/a_script/classify_proteins.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/05_Classify/01_Proteins`](../../../../02_HRvLR_groups/05_Classify/01_Proteins/README.md), with HR and LR read as Hi and Lo.
