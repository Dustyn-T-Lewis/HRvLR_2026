# 06_Classify/02_Sets

Does one set score separate the groups? A copy of `02_HRvLR_groups/05_Classify/02_Sets` run on the muscle strength clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_strength/01_Cluster/c_data/inputs/DAList_corrected.rds` |
| Writes | `c_data/02_sets.xlsx`, `b_reports/02_sets_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classify/02_Sets/a_script/classify_sets.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/05_Classify/02_Sets`](../../../../02_HRvLR_groups/05_Classify/02_Sets/README.md), with HR and LR read as Hi and Lo.
