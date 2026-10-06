# 04_Pathways

Which gene sets differ? A copy of `02_HRvLR_groups/03_Pathways` run on the muscle size clusters: groups Hi and Lo from
`01_Cluster`, read from the relabelled inputs in `01_Cluster/c_data/inputs/`.

| | |
|---|---|
| Reads | `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_size/03_DEP/c_data/fit.rds`, `03_Cluster_groups/Muscle_size/02_Design/c_data/design.rds`, `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/DAList_limpa.rds`, `03_Cluster_groups/Muscle_size/01_Cluster/c_data/inputs/singscore.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/03_contrasts.xlsx`, `c_data/set_tests.rds`, `b_reports/03_contrasts_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/04_Pathways/a_script/04_pathways.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/03_Pathways`](../../../02_HRvLR_groups/03_Pathways/README.md), with HR and LR read as Hi and Lo.
