# 02_Design_Contrasts

Which comparisons? A copy of `02_HRvLR_groups/01_Design_Contrasts` run on the muscle size clusters: groups Hi and Lo from
`01_Clustering`, read from the relabelled inputs in `01_Clustering/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/proteins.rds` |
| Writes | `c_data/01_design.xlsx`, `c_data/design.rds` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/02_Design_Contrasts/a_script/02_design.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md), with HR and LR read as Hi and Lo.
