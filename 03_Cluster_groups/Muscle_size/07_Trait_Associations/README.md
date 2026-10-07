# 07_Trait_Associations

Does a phenotype link differ by group? A copy of `02_HRvLR_groups/07_Trait_Associations` run on the muscle size clusters: groups Hi and Lo from
`01_Clustering`, read from the relabelled inputs in `01_Clustering/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/DAList_corrected.rds`, `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/modules.rds`, `00_Input/phenotype.csv`, `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/groups.rds` |
| Writes | `c_data/06_associate.xlsx`, `b_reports/06_associate_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/07_Trait_Associations/a_script/07_associate.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md), with HR and LR read as Hi and Lo.
