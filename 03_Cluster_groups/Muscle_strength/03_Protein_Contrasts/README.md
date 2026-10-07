# 03_Protein_Contrasts

Which proteins differ? A copy of `02_HRvLR_groups/02_Protein_Contrasts` run on the muscle strength clusters: groups Hi and Lo from
`01_Clustering`, read from the relabelled inputs in `01_Clustering/c_data/inputs/`.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/proteins.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/02_contrasts.xlsx`, `c_data/fit.rds`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/a_script/03_dep.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/02_Protein_Contrasts`](../../../02_HRvLR_groups/02_Protein_Contrasts/README.md), with HR and LR read as Hi and Lo.
