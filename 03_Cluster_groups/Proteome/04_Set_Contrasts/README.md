# 04_Set_Contrasts

Which gene sets differ? A copy of `02_HRvLR_groups/03_Set_Contrasts` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Proteome/03_Protein_Contrasts/<input>/c_data/fit.rds`, `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_limpa.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/singscore.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/03_contrasts.xlsx`, `<input>/c_data/set_tests.rds`, `<input>/b_reports/03_contrasts_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/04_Set_Contrasts/a_script/04_pathways.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/03_Set_Contrasts`](../../../02_HRvLR_groups/03_Set_Contrasts/README.md).
