# 03_Protein_Contrasts

Which proteins differ? A copy of `02_HRvLR_groups/02_Protein_Contrasts` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/proteins.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/02_contrasts.xlsx`, `<input>/c_data/fit.rds`, `<input>/b_reports/02_contrasts_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/03_Protein_Contrasts/a_script/03_dep.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/02_Protein_Contrasts`](../../../02_HRvLR_groups/02_Protein_Contrasts/README.md).
