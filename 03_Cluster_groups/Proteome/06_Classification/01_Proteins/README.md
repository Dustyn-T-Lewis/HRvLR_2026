# 06_Classification/01_Proteins

Does one protein separate the groups? A copy of `02_HRvLR_groups/06_Classification/01_Proteins` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_limpa.rds` |
| Writes | `<input>/c_data/t1.xlsx`, `<input>/b_reports/t1_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/06_Classification/01_Proteins/a_script/classify_proteins.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md).
