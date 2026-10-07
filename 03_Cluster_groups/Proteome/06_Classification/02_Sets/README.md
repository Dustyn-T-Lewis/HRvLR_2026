# 06_Classification/02_Sets

Does one set score separate the groups? A copy of `02_HRvLR_groups/06_Classification/02_Sets` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds` |
| Writes | `<input>/c_data/t1.xlsx`, `<input>/b_reports/t1_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/06_Classification/02_Sets/a_script/classify_sets.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/06_Classification/02_Sets`](../../../../02_HRvLR_groups/06_Classification/02_Sets/README.md).
