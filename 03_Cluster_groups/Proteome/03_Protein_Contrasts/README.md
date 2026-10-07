# 03_DEP

Which proteins differ? A copy of `02_HRvLR_groups/02_DEP` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Cluster/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/02_Design/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/<input>/proteins.rds`, `03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/<input>/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/02_contrasts.xlsx`, `<input>/c_data/fit.rds`, `<input>/b_reports/02_contrasts_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/03_DEP/a_script/03_dep.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/02_DEP`](../../../02_HRvLR_groups/02_DEP/README.md).
