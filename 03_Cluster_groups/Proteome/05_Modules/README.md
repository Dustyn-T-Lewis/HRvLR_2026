# 05_Modules

Which modules differ? A copy of `02_HRvLR_groups/04_Modules/01_Contrasts` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Cluster/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/<input>/modules.rds`, `03_Cluster_groups/Proteome/02_Design/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/03_DEP/<input>/c_data/fit.rds`, `03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/<input>/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/02_contrasts.xlsx`, `<input>/b_reports/02_contrasts_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/05_Modules/a_script/05_modules.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/04_Modules/01_Contrasts`](../../../02_HRvLR_groups/04_Modules/01_Contrasts/README.md).
