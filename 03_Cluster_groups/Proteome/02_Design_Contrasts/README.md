# 02_Design_Contrasts

Which comparisons? A copy of `02_HRvLR_groups/01_Design_Contrasts` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/proteins.rds` |
| Writes | `<input>/c_data/01_design.xlsx`, `<input>/c_data/design.rds`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/02_Design_Contrasts/a_script/02_design.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md).
