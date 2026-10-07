# 07_Trait_Associations

Does a phenotype link differ by group? A copy of `02_HRvLR_groups/07_Trait_Associations` run once per proteome clustering input (T1, T2, T3,
training, acute) on that input's clusters, C1 and C2, read from the relabelled inputs in
`01_Clustering/c_data/inputs/<input>/`. Each input writes to its own folder; `stack_inputs()`
joins the five workbooks into `c_data/all_inputs.xlsx` with `input` and `validity` columns
(see `../README.md`).

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/modules.rds`, `00_Input/phenotype.csv`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/groups.rds` |
| Writes | `<input>/c_data/06_associate.xlsx`, `<input>/b_reports/06_associate_figures.pdf`, `c_data/all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/07_Trait_Associations/a_script/07_associate.R` |

Method, figures and sheets are as in [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md).

`within_group` is not stacked: over five inputs it passes xlsx's row limit, so each input's own
workbook keeps it.
