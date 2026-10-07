# 05_Module_Contrasts

Which co-expression modules shift?

The HR/LR stage [`02_HRvLR_groups/04_Module_Contrasts`](../../../02_HRvLR_groups/04_Module_Contrasts/README.md) run once per clustering input (T1, T2, T3, training, acute) on that input's clusters, C1 and C2. Each input writes to its own folder; `stack_inputs()` (`../R/inputs.R`) joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/modules.rds`, `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/03_Protein_Contrasts/<input>/c_data/fit.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/module_contrasts.xlsx`, `<input>/b_reports/module_contrasts_report.pdf`, `c_data/module_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/05_Module_Contrasts/a_script/module_contrasts.R` |

```
05_Module_Contrasts/
├── a_script/module_contrasts.R
├── T1/ T2/ T3/ training/ acute/   one output set per clustering input
└── c_data/                        the five workbooks stacked
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/04_Module_Contrasts`](../../../02_HRvLR_groups/04_Module_Contrasts/README.md).
