# 04_Set_Contrasts

Which gene sets shift?

The HR/LR stage [`02_HRvLR_groups/03_Set_Contrasts`](../../../02_HRvLR_groups/03_Set_Contrasts/README.md) run once per clustering input (T1, T2, T3, training, acute) on that input's clusters, C1 and C2. Each input writes to its own folder; `stack_inputs()` (`../R/inputs.R`) joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Proteome/03_Protein_Contrasts/<input>/c_data/fit.rds`, `03_Cluster_groups/Proteome/02_Design_Contrasts/<input>/c_data/design.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_limpa.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/singscore.rds`, `00_Input/phenotype.csv` |
| Writes | `<input>/c_data/set_contrasts.xlsx`, `<input>/c_data/set_tests.rds`, `<input>/b_reports/set_contrasts_report.pdf`, `c_data/set_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/04_Set_Contrasts/a_script/set_contrasts.R` |

```
04_Set_Contrasts/
├── a_script/set_contrasts.R
├── T1/ T2/ T3/ training/ acute/   one output set per clustering input
└── c_data/                        the five workbooks stacked
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/03_Set_Contrasts`](../../../02_HRvLR_groups/03_Set_Contrasts/README.md).
