# 02_Design_Contrasts

Which comparisons?

The HR/LR stage [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md) run once per clustering input (T1, T2, T3, training, acute) on that input's clusters, C1 and C2. Each input writes to its own folder; `stack_inputs()` (`../R/inputs.R`) joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/proteins.rds` |
| Writes | `<input>/c_data/design.rds`, `<input>/c_data/design_contrasts.xlsx`, `c_data/design_contrasts_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/02_Design_Contrasts/a_script/design_contrasts.R` |

```
02_Design_Contrasts/
├── a_script/design_contrasts.R
├── T1/ T2/ T3/ training/ acute/   one output set per clustering input
└── c_data/                        the five workbooks stacked
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md).
