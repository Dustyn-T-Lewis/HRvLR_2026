# 07_Trait_Associations

Does a feature's link to a trait differ between the groups?

The HR/LR stage [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md) run once per clustering input (T1, T2, T3, training, acute) on that input's clusters, C1 and C2. Each input writes to its own folder; `stack_inputs()` (`../R/inputs.R`) joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/modules.rds`, `00_Input/phenotype.csv`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/groups.rds` |
| Writes | `<input>/c_data/trait_associations.xlsx`, `<input>/b_reports/trait_associations_report.pdf`, `c_data/trait_associations_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/07_Trait_Associations/a_script/trait_associations.R` |

```
07_Trait_Associations/
├── a_script/trait_associations.R
├── T1/ T2/ T3/ training/ acute/   one output set per clustering input
└── c_data/                        the five workbooks stacked
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md).
