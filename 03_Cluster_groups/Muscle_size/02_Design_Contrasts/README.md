# 02_Design_Contrasts

Which comparisons?

The HR/LR stage [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md) run on the muscle-size clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/proteins.rds` |
| Writes | `c_data/design.rds`, `c_data/design_contrasts.xlsx` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/02_Design_Contrasts/a_script/design_contrasts.R` |

```
02_Design_Contrasts/
├── a_script/design_contrasts.R
├── b_reports/   
└── c_data/      design.rds, design_contrasts.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/01_Design_Contrasts`](../../../02_HRvLR_groups/01_Design_Contrasts/README.md).
