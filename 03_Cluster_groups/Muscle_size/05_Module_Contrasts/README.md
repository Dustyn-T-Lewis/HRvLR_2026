# 05_Module_Contrasts

Which co-expression modules shift?

The HR/LR stage [`02_HRvLR_groups/04_Module_Contrasts`](../../../02_HRvLR_groups/04_Module_Contrasts/README.md) run on the muscle-size clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/modules.rds`, `03_Cluster_groups/Muscle_size/02_Design_Contrasts/c_data/design.rds`, `03_Cluster_groups/Muscle_size/03_Protein_Contrasts/c_data/fit.rds`, `03_Cluster_groups/Muscle_size/01_Clustering/c_data/inputs/DAList_limpa.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/module_contrasts.xlsx`, `b_reports/module_contrasts_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_size/05_Module_Contrasts/a_script/module_contrasts.R` |

```
05_Module_Contrasts/
├── a_script/module_contrasts.R
├── b_reports/   module_contrasts_report.pdf
└── c_data/      module_contrasts.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/04_Module_Contrasts`](../../../02_HRvLR_groups/04_Module_Contrasts/README.md).
