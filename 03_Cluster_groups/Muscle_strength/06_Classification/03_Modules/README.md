# 06_Classification/03_Modules

Can one module tell the groups apart?

The HR/LR stage [`02_HRvLR_groups/06_Classification/03_Modules`](../../../../02_HRvLR_groups/06_Classification/03_Modules/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/modules.rds` |
| Writes | `c_data/module_classification.xlsx`, `b_reports/module_classification_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classification/03_Modules/a_script/module_classification.R` |

```
03_Modules/
├── a_script/module_classification.R
├── b_reports/   module_classification_report.pdf
└── c_data/      module_classification.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/06_Classification/03_Modules`](../../../../02_HRvLR_groups/06_Classification/03_Modules/README.md).
