# 04_Set_Contrasts

Which gene sets shift?

The HR/LR stage [`02_HRvLR_groups/03_Set_Contrasts`](../../../02_HRvLR_groups/03_Set_Contrasts/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/c_data/fit.rds`, `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_limpa.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/singscore.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/set_contrasts.xlsx`, `c_data/set_tests.rds`, `b_reports/set_contrasts_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/04_Set_Contrasts/a_script/set_contrasts.R` |

```
04_Set_Contrasts/
├── a_script/set_contrasts.R
├── b_reports/   set_contrasts_report.pdf
└── c_data/      set_contrasts.xlsx, set_tests.rds
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/03_Set_Contrasts`](../../../02_HRvLR_groups/03_Set_Contrasts/README.md).
