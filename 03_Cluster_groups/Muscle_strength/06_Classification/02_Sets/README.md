# 06_Classification/02_Sets

Can one gene set score tell the groups apart?

The HR/LR stage [`02_HRvLR_groups/06_Classification/02_Sets`](../../../../02_HRvLR_groups/06_Classification/02_Sets/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_corrected.rds` |
| Writes | `c_data/set_classification.xlsx`, `b_reports/set_classification_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classification/02_Sets/a_script/set_classification.R` |

```
02_Sets/
├── a_script/set_classification.R
├── b_reports/   set_classification_report.pdf
└── c_data/      set_classification.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/06_Classification/02_Sets`](../../../../02_HRvLR_groups/06_Classification/02_Sets/README.md).
