# 06_Classification/01_Proteins

Can one protein tell the groups apart?

The HR/LR stage [`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_corrected.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_limpa.rds` |
| Writes | `c_data/protein_classification.xlsx`, `b_reports/protein_classification_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/06_Classification/01_Proteins/a_script/protein_classification.R` |

```
01_Proteins/
├── a_script/protein_classification.R
├── b_reports/   protein_classification_report.pdf
└── c_data/      protein_classification.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md).
