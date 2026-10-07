# 03_Protein_Contrasts

Which proteins change, and which differ between the groups?

The HR/LR stage [`02_HRvLR_groups/02_Protein_Contrasts`](../../../02_HRvLR_groups/02_Protein_Contrasts/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/proteins.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/protein_contrasts.xlsx`, `b_reports/protein_contrasts_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/a_script/protein_contrasts.R` |

```
03_Protein_Contrasts/
├── a_script/protein_contrasts.R
├── b_reports/   protein_contrasts_report.pdf
└── c_data/      fit.rds, protein_contrasts.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/02_Protein_Contrasts`](../../../02_HRvLR_groups/02_Protein_Contrasts/README.md).
