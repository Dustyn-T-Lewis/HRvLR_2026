# 07_Trait_Associations

Does a feature's link to a trait differ between the groups?

The HR/LR stage [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md) run on the strength clusters, groups Hi and Lo, read from the relabelled inputs in `01_Clustering/c_data/inputs/`. The code differs from the template only in its paths and group names.

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/DAList_corrected.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/modules.rds`, `00_Input/phenotype.csv`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/inputs/groups.rds` |
| Writes | `c_data/trait_associations.xlsx`, `b_reports/trait_associations_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/07_Trait_Associations/a_script/trait_associations.R` |

```
07_Trait_Associations/
├── a_script/trait_associations.R
├── b_reports/   trait_associations_report.pdf
└── c_data/      trait_associations.xlsx
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/07_Trait_Associations`](../../../02_HRvLR_groups/07_Trait_Associations/README.md).
