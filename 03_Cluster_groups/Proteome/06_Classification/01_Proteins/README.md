# 06_Classification/01_Proteins

Can one protein tell the groups apart?

The HR/LR stage [`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md) run once per clustering input (T1, T2, T3, training, acute) on that input's clusters, C1 and C2. Each input writes to its own folder; `stack_inputs()` (`../R/inputs.R`) joins the five workbooks into one, adding `input` and `validity` (circular, caution or clean) to every row.

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_corrected.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/DAList_limpa.rds` |
| Writes | `<input>/c_data/protein_classification.xlsx`, `<input>/b_reports/protein_classification_report.pdf`, `c_data/protein_classification_all_inputs.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/06_Classification/01_Proteins/a_script/protein_classification.R` |

```
01_Proteins/
├── a_script/protein_classification.R
├── T1/ T2/ T3/ training/ acute/   one output set per clustering input
└── c_data/                        the five workbooks stacked
```

The method, the sheets, the columns and the report pages are described in [`02_HRvLR_groups/06_Classification/01_Proteins`](../../../../02_HRvLR_groups/06_Classification/01_Proteins/README.md).
