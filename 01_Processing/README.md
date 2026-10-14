# 01_Processing

The steps every analysis reads, built once: from the raw protein report to cleaned, normalised,
complete protein values, gene-set scores and co-expression modules. No step reads a phenotype or a
contrast. One reads the HR/LR label: the detection filter counts detections per group-by-time
cell, never abundance.

| Step | Does | Writes for later steps |
|---|---|---|
| [`01_Protein_Filtering`](01_Protein_Filtering/README.md) | removes contaminants, outlier runs and rarely detected proteins | `DAList_filtered.rds` |
| [`02_Normalisation`](02_Normalisation/README.md) | cyclic loess across runs | `DAList_normalized.rds` |
| [`03_Quantification`](03_Quantification/README.md) | limpa values and standard errors for every protein in every run | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | decides set membership once | `gene_sets.rds` |
| [`05_Set_Scores`](05_Set_Scores/README.md) | singscore per set per run | `singscore.rds` |
| [`06_Module_Building`](06_Module_Building/README.md) | WGCNA modules on within-person variation | `modules.rds` |
| [`07_Module_Characterisation`](07_Module_Characterisation/README.md) | gene sets, top members and STRING edges per module | `module_characterisation.xlsx` |

The three matrices from `03_Quantification` differ in one way each:

| File | Values | Read by |
|---|---|---|
| `proteins.rds` | limpa's values, each with a standard error | protein contrasts, which weight each value by its precision |
| `DAList_limpa.rds` | the same values, no blood correction | fry set tests, which carry blood in their design |
| `DAList_corrected.rds` | each protein's blood slope removed against time point | methods that take no covariate: ranks, WGCNA, correlations, PCA |

## Run

```sh
Rscript run_all.R 01_Processing
```

Each step's README gives its own render command. The QC supplement drawn from these outputs is
`05_Figures/Processing/S01_QC`.
