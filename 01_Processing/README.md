# 01_Processing

The shared steps every analysis folder reads: from the raw protein report to cleaned, normalised,
complete protein values, gene-set scores and co-expression modules. Building them once means a
difference between the analysis folders comes from the question, not from the preprocessing.

No step here reads a phenotype or a contrast. One step reads the HR/LR label: the detection
filter in `01_Protein_Filtering` keeps a protein seen in at least 5 samples of at least one of the
six group-by-time cells. It counts detections, never abundance.

| Step | Question | Writes for later steps |
|---|---|---|
| [`01_Protein_Filtering`](01_Protein_Filtering/README.md) | Which proteins and which runs are fit to analyse? | `DAList_filtered.rds` |
| [`02_Normalisation`](02_Normalisation/README.md) | How are the samples made comparable? | `DAList_normalized.rds` |
| [`03_Quantification`](03_Quantification/README.md) | How does every protein get a value in every sample? | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | Which groups of proteins count as one process? | `gene_sets.rds` |
| [`05_Set_Scores`](05_Set_Scores/README.md) | How strongly does each biopsy express each set? | `singscore.rds` |
| [`06_Module_Building`](06_Module_Building/README.md) | Which proteins move together within a person? | `modules.rds` |
| [`07_Module_Characterisation`](07_Module_Characterisation/README.md) | What is each module? | `module_characterisation.xlsx` |
| [`Figures/S01_QC`](Figures/S01_QC/README.md) | Quality-control supplement | `s01_qc.pdf` |

The three matrices from `03_Quantification` differ in one way each:

| File | Values | Read by |
|---|---|---|
| `proteins.rds` | limpa's complete values, each with a standard error | the protein contrasts, which weight each value by its precision |
| `DAList_limpa.rds` | the same values, no blood correction | fry set tests, which carry blood in their design |
| `DAList_corrected.rds` | each protein's blood slope removed against time point | methods that take no covariate: ranks, WGCNA, correlations, PCA |

`Figures/S01_QC` also reads `03_Group_comparisons/HRvLR/02_Protein_Contrasts`, so `run_all.R` runs it last.

```sh
Rscript run_all.R 01                                                  # every step here
Rscript 01_Processing/03_Quantification/a_script/quantification.R     # one step
```
