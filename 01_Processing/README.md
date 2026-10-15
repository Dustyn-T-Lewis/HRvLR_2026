# 01_Processing

The steps every analysis reads, built once: from DIA-NN's precursor report to complete protein
values with standard errors, following limpa's own pipeline. No step reads a phenotype or a
contrast.

| Step | Does | Writes for later steps |
|---|---|---|
| [`01_Precursors`](01_Precursors/README.md) | reads the report at q ≤ 0.01, keeps one injection per biopsy, applies limpa's identity filters | `precursors.rds`, `run_map.csv` |
| [`02_Quantification`](02_Quantification/README.md) | detection probability curve, protein values and standard errors, blood index, contaminants removed, proteins detected in fewer than 8 biopsies dropped | `proteins.rds` |
| [`03_Blood_Correction`](03_Blood_Correction/README.md) | removes each protein's blood slope from its observed values, for methods that cannot carry a covariate | `proteins_corrected.rds` |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | one protein per gene, MSigDB and GO Slim sets cut to measured genes, 15 to 500 genes, disease sets removed by name | `gene_sets.rds` |
| [`05_Set_Scores`](05_Set_Scores/README.md) | singscore for every tested set in every biopsy, on the corrected matrix | `singscore.rds` |
| [`06_Module_Building`](06_Module_Building/README.md) | WGCNA modules (signed, bicor) on person-centred corrected values; eigengenes and kME | `modules.rds` |
| [`07_Module_Characterisation`](07_Module_Characterisation/README.md) | gene-set ORA, top members and STRING interaction test per module | `module_characterisation.xlsx` |

What changed from the protein-level route on `main`: the input is the precursor report rather than
DIA-NN's protein matrix; limpa quantifies from precursors (`dpcQuant()`) rather than from proteins
(`dpcQuantByRow()`); there is no normalisation step (limpa reads DIA-NN's normalised intensities),
no precursor-level missing-value filter (limpa models missing values; its own `filterByDetection()`
runs on proteins after quantification) and no outlier removal (weak biopsies are kept
and down-weighted in the comparisons with `dpcDE(sample.weights = TRUE)`).

## Run

```sh
Rscript run_all.R 01_Processing
```
