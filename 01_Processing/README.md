# 01_Processing

The steps every analysis reads, built once: from DIA-NN's precursor report to complete protein
values with standard errors, following limpa's own pipeline. No step reads a phenotype or a
contrast.

| Step | Does | Writes for later steps |
|---|---|---|
| [`01_Precursors`](01_Precursors/README.md) | reads the report at q ≤ 0.01, keeps one injection per biopsy, applies limpa's identity filters | `precursors.rds`, `run_map.csv` |
| [`02_Quantification`](02_Quantification/README.md) | detection probability curve, protein values and standard errors, blood index, contaminants removed | `proteins.rds` |
| `04_Gene_Sets` to `07_Module_Characterisation` | carried over from `main`; still read the old processing files and are not run until re-pointed | |

What changed from the protein-level route on `main`: the input is the precursor report rather than
DIA-NN's protein matrix; limpa quantifies from precursors (`dpcQuant()`) rather than from proteins
(`dpcQuantByRow()`); there is no normalisation step (limpa reads DIA-NN's normalised intensities),
no detection filter (limpa models missing values) and no outlier removal (weak biopsies are kept
and down-weighted in the comparisons with `dpcDE(sample.weights = TRUE)`).

## Run

```sh
Rscript run_all.R 01_Processing
```
