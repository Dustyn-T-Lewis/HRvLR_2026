# 02_Quantification

Fits limpa's detection probability curve, quantifies every protein in every biopsy with a standard
error, builds the blood index and removes contaminant proteins.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/01_Precursors/c_data/precursors.rds` | the precursor matrix |
| `references/contaminants.csv` | 222 contaminant proteins, frozen from `main`; see `references/README.md` |

## Run

```sh
quarto render 01_Processing/02_Quantification/a_script/quantification.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/proteins.rds` | limma `EList` from `dpcQuant()`: `E` (complete log2 values), `other$standard.error`, `genes` (precursor count, share observed), `targets` (with the blood index) |
| `c_data/quantification.xlsx` | curve coefficients, every protein, blood index, contaminants removed, genes kept for review |
| `b_reports/quantification.html` | the rendered walkthrough |

## Read with care

- The curve and the protein values are fitted on all proteins; contaminants are removed afterwards.
  limpa's documentation says nothing on contaminants; this order is the project's choice.
- 355 of 2,179 proteins rest on one precursor, and values observed in few biopsies carry wide
  standard errors. limpa keeps them all; `dpcDE()` weights each value by its precision.
- Proteins found only in the precursor data cannot be on the frozen contaminant list; the step
  prints any keratin, immunoglobulin or haemoglobin gene left (KRT19 at present) and keeps it.
