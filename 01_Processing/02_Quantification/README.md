# 02_Quantification

Fits limpa's detection probability curve, quantifies every protein in every biopsy with a standard
error, builds the blood index, removes contaminant proteins and keeps proteins detected in at
least 8 biopsies (`filterByDetection()`).

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
| `c_data/quantification.xlsx` | curve coefficients, every protein kept, blood index, contaminants removed, genes kept for review, proteins dropped for detection |
| `b_reports/quantification.html` | the rendered walkthrough |

## Read with care

- The curve and the protein values are fitted on all proteins; contaminants are removed afterwards.
  limpa's documentation says nothing on contaminants; this order is the project's choice.
- 2,059 proteins remain; 276 rest on one precursor. Values a biopsy did not detect come from the
  detection curve and carry wide standard errors; `dpcDE()` weights each value by its precision.
- The thin biopsies (S28-T3, S28-T2, S10-T2) detect few proteins, so most of their values come from
  the curve and sit low (median 16.3 and 16.8 against 17.8). The protein density plot shows them;
  the comparisons' sample weights should down-weight them, which is checked there.
- Proteins found only in the precursor data cannot be on the frozen contaminant list; the step
  prints any keratin, immunoglobulin or haemoglobin gene left (KRT19 at present) and keeps it.
