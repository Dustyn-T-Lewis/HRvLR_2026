# 01_Preprocess / 01_Filtering

Removes contaminant proteins by identity, drops outlier samples, then filters on detection.

| | |
|---|---|
| **Script** | `a_script/01_filter.qmd` |
| **Reads** | `00_Input/HRvLR_raw.xlsx`, `HRvLR_meta.csv`, `blood_contaminants.csv`, `HPA_annotations_full.tsv`, `RBC_proteome_reference.tsv` |
| **Writes** | `c_data/DAList_filtered.rds`, `c_data/01_filter.xlsx` |

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```

Only `DAList_filtered.rds` is read downstream. The workbook records the filter log, every
protein's call, the contamination panels and the outlier flags.

## What it does

1. Read the report and keep the highest-mean row per accession.
2. Compute the blood index (mean log2 haemoglobin per sample) and each protein's correlation with it.
3. Call every protein: curated contaminant, HPA plasma, immunoglobulin or erythrocyte protein, or keep.
   A muscle rescue (myonuclei ≥ 20 nCPM and plasma concentration < 1e9 pg/L) overturns HPA tags.
4. Build the DAList, drop samples three of four outlier methods flag.
5. Keep proteins detected in at least 5 samples of at least one Group_Time cell.

| Step | Removed | Left |
|---|---:|---:|
| raw input | | 2,400 |
| plasma | 120 | |
| immunoglobulin | 41 | |
| erythrocyte | 19 | |
| keratin | 18 | |
| complement | 12 | |
| globin | 7 | |
| leukocyte | 5 | 2,178 |
| outlier samples (S29_T1, S28_T2, S28_T3; 48 → 45) | 0 | 2,178 |
| detection | 278 | 1,900 |

Each verdict pools the curated list and the HPA rules that reach the same class.

## Notes

**Removal is by identity, not covariation.** `blood_cor` is reported beside every call and gates
nothing. Its old cut (0.45) rested on a permutation null that recomputes to 0.59 and depends on
how many samples a protein was seen in. The last section of the notebook shows this.

**The T3 biopsies carry more blood in LR.** On the log2 haemoglobin index the arm-by-T3 term is
b = −1.21 (p = 0.032). The confound does not cancel in the interaction contrasts.

**One proteoDA function is bypassed.** `filter_proteins_by_annotation()` errors on every real
DAList, so contaminant removal is a plain subset.
