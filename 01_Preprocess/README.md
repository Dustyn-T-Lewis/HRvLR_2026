# 01 · Preprocess

Protein report to a normalised protein table, with an imputed copy for the methods that need one.

```
00_Input/HRvLR_raw.xlsx + HRvLR_meta.csv
  01_Filtering      contaminants, outlier samples, detection  -> DAList_filtered.rds
  02_Normalization  cyclic loess                              -> DAList_normalized.rds
  03_Imputation     missForest                                -> DAList_imputed.rds
```

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Normalization/a_script/02_normalize.qmd --output-dir ../b_reports
quarto render 01_Preprocess/03_Imputation/a_script/03_impute.qmd       --output-dir ../b_reports
```

Run them in order: each reads the previous one's `.rds`. All three render in under a minute.

## What comes out

`DAList_normalized.rds` is the matrix the model is fitted to: 1,900 proteins by 45 samples,
12.3% missing. It stays unimputed, because limma fits each protein on the samples where it was
seen. `DAList_imputed.rds` fills those gaps for `fry`, singscore and WGCNA, and is read by
nothing else.

## Two orderings that matter

Contaminants go **before** normalisation, because cyclic loess estimates each sample's reference
from its own intensity distribution.

Detection filtering runs **after** outlier removal, so no discarded sample counts toward a
protein's detections.
