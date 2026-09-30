# 01 · Preprocess

Turns the protein report into a normalised protein table, plus an imputed copy for the methods that
need a complete matrix.

| Step | Runs | Writes |
|---|---|---|
| [`01_Filtering`](01_Filtering/README.md) | contaminant removal by identity, outlier samples, detection filter | `DAList_filtered.rds` |
| [`02_Normalization`](02_Normalization/README.md) | cyclic loess, then each protein's blood slope removed | `DAList_normalized.rds`, `DAList_unadjusted.rds` |
| [`03_Imputation`](03_Imputation/README.md) | missForest | `DAList_imputed.rds` |

Each step reads the previous step's `.rds`. Contaminants leave before normalisation, because cyclic
loess estimates each sample's reference from that sample's own intensity distribution. Every limma
fit downstream reads `DAList_normalized.rds`; only fry, singscore and WGCNA read
`DAList_imputed.rds`.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Summary/01_Candidates` gathers the nominal hits across levels.
