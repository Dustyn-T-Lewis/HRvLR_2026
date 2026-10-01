# 01 · Preprocess

Turns the protein report into a normalised protein table and models its missing values with limpa,
so every protein has a value in every sample with a standard error.

| Step | Runs | Writes |
|---|---|---|
| [`01_Filtering`](01_Filtering/README.md) | contaminant removal by identity, outlier samples, detection filter | `DAList_filtered.rds` |
| [`02_Normalization`](02_Normalization/README.md) | cyclic loess | `DAList_normalized.rds` |
| [`03_Quantification`](03_Quantification/README.md) | limpa detection curve and row-wise quantification; blood-corrected copy | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |

Each step reads the previous step's `.rds`. Contaminants leave before normalisation, because cyclic
loess estimates each sample's reference from that sample's own intensity distribution. Blood is
not subtracted here: every limma model carries it as a covariate, and the methods that cannot take
one read `DAList_corrected.rds`.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Summary/02_Candidates` gathers the nominal hits across levels.
