# 01 · Preprocess

Turns the protein report into a normalised protein table, plus an imputed copy for the methods
that need a complete matrix.

| Step | Runs | Writes |
|---|---|---|
| [`01_Filtering`](01_Filtering/README.md) | contaminant removal, outlier samples, detection filter | `DAList_filtered.rds` |
| [`02_Normalization`](02_Normalization/README.md) | cyclic loess | `DAList_normalized.rds` |
| [`03_Imputation`](03_Imputation/README.md) | missForest | `DAList_imputed.rds` |

Each step reads the previous step's `.rds`.

## Contaminants leave by identity before normalisation

Contaminants go before normalisation, because cyclic loess estimates each sample's reference from
that sample's own intensity distribution. Removal is by identity: a curated list plus HPA plasma,
immunoglobulin and erythrocyte tags, with a muscle rescue (myonuclei ≥ 20 nCPM and plasma
concentration < 1e9 pg/L) that overturns an HPA tag. `blood_cor` is reported beside every call and
gates nothing; `00_Input/PROVENANCE.md` records why its old cut was retired.

A sample goes when three of four outlier methods flag it. The detection filter runs last, so no
discarded sample counts toward a protein's detections: a protein needs 5 detections in at least
one group cell.

| Step (as `filter_log` orders them) | Removed | Left |
|---|---:|---:|
| raw input | | 2,400 |
| duplicate accession | 0 | 2,400 |
| complement | 12 | 2,388 |
| erythrocyte | 19 | 2,369 |
| globin | 7 | 2,362 |
| immunoglobulin | 41 | 2,321 |
| keratin | 18 | 2,303 |
| leukocyte | 5 | 2,298 |
| plasma | 120 | 2,178 |
| outlier samples (S29_T1, S28_T2, S28_T3; 48 → 45) | 0 | 2,178 |
| detection | 278 | 1,900 |

The T3 biopsies carry more blood in LR. In a mixed model on the log2 haemoglobin index, LR's rise
at T3 exceeds HR's by 1.21 (t = 2.22, subject-label permutation p = 0.023). On the blood panel's
share of signal the same term gives p = 0.356. This confound does not cancel in the interaction
contrasts.

## The fitted matrix stays unimputed

`normalize_data("cycloess")` is `limma::normalizeCyclicLoess(method = "fast")`. With
`adaptive.span = TRUE` the span comes from `chooseLowessSpan(nrow)`, not the 0.7 in the formals, so
it moves whenever filtering changes the protein count.

`DAList_normalized.rds` holds 1,900 proteins by 45 samples, 12.3% missing. `DAList_imputed.rds`
fills the gaps with missForest (ranger backend, 10 iterations, 100 trees, seed 42, rows sorted
before the fit; out-of-bag NRMSE 0.151) for `fry`, singscore and WGCNA.
