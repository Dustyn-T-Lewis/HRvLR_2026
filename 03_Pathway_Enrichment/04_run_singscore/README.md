# 04_run_singscore

One singscore per set per sample.

| | |
|---|---|
| **Reads** | `gene_sets.rds`, `DAList_imputed.rds` |
| **Writes** | `singscore.rds`, `set_scores.csv`, `04_run_singscore.xlsx`, 2 figures |

singscore ranks proteins within each sample, so a score does not depend on the cohort. Ranks need
every protein in every sample, so this reads the imputed matrix.

Subject dominates raw scores: PC1 carries 25.6% of the variance, of which subject explains 0.68.
`05_classify_and_associate_sets` therefore reads within-subject change as well as levels.

Figures: mean score against mean dispersion per set, coloured by collection; score distribution by
Group_Time.

```r
ss <- readRDS("03_Pathway_Enrichment/04_run_singscore/c_data/singscore.rds")
ss$scores            # 1,378 sets x 45 samples
ss$collection_spread # median score and dispersion per collection
ss$structure_check   # variance per component, subject share
```
