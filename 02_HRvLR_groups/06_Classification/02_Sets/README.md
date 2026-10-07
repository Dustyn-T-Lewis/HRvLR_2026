# 06_Classification/02_Sets

Can one gene set's score tell one group from the other?

## The idea

The method, the five tasks, the AUC and its two intervals are explained in
[`06_Classification`](../README.md). This substage runs them on every gene set that qualifies
for testing (15 to 500 measured proteins).

The value is the set's singscore in each biopsy, from `01_Trunk/05_Set_Scores`. singscore ranks
all proteins within one biopsy and scores a set by where its members sit in that ranking: high
when they sit near the top, low near the bottom. It needs no other biopsy and no group label, so
every biopsy gets one score per set, computed on the blood-corrected complete matrix.

The five collections (Hallmark, Reactome, KEGG_Legacy, GOBP, GO_Slim) differ in size and overlap
heavily with each other, so each collection is its own screen: BH and Storey's q run within task
and collection.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Values | `readRDS()` | `scores` from `singscore.rds`, rows limited to the qualifying sets in `gene_sets.rds`; biopsy metadata from `DAList_corrected.rds`. The script stops if the columns are out of order. |
| Labels | `enrichVolcano::ev_clean_label()` | Readable set names; report headers read "collection: set name", wrapped at 50 characters. |
| Row order | `hclust(as.dist(1 - cor(t(...))), "average")` | On singscore, within each collection, for the overviews only. |
| Scoring | `classify_step(collection_of = ...)` in `R/classify.R` | As described in the parent README, with collection added to the BH scope, the counts and the plot index. |

## How to read the outputs

- `counts` has one row per task and collection. Compare each collection's `n_nominal` with 5% of
  its own `n_tested`.
- Overlapping sets share members, so several sets reaching p < 0.05 together can be one signal
  counted several times.
- A set's AUC reflects its members' joint position in each biopsy, not their average change.
  `03_Set_Contrasts` asks the second question.

## Files

| | |
|---|---|
| Reads | `01_Trunk/05_Set_Scores/c_data/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` (set catalogue), `01_Trunk/03_Quantification/c_data/DAList_corrected.rds` (biopsy metadata) |
| Writes | `c_data/set_classification.xlsx` (read by `Figures/F3_People`), `b_reports/set_classification_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/02_Sets/a_script/set_classification.R` |

```
02_Sets/
├── a_script/set_classification.R
├── b_reports/set_classification_report.pdf
└── c_data/set_classification.xlsx
```

## Workbook: `set_classification.xlsx`

Sheets as in the parent README, each with a `collection` column.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a task and collection | `n_positive`, `n_negative`, `min_p`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `auc` | a set in a task | `collection`, `label`, `auc`, `auc_lower`, `auc_upper`, `p`, `fdr`, `qvalue`, `boot_lower`, `boot_upper` |
| `plot_index` | an overview page or a drawn set | `task`, `collection`, `feature`, `pdf_page`, `position` |

## Report: `set_classification_report.pdf`

1. One overview page per collection: AUC for every set in the five tasks beside a p/FDR status
   map.
2. ROC curves for the sets at p < 0.05 in each task, at most 24 per task across collections.
