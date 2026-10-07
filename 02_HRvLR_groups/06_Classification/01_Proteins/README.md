# 06_Classification/01_Proteins

Can one protein's abundance tell one group from the other?

## The idea

The method, the five tasks, the AUC and its two intervals are explained in
[`06_Classification`](../README.md). This substage runs them on every protein.

The values are limpa's complete log2 abundance with the blood slope removed
(`DAList_corrected.rds`). limpa fills every protein in every biopsy, estimating the values the
mass spectrometer did not record from how detection depends on abundance. A rank test cannot tell
an estimated value from a measured one, so each row carries `observed_share`, the share of the
protein's values that were measured. A protein that separates the groups mostly through estimated
values deserves less trust than one measured throughout.

BH and Storey's q run within task, over all proteins.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Values | `readRDS()` | `dal$data` from `DAList_corrected.rds`; `observed_share` is `rowMeans(dal$observed)`. |
| Labels | gene symbol per UniProt accession | The script stops if a symbol is missing or repeated. |
| Row order | `hclust(as.dist(1 - cor(t(...))), "average")` | On limpa's complete matrix (`DAList_limpa.rds`), for the overview only. |
| Scoring | `classify_step()` in `R/classify.R` | As described in the parent README. |

## How to read the outputs

- Check `observed_share` before reading a protein's AUC.
- With about 1,900 proteins in each task, the count expected at p < 0.05 by chance alone is about
  95. Compare `n_nominal` in `counts` with that.
- A protein's AUC here and its logFC in `02_Protein_Contrasts` answer different questions: the
  logFC compares group averages under the model, the AUC asks how often individuals are ordered
  by group.

## Files

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds` (values, metadata, `observed`), `01_Processing/03_Quantification/c_data/DAList_limpa.rds` (overview row order) |
| Writes | `c_data/protein_classification.xlsx` (read by `Figures/F3_People`), `b_reports/protein_classification_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/06_Classification/01_Proteins/a_script/protein_classification.R` |

```
01_Proteins/
├── a_script/protein_classification.R
├── b_reports/protein_classification_report.pdf
└── c_data/protein_classification.xlsx
```

## Workbook: `protein_classification.xlsx`

Sheets as in the parent README. `auc` adds `label` (gene symbol) and `observed_share`.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a task | `n_positive`, `n_negative`, `min_p`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `auc` | a protein in a task | `label`, `observed_share`, `auc`, `auc_lower`, `auc_upper`, `p`, `fdr`, `qvalue`, `boot_lower`, `boot_upper` |
| `plot_index` | an overview page or a drawn protein | `task`, `feature`, `pdf_page`, `position` |

## Report: `protein_classification_report.pdf`

1. One overview page: AUC for every protein in the five tasks beside a p/FDR status map, rows in
   clustering order.
2. ROC curves for the proteins at p < 0.05 in each task, at most 24 per task.
