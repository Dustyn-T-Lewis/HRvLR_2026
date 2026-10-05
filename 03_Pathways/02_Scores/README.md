# 01_Scores

Scores every sample on every set with singscore.

| | |
|---|---|
| Reads | `01_Gene_Sets/c_data/gene_sets.rds`, `01_Preprocess/03_Quantify/c_data/DAList_corrected.rds` |
| Writes | `c_data/singscore.rds`, `c_data/02_scores.xlsx`, `b_reports/02_scores_figures.pdf` |
| Run | `Rscript 03_Pathways/02_Scores/a_script/02_scores.R` |

## Method

`singscore::rankGenes()` then `multiScore()` on limpa's complete matrix with the blood slope
removed: one rank-based score per set and
sample. It carries no p-value and never sees a contrast.

## Figures

Mean score against dispersion per set, and the score distribution by group.

## Workbook

| Sheet | Holds |
|---|---|
| `score_summary` | Range and spread of the score matrix. |
| `collection_spread` | Median score and dispersion per collection. |
| `set_spread` | Mean score and dispersion per set across samples. |
| `structure_check` | Variance explained by the first four components, and the subject share of PC1 and PC2. |
| `set_scores` | The set by sample score matrix, as saved in singscore.rds. |
