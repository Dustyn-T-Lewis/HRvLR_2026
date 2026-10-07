# 05_Set_Scores

How strongly is each gene set expressed in each biopsy?

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/set_scores.xlsx`, `c_data/singscore.rds`, `b_reports/set_scores_report.pdf` |
| Run | `Rscript 01_Processing/05_Set_Scores/a_script/set_scores.R` |

## Key result

- 1,378 sets scored in 45 samples.
- Subject identity explains 58% of PC1 and 45% of PC2 of the scores.

## Method

- `singscore::rankGenes()` ranks the proteins within each sample; `singscore::multiScore()` gives one score per qualifying set and sample. Every set is scored as an up-set.
- `structure_check`: PCA of the samples on their scores and the share of each component that subject identity explains.

## Read with care

- Scores follow the person more than the time point (`subject_share`), which is why the analyses block on subject and use within-person change.

## Sheets

`set_scores.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `score_summary` | Range and spread of the score matrix. |
| `collection_spread` | Median score and dispersion per collection. |
| `set_spread` | Mean score and dispersion per set across samples. |
| `structure_check` | Variance explained by the first four components, and the subject share of PC1 and PC2. |
| `set_scores` | The set by sample score matrix, as saved in singscore.rds. |
