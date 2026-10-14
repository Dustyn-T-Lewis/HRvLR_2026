# 05_Set_Scores

Scores every run on every gene set with singscore: the average within-run rank of the set's
members. Descriptive only; no test.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the 1,378 tested sets and the protein-to-gene map |
| `01_Processing/03_Quantification/c_data/DAList_corrected.rds` | limpa's complete matrix with blood slopes removed |

## Run

```sh
quarto render 01_Processing/05_Set_Scores/a_script/set_scores.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/set_scores.html` | the walkthrough: scores, dispersion, person share of the variation |
| `c_data/singscore.rds` | the set-by-run score matrix; read by set contrasts, classification and trait steps |
| `c_data/set_scores.xlsx` | score summaries, per-set spread, PCA check, the full matrix |
