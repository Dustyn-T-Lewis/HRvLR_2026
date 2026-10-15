# 05_Set_Scores

Scores every biopsy on every tested gene set with singscore: the average within-biopsy rank of the
set's members. Descriptive only; no test.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the 1,422 tested sets and the protein-to-gene map |
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | limpa's complete matrix with blood slopes removed from observed values |

## Run

```sh
quarto render 01_Processing/05_Set_Scores/a_script/set_scores.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/singscore.rds` | `scores`: the set-by-biopsy score matrix (1,422 × 48) |
| `c_data/set_scores.xlsx` | score summaries, per-set spread, PCA check, the full matrix |
| `b_reports/set_scores.html` | the rendered walkthrough |

## Read with care

- Person identity explains 70% of PC1 and 41% of PC2 of the score matrix. Comparisons of scores
  must account for the person.
- Values limpa estimated from the detection curve enter the ranking as they are; a set whose
  members are rarely observed is scored partly on estimates.
- Scores of overlapping sets (GO Slim over GO:BP, Reactome parents over children) are not
  independent.
