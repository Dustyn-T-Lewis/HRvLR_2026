# 05 · Integration

Works across the protein, pathway and module levels.

| Step | Runs | Writes |
|---|---|---|
| [`01_Chance_Line`](01_Chance_Line/README.md) | chance line for the candidate counts: 200 shuffles of the subjects' phenotype rows | workbook, PDF |
| [`02_Candidates`](02_Candidates/README.md) | every nominal hit on the HR-against-LR tests, per level, training and acute | workbook |
| [`03_Panel`](03_Panel/README.md) | elastic-net panels and DIABLO, with permutation nulls | workbook, PDF |
| [`04_Continuous`](04_Continuous/README.md) | trait-ranked pathways, sPLS, module membership, STRING, the composition check, all against phenotype shuffles | workbook, PDF, `shuffled_t.rds` |
| [`05_Clustering`](05_Clustering/README.md) | proteins by association profile, per-trait groups, people, within-pathway splits, biclusters, each against its shuffle null | workbook, PDF |

`04_Continuous` saves the protein t from 100 phenotype shuffles; `05_Clustering` reads it, so the
two share one null.
