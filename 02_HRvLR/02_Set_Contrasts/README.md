# 02_Set_Contrasts

Tests every gene set in the nine contrasts with fgsea and fry, FDR within each database, and picks
the sets the pathway figure draws: candidates per window, grouped by shared proteins into a tree,
one set per cluster.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | sets (disease sets excluded), catalogue, one protein per gene |
| `02_HRvLR/01_Protein_Contrasts/c_data/fit.rds` | the protein fit and contrasts |

## Run

```sh
quarto render 02_HRvLR/02_Set_Contrasts/a_script/set_contrasts.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/set_tests.rds` | `set_tests` (every test, with leading edges) and `selection` (candidates, tree clusters, representatives) |
| `c_data/set_contrasts.xlsx` | the same two tables |
| `b_reports/set_contrasts.html` | the rendered walkthrough |

## Read with care

- fgsea ranks proteins against each other; fry tests whether a set's proteins move together across
  biopsies. A set can pass one and not the other.
- The tree groups sets sharing at least about 15% of proteins (Jaccard height 0.85); a cluster's
  representative is its top-ranked set, not a merged pathway.
- GO Slim is tested but left out of the figure selection, since its terms repeat GO:BP.
- Training_LR protein t values lean positive (median 0.65), so fgsea cannot estimate p for 11
  sets; their p, NES and FDR stay missing and the report lists them. fry tests every set.
