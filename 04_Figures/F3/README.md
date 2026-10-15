# F3: pathways

One panel, drawn with the shared pathway dot-matrix (`style/pathway_grid.R`).

| Panel | Shows |
|---|---|
| a_pathways | Per window (Training, Acute, Between groups, Interaction), the sets `02_HRvLR/02_Set_Contrasts` selects as candidates, grouped by a tree of shared proteins and cut at Jaccard height 0.85; one set per cluster (`representative_0.85`), up to ten per window by rank. Each set is drawn in all nine contrasts: fill fgsea NES, size -log10 FDR, square the database. `*` marks a set whose FDR < 0.05 calls in its window disagree in sign. |

Disease sets are removed upstream in `01_Processing/04_Gene_Sets`; GO Slim is left out of the
selection (see S3d).

## Inputs

| File | Used for |
|---|---|
| `02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds` | `set_tests` (fgsea NES and FDR), `selection` (rank, tree cuts, representatives) |

## Run

```sh
Rscript 04_Figures/F3/a_script/F3_composite.R            # figure and workbook
Rscript 04_Figures/F3/a_script/panels/a_pathways.R       # the panel alone
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/a_pathways.pdf` | the panel |
| `b_reports/F3_composite.pdf` | the figure, 180 mm wide |
| `c_data/F3.xlsx` | sheet `a_pathways`: every drawn cell (window, rank, set, contrast, NES, FDR) |
