# F3: pathways

One panel, drawn with the shared pathway dot-matrix (`style/pathway_grid.R`).

| Panel | Shows |
|---|---|
| a_pathways | Seven row blocks: Training and Acute each split into shared, HR only and LR only (fgsea FDR < 0.05 in both groups, HR alone or LR alone; candidate count in the strip), then Between groups (two of three time points). Within a block, sets are grouped by a tree of shared proteins cut at Jaccard height 0.85, each cluster shown by its most specific set (`representative_0.85`), up to eight per block by rank. Each set is drawn in all nine contrasts: fill fgsea NES, size -log10 FDR, square the database. `*` marks calls of opposite sign within the window, `†` a set whose interaction also reaches FDR < 0.05. The interaction window is in S3. |

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
