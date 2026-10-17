# S3: pathways, full selection

Four panels, drawn with the shared pathway dot-matrix (`style/pathway_grid.R`). Together they run
well past 240 mm, so `S3_composite.pdf` has one page per panel, each page as tall as its panel.

| Panel | Shows |
|---|---|
| a_tree_matrix | Every candidate set per block from `02_HRvLR/02_Set_Contrasts`, in its tree's leaf order, in all nine contrasts. Grey and white bands alternate between clusters at height 0.85; each cluster's most specific set is in bold (F3 draws the first eight per block). About 1.2 m tall: it lists all 434 candidates. |
| b_cut_075 | One set per cluster with the tree cut at 0.75 (tighter clusters, more rows), every cluster of every block, the interaction window included. |
| c_cut_09 | The same at 0.9 (looser clusters, fewer rows). |
| d_go_slim | Every GO Slim theme in the nine contrasts, rows ordered by average-linkage clustering of their NES. GO Slim is tested upstream but kept out of the F3 selection. |

Fill is fgsea NES, size -log10 FDR, the square the database; `*` marks a set whose FDR < 0.05 calls
in its window disagree in sign, `†` a training or acute set whose interaction also reaches
FDR < 0.05. A cell is empty where fgsea returned no score.

## Inputs

| File | Used for |
|---|---|
| `02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds` | `set_tests` (fgsea NES and FDR, GO Slim included), `selection` (rank, tree order, cuts, representatives) |

## Run

```sh
Rscript 04_Figures/S3/a_script/S3_composite.R              # all pages and workbook
Rscript 04_Figures/S3/a_script/panels/a_tree_matrix.R      # one panel
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/{a_tree_matrix,b_cut_075,c_cut_09,d_go_slim}.pdf` | one PDF per panel |
| `b_reports/S3_composite.pdf` | four pages, 180 mm wide, one panel each with its letter |
| `c_data/S3.xlsx` | one sheet per panel with every drawn cell |
