# S3: pathways, full selection

Five panels, drawn with the shared pathway dot-matrix (`style/pathway_grid.R`). Together they run
well past one page, so `S3_composite.pdf` has a page per panel and panel a a page per window (four
pages), each 250 mm wide and as tall as its content.

| Panel | Shows |
|---|---|
| a_tree_matrix | Every candidate set per block from `02_HRvLR/02_Set_Contrasts`, in its tree's leaf order, in every contrast, one page per window (Training, Acute, Between groups, Interaction). Grey and white bands alternate between clusters at height 0.85; each cluster's most specific set is in bold (F3 draws the first few per block). It lists all 434 candidates. |
| b_cut_075 | One set per cluster with the tree cut at 0.75 (tighter clusters, more rows), every cluster of every block in all four windows. |
| c_cut_09 | The same at 0.9 (looser clusters, fewer rows). |
| d_go_slim | Every GO Slim theme in every contrast, rows ordered by average-linkage clustering of their NES. GO Slim is tested upstream but kept out of the F3 selection. |
| e_hallmark_heatmap | The F5 a heatmap for the 35 Hallmark sets with 15 to 500 measured genes, as singscore set scores: limma contrasts on the scores (`score_tests` in `02_Set_Contrasts`) and the trait links of `07_Group_Trait_Slopes`. Hallmark is chosen in advance, so the rows are not selected by the tests drawn; rows are clustered by their 36 values. |

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
| `b_reports/S3_composite.pdf` | eight pages, 250 mm wide: panel a by window, then b, c, d and e |
| `c_data/S3.xlsx` | one sheet per panel with every drawn cell |
