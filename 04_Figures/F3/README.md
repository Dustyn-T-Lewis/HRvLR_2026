# F3: pathways

The pathway dot-matrix (`style/pathway_grid.R`). The gene-set counts behind it are F1 e.

| Panel | Shows |
|---|---|
| a_pathways | Seven row blocks: Training and Acute each split into shared, HR only and LR only (fgsea FDR < 0.05 in both groups, HR alone or LR alone; candidate count in the strip), then Between groups (two of three time points). Within a block, sets are grouped by a tree of shared proteins cut at Jaccard height 0.85, each cluster shown by its most specific set (`representative_0.85`), up to eight per block by rank. Each set is drawn in all nine contrasts: fill fgsea NES, size -log10 FDR, square the database. `*` marks calls of opposite sign within the window, `†` a set whose interaction also reaches FDR < 0.05. The interaction window is in S3. |

Disease sets are removed upstream in `01_Processing/04_Gene_Sets`; GO Slim is left out of the
selection and the counts (see S3d).

## Read with care


- "HR only" means fgsea FDR < 0.05 in HR and not in LR; it is not a test of HR against LR. The
  `†` sets are the ones the interaction also supports (22 training, 53 acute sets in F1 e).
- Sets overlap: one mitochondrial programme reaches many sets, which is why panel a shows one set
  per tree cluster and F1 e's counts are not independent programmes.
- The T3 HR − LR count (263 sets) is the largest in F1 e, while T1 has 65 and T2 none; panel a's
  between-group block needs two of the three time points.

## Inputs


| File | Used for |
|---|---|
| `02_HRvLR/02_Set_Contrasts/c_data/set_contrasts.rds` | `set_tests` (fgsea NES and FDR), `selection` (rank, tree cuts, representatives) |

## Supplementary figure S3: pathways, full selection

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

### Inputs


| File | Used for |
|---|---|
| `02_HRvLR/02_Set_Contrasts/c_data/set_contrasts.rds` | `set_tests` (fgsea NES and FDR, GO Slim included), `selection` (rank, tree order, cuts, representatives) |

## Run

```sh
Rscript run_all.R 04_Figures/F3                       # F3 and S3
Rscript 04_Figures/F3/a_script/F3_composite.R        # F3 alone
```

## Outputs

| File | Holds |
|---|---|
| `a_script/`, `b_reports/` | F3's panel scripts and PDFs (`a_*.R`, `a_*.pdf`, ...) and its composite (`F3_composite.pdf`: the figure, 250 mm wide landscape) |
| `a_script/supp/`, `b_reports/supp/` | S3: `S3_*.R`, `S3_*.pdf` and `S3_composite.pdf` (eight pages, 250 mm wide: panel a by window, then b, c, d and e) |
| `c_data/F3.xlsx`, `c_data/S3.xlsx` | one sheet per panel, the numbers drawn |
