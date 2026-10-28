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
| `02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds` | `set_tests` (fgsea NES and FDR), `selection` (rank, tree cuts, representatives) |

## Run

```sh
Rscript 04_Figures/F3/a_script/F3_composite.R            # figure and workbook
Rscript 04_Figures/F3/a_script/panels/a_pathways.R       # the panel alone
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel alone |
| `b_reports/F3_composite.pdf` | the figure, 250 mm wide landscape |
| `c_data/F3.xlsx` | sheet `a_pathways`: every drawn cell (window, rank, set, contrast, NES, FDR) |
