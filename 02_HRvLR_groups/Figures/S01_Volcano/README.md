# S01_Volcano

Supplement: which proteins move in each within-group and interaction contrast, and which fgsea
sets do they belong to?

## The idea

A volcano plot puts every protein at its logFC (x: left down, right up) against -log10 p (y:
higher is smaller p). Proteins that change a lot with strong evidence sit in the upper corners.
Here point colour marks proteins at BH FDR < 0.05, and a badge counts them.

Around the volcano, a ring carries up to eight gene sets: the fgsea sets at FDR < 0.05 that
survived `collapsePathways()` (non-redundant), the eight with the smallest FDR. Each set is an arc
coloured by NES (positive: the set sits toward the up end of the ranking), with ticks running to
its leading-edge proteins, the members driving its score. The panel therefore carries two claims
at two levels: point colour and badges read protein FDR; the ring reads set FDR.

Eight volcanoes, in this order:

| Page | Contrast | Labels |
|---|---|---|
| 1 | `Training_Interaction` | up to five proteins at FDR < 0.05, smallest FDR first |
| 2 | `Acute_Interaction` | as above |
| 3 to 6 | `Training_HR`, `Training_LR`, `Acute_HR`, `Acute_LR` | as above |
| 7, 8 | `Training_HR`, `Training_LR` again | the five proteins with the smallest Π |

The interaction pages carry their algebra as a subtitle, for example
(HR_T2 - HR_T1) - (LR_T2 - LR_T1). The between-group level contrasts (`Baseline_HRvLR`,
`Trained_HRvLR`, `Acute_HRvLR`) are not drawn; `03_Set_Contrasts` reports them.

Π is p raised to the power \|logFC\| (Xiao et al. 2014), which rewards a small p and a large fold
change together. Pages 7 and 8 use it to choose labels only: same points, colours and rings as
pages 3 and 4, and a protein named there is not a hit.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Points | `protein_results` from `set_tests.rds` | Untested proteins dropped; p bounded at the smallest double for the log only. |
| Ring | `set_tests`, `method == "fgsea"`, `main`, FDR < 0.05 | `slice_min(fdr, n = 8)`. Two sets whose names match after the collection prefix (GO and Reactome muscle contraction, for example) get the collection appended so their arcs differ. |
| Leading edges | gene symbol to point label | Edges are translated into the label space the points use, so every tick finds its protein. |
| Drawing | `enrichVolcano::volcano_ring()`, `volcano_ring_theme()` | `p_threshold = 0.05` on FDR, `logfc_threshold = 0`; default palette (red up, blue down). |
| Output | `pdf()` | 11 x 8.5 in, one volcano per page. |

## How to read the outputs

- A protein's colour and a set's arc are separate tests. A protein under a ringed set's tick need
  not be significant itself, and a coloured protein need not sit in any ringed set.
- fgsea ranks against the other proteins and treats proteins as independent draws; a ringed set
  is a ranking result, not a fry result.
- A volcano for one group's training and one for the other's, side by side, do not compare the
  groups. The interaction pages do.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` (`protein_results`, `set_tests`) |
| Writes | `b_reports/s01_volcano.pdf`, `c_data/s01_volcano.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/S01_Volcano/a_script/s01_volcano.R` |

```
S01_Volcano/
├── a_script/s01_volcano.R
├── b_reports/s01_volcano.pdf
└── c_data/s01_volcano.xlsx
```

## Workbook: `s01_volcano.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook. The points themselves are `protein_results` in
`set_tests.rds`.

| Sheet | One row is | Key columns |
|---|---|---|
| `ringed_sets` | a set ringed on one volcano | `contrast`, `rank_by`, `set_id`, `collection`, `pathway`, `nes`, `fdr` |
| `labelled_proteins` | a protein named on one volcano | `contrast`, `rank_by` (fdr or pi), `label` |
