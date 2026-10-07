# S02_NES_Scatter

Supplement: do HR and LR shift the same gene sets, set by set, in training and in the bout?

## The idea

fgsea gives every gene set one normalised enrichment score (NES) per contrast: positive when the
set sits toward the up end of the protein ranking, negative toward the down end. Plotting one
group's NES against the other's for the same window puts the two groups' set-level responses on
one pair of axes. A set on the dashed identity line moved the same way and as far in both; a set
in the top-right or bottom-left quadrant moved the same way; a set in the other two quadrants
moved in opposite directions.

Two pairs are drawn: `Training_HR` against `Training_LR`, and `Acute_HR` against `Acute_LR`. Only
sets scored in both contrasts of a pair are compared (a set near the 15-protein floor can drop out
of one).

Terms used on the pages and in the workbook:

| Term | Meaning |
|---|---|
| significance | where the set reaches fgsea FDR < 0.05: both, one contrast, or neither (NS) |
| discordant | the two NES have opposite signs; the discordant panel and counts take only sets significant in at least one group |
| survivor | kept as non-redundant by `collapsePathways()` in either contrast |
| concordance rho | Spearman correlation of the two NES across the sets in a population |

Spearman's rho is the correlation of ranks: +1 when the sets order identically in both groups, 0
when the orders are unrelated.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Pairing | `tidyr::pivot_wider()` on `set_tests`, fgsea rows | One row per set with both NES, both FDR and both `main` flags. Point size is the larger of the two measured sizes. |
| Concordance | `cor.test(method = "spearman", exact = FALSE)` | Per population: all collections, collapse survivors, and Hallmark with GO Slim. Needs at least three sets; a page panel prints rho only with 30 or more. |
| Composite 1 | `patchwork::wrap_plots()` | All collections; collapse survivors; the discordant significant sets on their own axes, each named. |
| Composite 2 | as above | Hallmark and GO Slim, the two collections whose members do not nest, with the eight significant sets of smallest summed FDR named; then the significant same-sign sets split into "Up in both" and "Down in both", each rescaled so every set can be named. |
| Output | `pdf()` | 11 x 8.5 in; two composites per pair, four pages. |

## How to read the outputs

- Sets overlap and share proteins, so they are not independent points. The concordance rho
  describes agreement; its p treats sets as independent and overstates the evidence. The Hallmark
  and GO Slim population limits the nesting but does not remove overlap.
- NES compares a set with the other proteins in the same contrast. If one group's whole proteome
  moved more, its NES can differ from the other's without the set itself behaving differently.
- Agreement between the two groups' NES is not a test of a group difference; the interaction
  contrasts in `03_Set_Contrasts` are.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` (`set_tests`) |
| Writes | `b_reports/s02_nes_scatter.pdf`, `c_data/s02_nes_scatter.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/S02_NES_Scatter/a_script/s02_nes_scatter.R` |

```
S02_NES_Scatter/
├── a_script/s02_nes_scatter.R
├── b_reports/s02_nes_scatter.pdf
└── c_data/s02_nes_scatter.xlsx
```

## Workbook: `s02_nes_scatter.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `concordance` | a pair and population | `pair`, `population`, `n_sets`, `n_significant`, `n_discordant`, `rho`, `p`, `same_sign` |
| `nes_scatter` | a set scored in both contrasts of a pair | `pair`, `x_contrast`, `y_contrast`, `collection`, `pathway`, `n_proteins`, `nes_x`, `nes_y`, `fdr_x`, `fdr_y`, `significance`, `discordant`, `survivor` |
