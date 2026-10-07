# 02_Trait_Pathways

Which gene sets lead each trait's protein ranking, with every person pooled?

## The idea

`01_Trait_Associations/01_Proteins` gives every protein a limma t for each window and trait: large
and positive when the protein rises with the trait across people, large and negative when it falls.
Sorting the proteins by that t gives one ranking per window and trait.

Gene set enrichment (fgsea) asks whether a set's members bunch toward one end of a ranking rather
than scattering through it. It walks down the ranking and builds a running score that climbs at
each member, by its t, and drops at each non-member; the enrichment score is the furthest that
score gets from zero. The normalised enrichment score (NES) scales it against random sets of the
same size, so sets of different sizes compare. Positive NES: the set sits among proteins that
rise with the trait. Negative NES: among proteins that fall.

The p-value comes from comparing the set with random sets of the same size drawn from the
ranking. That null treats proteins as independent draws. They are not: proteins of one pathway
rise and fall together across people, so their t values move together, and a real set lands at
one end of the ranking more readily than any random set does. fgsea p-values therefore come out
too small, and set counts at FDR < 0.05 run high on any ranking. Every t also rests on the same
16 people or fewer, which the protein-level null does not see. Read which sets lead and in which
direction; do not read the count of sets as evidence.

This stage differs from `01_Trait_Associations/02_Sets`. There each set's singscore is a value per
person tested against the trait across people, which keeps the people as the unit. Here a set is
tested against the other proteins in one ranking.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Rankings | `readRDS()` of `associate.rds`, `change_score` | The limma t per protein, window and trait: six windows, the change traits in `training`, `baseline` and `acute`, the trait levels in the level windows. |
| One protein per gene | `gene_sets.rds` `protein_map`, `selected` | Sets are keyed on gene symbols, so only the representative protein of each symbol enters, and the ranking is named by symbol. |
| Test | `fgsea::fgsea(minSize = 15, maxSize = 500)` | Every qualifying set (Hallmark, Reactome, KEGG_Legacy, GO BP, GO Slim) on each ranking. `set.seed(1)` before each run, so reruns match. |
| Multiple testing | fgsea's `padj` (BH) | Within one window and trait, across all sets and collections together. |
| Counts | `summarise()` | Per window and trait: sets tested, at p < 0.05, at FDR < 0.05, and up or down by the sign of NES. |

## How to read the outputs

- NES says direction and how tightly the set bunches, not how much its proteins change.
- The set counts run high for the reason above. A trait with many sets at FDR < 0.05 is not, on
  that basis, more linked to the proteome than one with few.
- Sets overlap heavily (a Reactome pathway and its parent share most proteins). Several hits can
  be one signal.
- `baseline` and `baseline_level` rank the same T1 values against different traits; `training`,
  `baseline` and `acute` rank different values against the same traits.
- The leading edge (the members that carry the score) is kept in `pathways.rds`, not in the
  workbook.
- No sample-level test (such as fry, which keeps the correlation between proteins) runs here.

## Files

| | |
|---|---|
| Reads | `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/associate.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (sets, catalogue, protein-to-gene map) |
| Writes | `c_data/trait_pathways.xlsx`, `c_data/pathways.rds` (the `pathways` table with `leading_edge`), `b_reports/trait_pathways_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/02_Trait_Pathways/a_script/trait_pathways.R` |

```
02_Trait_Pathways/
├── a_script/trait_pathways.R
├── b_reports/trait_pathways_report.pdf
└── c_data/
    ├── trait_pathways.xlsx
    └── pathways.rds
```

## Workbook: `trait_pathways.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | a window and trait | `n_tested`, `n_nominal`, `n_fdr_05`, `n_up`, `n_down` |
| `pathways` | a set in a window and trait | `set_id`, `collection`, `pathway`, `size` (members in the ranking), `nes`, `p`, `fdr` |

## Report: `trait_pathways_report.pdf`

1. Sets at FDR < 0.05 per trait, up and down, one facet per window.
2. One page per window that has any set at FDR < 0.05: up to ten such sets per trait, smallest p
   first, as NES bars.
