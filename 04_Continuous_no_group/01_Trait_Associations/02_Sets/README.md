# 01_Trait_Associations / 02_Sets

Which gene-set scores track each size and strength trait, with every person pooled?

## The idea

A gene set is a published list of genes that share a pathway or process (Hallmark, Reactome,
KEGG_Legacy, GO BP, or a GO Slim term). `01_Processing/05_Set_Scores` gives each set one number per
biopsy with singscore: the proteins of that biopsy are ranked by abundance, and the set's score is
where its members sit in that ranking, high when they sit near the top. The score rests on ranks
within one biopsy, so it needs no covariate and does not move when other biopsies change.

Each score then meets each trait through the shared association step, `associate_step()` in
`R/associate.R`: the sample model, the change score in six windows, and the growth + strength
and growth x strength models. `../README.md` explains each from the ground up.

Only sets that qualify in `01_Processing/04_Gene_Sets` enter: 15 to 500 members, and at least 15 of
them measured here. Each collection is its own screen, so BH and Storey's q run within collection
as well as within trait and term, or window and trait. A small collection is not corrected
against a large one.

This is a different question from `02_Trait_Pathways`. There a set is tested on how its proteins
rank among all proteins in one trait's t ranking. Here a set's score is a value per person, and
the test asks whether that value follows the trait across people, as for a single protein.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Sets | `readRDS()` of `gene_sets.rds`, `filter(qualifies)` | The set catalogue: `set_id`, `collection`, `pathway`. |
| Values | `readRDS()` of `singscore.rds` | `scores`, one row per qualifying set, one column per biopsy; checked to match the sample sheet's order. |
| Sample sheet | `DAList_corrected.rds` `metadata` | Subject and time point of each biopsy. |
| Labels | `enrichVolcano::ev_clean_label()` | Readable set names; report headers prefix the collection and wrap at 50 characters. |
| Row order | `hclust(as.dist(1 - cor(t(...))), "average")` on singscore, within collection | Orders the overview heatmaps only. |
| Associations | `associate_step(collection_of = ...)` | As in `../README.md`, with BH and q within collection. |

## How to read the outputs

- A set score summarises many proteins, so it can follow a trait when no single member does.
  It can also follow a trait because a few strong members move while the rest do not; check the
  members in `01_Proteins`.
- Sets overlap. Two sets sharing most of their proteins give near-identical scores and will hit
  together; that is one signal.
- The score is relative within the biopsy. A set can rise in rank because other proteins fell.
- BH is within collection, so `fdr` from different collections answers different screens.

## Files

| | |
|---|---|
| Reads | `01_Processing/05_Set_Scores/c_data/singscore.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (catalogue and names), `01_Processing/03_Quantification/c_data/DAList_corrected.rds` (sample sheet only), `00_Input/phenotype.csv` |
| Writes | `c_data/set_trait_associations.xlsx`, `b_reports/set_trait_associations_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Trait_Associations/02_Sets/a_script/set_trait_associations.R` |

```
02_Sets/
├── a_script/set_trait_associations.R
├── b_reports/set_trait_associations_report.pdf
└── c_data/set_trait_associations.xlsx
```

## Workbook: `set_trait_associations.xlsx`

`read_me` lists the sheets; `columns` defines every column. `feature` is the set ID
(`collection|set name`), `collection` its collection and `label` its cleaned name.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | an analysis, comparison and collection | `analysis`, `comparison`, `collection`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `model` | a set, trait and term | `trait`, `term`, `n_samples`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `change_score` | a set, window and trait | `window`, `trait`, `n_subjects`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `rho`, `rho_lower`, `rho_upper`, `p_rho`, `bicor`, `p_bicor`, `fdr`, `qvalue` |
| `two_predictor` | a set, window, 1RM and coefficient | `strength_test`, `coefficient`, `trait`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `interaction` | a set, window and 1RM | `strength_test`, `interaction`, `interaction_lower`, `interaction_upper`, `t`, `p`, `fdr`, `qvalue` |
| `plot_index` | a report page or a drawn hit | `analysis`, `collection`, `pdf_page`, `position` |

`slope` is change in singscore per unit of trait.

## Report: `set_trait_associations_report.pdf`

1. Overviews, one block per collection: a sample-model page (within- and between-person t per
   trait beside the p and BH status maps), then one change-score page per window (Spearman rho
   beside the limma status).
2. Sample-model hits, up to 24 per trait and term across all collections, twelve panels to a page.
3. Change-score hits, up to 24 per window.
4. Two-predictor hits, up to 24 per window.

Panel headers name the collection and set, then give the slope or rho, p, BH q and the number of
biopsies or people.
