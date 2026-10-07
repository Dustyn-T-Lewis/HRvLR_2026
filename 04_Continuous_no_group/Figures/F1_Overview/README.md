# F1_Overview

Figure 1 for the continuous arc: how people responded, how many features track each trait, and
which sets and protein groups lead, on one page.

## The idea

The figure gathers the arc's saved outputs into five panels. It reads workbooks and RDS files
written by the stages and recomputes no test. Its one computation is a proteome PCA, a projection
of the biopsies onto the two directions of greatest spread, which carries no p-value.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Each person's training response: % change in each size and strength trait and the composite, scaled to a z score within the trait. People run left to right by composite hypertrophy. Grey: not measured. Volume is left out. | `phenotype.csv`, through `phenotype_tables()` |
| b | Share of features at p < 0.05 (limma, one value per person) for each trait in the `training`, `baseline` and `acute` windows, one facet per level (proteins, gene sets, modules). The colour scale is white at 5%, the share chance gives; red above, blue below. For gene sets the counts are summed over collections. | `counts` sheet of each level's workbook |
| c | Proteome PCA: every biopsy on PC1 and PC2 of the complete blood-corrected matrix, proteins scaled, shape by time point, coloured by the person's composite hypertrophy (left) and leg-press gain (right). All of one person's biopsies share a colour. | `DAList_corrected.rds`, `prcomp(scale. = TRUE)` |
| d | Gene sets leading each trait's protein ranking over training: fgsea NES as fill, -log10 FDR (capped at 8) as dot size, each selected set shown in every trait. | `trait_pathways.xlsx`, `pathways` |
| e | Profile groups: each group's mean limma t per window and trait, rows labelled with group number, size and main module. | `profile_clusters.xlsx`, `profile_means` and `profile_groups` |

Panel d selects its rows in four steps: training window only, GO Slim left out; keep sets at
FDR < 0.05 for at least two traits; rank them by that count, then by total |NES| across traits;
drop any set that shares half or more of its measured proteins with a stronger one
(`prune_overlapping()` in `R/figure_panels.R`) and keep the first 20.

Before writing, the script checks that every level, window and trait in panel b tested at least
one feature, that panel d has one row per selected set and trait, and that no share exceeds 1.

## How to read the figure

- Panel b compares each share with 5%. A share near 5% is what chance gives with no link at all;
  a share well above it says the trait has more links than chance, not which ones are real.
- Panel c tests nothing. Colour shows whether growth or strength lines up with the main axes of
  the proteome; the absence of a gradient does not rule out a link along a smaller axis.
- Panel d inherits fgsea's high counts (see `02_Trait_Pathways/README.md`). Read it for which sets
  lead and their direction.
- Panel e describes groups built from the t values it shows; see `03_Profile_Clusters/README.md`.

## Files

| | |
|---|---|
| Reads | `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/protein_trait_associations.xlsx`, `04_Continuous_no_group/01_Trait_Associations/02_Sets/c_data/set_trait_associations.xlsx`, `04_Continuous_no_group/01_Trait_Associations/03_Modules/c_data/module_trait_associations.xlsx` (each `counts`), `04_Continuous_no_group/02_Trait_Pathways/c_data/trait_pathways.xlsx` (`pathways`), `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` (set members for pruning), `04_Continuous_no_group/03_Profile_Clusters/c_data/profile_clusters.xlsx` (`profile_means`, `profile_groups`), `01_Trunk/03_Quantification/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `b_reports/f1_overview.pdf` (300 x 175 mm, `cairo_pdf` through `save_figure()`), `c_data/f1_overview.xlsx` |
| Run | `Rscript 04_Continuous_no_group/Figures/F1_Overview/a_script/f1_overview.R`, after every stage of the arc. `run_all.R` does not list it. |

```
F1_Overview/
├── a_script/f1_overview.R
├── b_reports/f1_overview.pdf
└── c_data/f1_overview.xlsx
```

Theme, colours and thresholds (`alpha_nominal`, `alpha_fdr`, both 0.05) come from
`R/figure_theme.R`; `trait_label()` turns `pct_leg_press_1rm` into "leg press 1RM (% change)".
`pathway_label()` and `prune_overlapping()` come from `R/figure_panels.R`.

## Workbook: `f1_overview.xlsx`

`read_me` lists the sheets; `columns` defines every column. Each sheet holds the numbers one panel
draws. `trait` and, in `links`, `window` hold the labels the panels print ("leg press 1RM
(% change)", "T1 level"), not the column codes.

| Sheet | Panel | One row is | Key columns |
|---|---|---|---|
| `people` | a | a person and trait | `subject`, `trait`, `value` (% change), `z` |
| `links` | b | a level, window and trait | `level`, `window`, `trait`, `n_tested`, `n_nominal`, `share` |
| `pca` | c | a biopsy and colouring trait | `sample_id`, `subject`, `timepoint`, `PC1`, `PC2`, `trait`, `value` |
| `sets` | d | a selected set and trait | `set_id`, `collection`, `pathway`, `trait`, `nes`, `fdr` |
| `profiles` | e | a profile group and screen | `group`, `n_proteins`, `main_module`, `screen`, `mean_t` |
