# F1_Overview

Figure 1: what does the design look like, and how does the whole proteome respond before any
single result?

## The idea

The figure moves from the design to the whole proteome, then to gene sets, then to whether one
group's response carries into the other. Every p-value it
draws comes from `02_Protein_Contrasts` or `03_Set_Contrasts`; its only computation is the
phenotype PCA, a projection with no test.

| Panel | Shows | Source |
|---|---|---|
| a | The design: measurements at each time point and the four contrast families, drawn in BioRender. | `assets/design.png` |
| b | Proteins per contrast, down to the left and up to the right. The light bar counts proteins at p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives proteins tested and the count at FDR < 0.05. | `set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and phenotype PCA (right), each with a shaded 95% ellipse per time point and its PERMANOVA. | `protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv` |
| d | Gene sets at fgsea FDR < 0.05 shared within each family of contrasts, each shown in all nine contrasts. | `set_tests.rds`, `set_tests` |
| e, f | Each group's p < 0.05 proteins for one window (e training, f acute) placed in the other group's ranking for that window, as limma's barcode plot. | rankings from `protein_results`; `protein_contrasts.xlsx`: `fry_concordance` |

Π (panel b) is p raised to the power \|logFC\|. It is Xiao et al.'s (2014) score, π = \|logFC\| x
-log10 p, turned back into p form: Π = 10^-π, so Π < 0.05 means π > 1.3. It rewards a small p and
a large fold change together. Here it ranks and selects; it is not a test.

PCA (panel c) rotates the data onto new axes ordered by how much variation each carries; PC1 and
PC2 are the two largest, and the axis label gives each one's share. Points that sit together have
similar profiles. The proteome PCA uses the complete blood-corrected matrix with proteins scaled.
The phenotype PCA uses the six size and strength traits (fCSA mixed, type I, type II, mCSA, leg
press 1RM, leg extension 1RM), scaled, on every T1 and T2 biopsy with all six measured; a biopsy
missing one is left out, not filled. PERMANOVA splits the distances between biopsies into the
share a term explains (R²) and gets its p from permutations. The proteome gets a time test and a
group test (on one mean profile per person); the phenotype gets a time test only, because HR and
LR were split on how much these traits changed.

The pathway grid (panel d) picks rows by family, from Hallmark, Reactome, KEGG and GO BP:

| Family | A set is shown when it reaches fgsea FDR < 0.05 in |
|---|---|
| Training | both groups' training contrasts |
| Acute | both groups' acute contrasts |
| HR - LR | at least two of the three time points |
| Interaction | either interaction; strongest five only |

Within a family, sets with one sign across their significant contrasts come first, then by total
\|NES\| across the family's contrasts. A set sharing half or more of its measured proteins with a
stronger set already shown is dropped (`prune_overlapping()`). `*` marks a set whose sign differs
among the contrasts where it is significant. A set is labelled "Proteasome:" when half or more of
its significant leading-edge proteins are proteasome subunits (gene symbols starting PSM), because
several Reactome signalling sets carry only that signal. Dot fill is NES, dot size -log10 FDR
(capped at 12), and the square beside each name marks its database.

The barcode (panels e, f) asks whether one group's response carries over to the other. For a
window, take the proteins at p < 0.05 in one group's contrast, split into up and down. Rank every
protein by its moderated t in the other group's contrast for the same window, top (up) to bottom
(down). Each member is a tick at its rank, the up set above the axis and the down set below; the
curve is a tricube moving average of membership along the rank (`limma::tricubeMovingAverage()`),
so a bulge at one end shows where members gather. Each lane's label gives the set size, the share
of members moving the same way in both groups, and the fry FDR from `02_Protein_Contrasts`. The two
groups are different people, so set and ranking share no cell mean, and the fry test is not
reusing the data that chose the set. fry's FDR there runs over the sets tested together in one
call (Π, BH and nominal selections, up and down, each with at least five members).

## What the code does

| Step | Function | Detail |
|---|---|---|
| Panel a | `png::readPNG()`, `grid::rasterGrob()` | The BioRender export, pinned under its title. |
| Panel b | `panel_protein_counts()` in `R/figure_overview.R` | Counts per contrast and direction from `protein_results`. |
| Panel c | `panel_pca_pair()`, `pca_panel()`; `prcomp(scale(...))` for the phenotype | Ellipses are `stat_ellipse(type = "t", level = 0.95)`. |
| Panel d | `panel_pathway_grid()`, `prune_overlapping()`, `pathway_label()` | Rule as above. |
| Panels e, f | `signature_barcode()` in `R/figure_panels.R` | `span = 0.12`; fry numbers read from `fry_concordance`. |
| Assembly | `figure_overview()`, `patchwork`, `save_figure()` | 300 x 175 mm, `cairo_pdf`. |
| Checks | `stopifnot()` | FDR counts per contrast equal `contrast_summary`; one PERMANOVA row per test; every pathway row has all nine contrasts; all eight barcode FDRs equal `fry_concordance`. |

## How to read the outputs

- Panel b counts are per contrast. With seven or eight people per group, read the count at
  `Baseline_HRvLR` (the floor) before reading any other.
- fgsea in panel d ranks a set against the other proteins and treats proteins as independent
  draws, so a large co-regulated set is easy to call. `03_Set_Contrasts` reports fry beside it.
- In e and f, carry-over in one direction and not the other can come from unequal power (one
  group's set is larger, or its ranking noisier), not only from biology. A time-linked technical
  effect would also move both groups together, and the design cannot separate it from shared
  biology.
- Panel a reads no data. To change it, edit it in BioRender, export a PNG to `assets/design.png`
  and rerun the script.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`contrast_summary`, `pca`, `pca_variance`, `permanova`, `fry_concordance`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (set members, for pruning), `00_Input/phenotype.csv`, `assets/design.png` |
| Writes | `b_reports/f1_overview.pdf`, `c_data/f1_overview.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F1_Overview/a_script/f1_overview.R` |

```
F1_Overview/
├── a_script/f1_overview.R
├── assets/design.png
├── b_reports/f1_overview.pdf
└── c_data/f1_overview.xlsx
```

## Workbook: `f1_overview.xlsx`

`read_me` lists the sheets and carries a note that panel a reads no data; `columns` defines every
column; `input_manifest` and `package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `protein_counts` | a contrast and direction (panel b) | `group`, `label`, `pi`, `nominal`, `fdr` |
| `pca_proteome` | a biopsy (panel c, left) | `PC1`, `PC2` |
| `pca_phenotype` | a biopsy with all six traits (panel c, right) | `subject`, `arm`, `timepoint`, `PC1`, `PC2` |
| `permanova` | a test (panel c) | `data`, `term`, `window`, `arms`, `r2`, `f`, `p` |
| `pathways` | a selected set in a contrast (panel d) | `selected_in`, `window`, `column`, `collection`, `label`, `nes`, `fdr` |
| `barcode` | a lane and direction (panels e, f) | `window`, `chosen`, `ranked_in`, `direction`, `n`, `same_sign_share`, `fdr` |
