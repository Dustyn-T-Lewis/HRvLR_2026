# F1_Overview

Figure 1: what does the design look like, and how does the whole proteome respond before any
single result? Every p-value it draws comes from `02_Protein_Contrasts` or `03_Set_Contrasts`; its
only computation is the phenotype PCA, a projection with no test.

| Panel | Shows | Source |
|---|---|---|
| a | The design: measurements at each time point and the four contrast families, drawn in BioRender. | `assets/design.png` |
| b | Proteins per contrast, down to the left and up to the right. The light bar counts proteins at p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives proteins tested and the count at FDR < 0.05. | `set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and phenotype PCA (right), each with a shaded 95% ellipse per time point and its PERMANOVA. | `protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv` |
| d | Gene sets at fgsea FDR < 0.05 shared within each family of contrasts, each shown in all nine contrasts. | `set_tests.rds`, `set_tests` |
| e, f | Each group's p < 0.05 proteins for one window (e training, f acute) placed in the other group's ranking for that window, as limma's barcode plot. | rankings from `protein_results`; `protein_contrasts.xlsx`: `fry_concordance` |

| | |
|---|---|
| Reads | `03_Group_comparisons/HRvLR/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`contrast_summary`, `pca`, `pca_variance`, `permanova`, `fry_concordance`), `03_Group_comparisons/HRvLR/03_Set_Contrasts/c_data/set_tests.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (set members, for pruning), `00_Setup/Input/phenotype.csv`, `assets/design.png` |
| Writes | `b_reports/f1_overview.pdf`, `c_data/f1_overview.xlsx` |
| Run | `Rscript 05_Figures/HRvLR/F1_Overview/a_script/f1_overview.R` |

## Reading the panels

- Π (panel b) is p raised to the power \|logFC\| (Xiao et al. 2014, in p form): it rewards a small
  p and a large fold change together. It ranks and selects; it is not a test.
- Panel c: the phenotype PCA uses the six size and strength traits on every T1 and T2 biopsy with
  all six measured. The proteome gets a time and a group PERMANOVA; the phenotype a time test only,
  because HR and LR were split on how these traits changed.
- Panel d picks rows by family: training and acute sets at fgsea FDR < 0.05 in both groups' contrasts,
  HR - LR sets at two of three time points, and the strongest five interaction sets. A set sharing
  half its proteins with a stronger one is dropped; `*` marks a set whose sign differs between its
  significant contrasts.
- Panels e and f: one group's p < 0.05 proteins for a window, ticked at their rank in the other
  group's moderated t for the same window, with fry's FDR from `fry_concordance`. The groups are
  different people, so the set and the ranking share no cell mean.
- The script stops if a plotted count, PERMANOVA row or barcode FDR differs from its source table.

## Read with care

- Read the floor (`Baseline_HRvLR`) in panel b before any other count.
- fgsea in panel d treats proteins as independent draws, so a large co-regulated set is easy to
  call; `03_Set_Contrasts` reports fry beside it.
- Carry-over in one direction only (e, f) can come from unequal power as well as biology.
- Panel a is a BioRender drawing (`assets/design.png`) and reads no data.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `protein_counts` | a contrast and direction (panel b) | `group`, `label`, `pi`, `nominal`, `fdr` |
| `pca_proteome` | a biopsy (panel c, left) | `PC1`, `PC2` |
| `pca_phenotype` | a biopsy with all six traits (panel c, right) | `subject`, `arm`, `timepoint`, `PC1`, `PC2` |
| `permanova` | a test (panel c) | `data`, `term`, `window`, `arms`, `r2`, `f`, `p` |
| `pathways` | a selected set in a contrast (panel d) | `selected_in`, `window`, `column`, `collection`, `label`, `nes`, `fdr` |
| `barcode` | a lane and direction (panels e, f) | `window`, `chosen`, `ranked_in`, `direction`, `n`, `same_sign_share`, `fdr` |
