# F1_Overview

Figure 1 for the continuous analysis: how people responded, how many features track each trait, and which sets and protein groups lead. It reads the folder's saved outputs; its one computation is a proteome PCA, which carries no p-value.

| Panel | Shows | Source |
|---|---|---|
| a | Each person's training response: % change in each size and strength trait and the composite, scaled to a z score within the trait. People run left to right by composite hypertrophy. Grey: not measured. Volume is left out. | `phenotype.csv`, through `phenotype_tables()` |
| b | Share of features at p < 0.05 (limma, one value per person) for each trait in the `training`, `baseline` and `acute` windows, one facet per level (proteins, gene sets, modules). The colour scale is white at 5%, the share chance gives; red above, blue below. For gene sets the counts are summed over collections. | `counts` sheet of each level's workbook |
| c | Proteome PCA: every biopsy on PC1 and PC2 of the complete blood-corrected matrix, proteins scaled, shape by time point, coloured by the person's composite hypertrophy (left) and leg-press gain (right). All of one person's biopsies share a colour. | `DAList_corrected.rds`, `prcomp(scale. = TRUE)` |
| d | Gene sets leading each trait's protein ranking over training: fgsea NES as fill, -log10 FDR (capped at 8) as dot size, each selected set shown in every trait. | `trait_pathways.xlsx`, `pathways` |
| e | Profile groups: each group's mean limma t per window and trait, rows labelled with group number, size and main module. | `profile_clusters.xlsx`, `profile_means` and `profile_groups` |

| | |
|---|---|
| Reads | `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/protein_trait_associations.xlsx`, `04_Continuous_no_group/01_Trait_Associations/02_Sets/c_data/set_trait_associations.xlsx`, `04_Continuous_no_group/01_Trait_Associations/03_Modules/c_data/module_trait_associations.xlsx` (each `counts`), `04_Continuous_no_group/02_Trait_Pathways/c_data/trait_pathways.xlsx` (`pathways`), `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` (set members for pruning), `04_Continuous_no_group/03_Profile_Clusters/c_data/profile_clusters.xlsx` (`profile_means`, `profile_groups`), `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | `b_reports/f1_overview.pdf` (300 x 175 mm, `cairo_pdf` through `save_figure()`), `c_data/f1_overview.xlsx` |
| Run | `Rscript 05_Figures/Continuous/F1_Overview/a_script/f1_overview.R`, after every step of the folder; `run_all.R` runs it after `03_Profile_Clusters`. |

## Reading the panels

- Panel b compares each share with 5%, what chance gives with no link; a share well above it says the trait has more links than chance, not which ones are real.
- Panel d keeps training-window sets at FDR < 0.05 for at least two traits (GO Slim out), ranked by that count and total |NES|, drops any set sharing half its proteins with a stronger one, and shows the first 20. It inherits fgsea's high counts.
- Panel c tests nothing: no colour gradient does not rule out a link along a smaller axis.
- The script stops if a level, window or trait tested nothing, or a share exceeds 1.

## Sheets

| Sheet | Panel | One row is | Key columns |
|---|---|---|---|
| `people` | a | a person and trait | `subject`, `trait`, `value` (% change), `z` |
| `links` | b | a level, window and trait | `level`, `window`, `trait`, `n_tested`, `n_nominal`, `share` |
| `pca` | c | a biopsy and colouring trait | `sample_id`, `subject`, `timepoint`, `PC1`, `PC2`, `trait`, `value` |
| `sets` | d | a selected set and trait | `set_id`, `collection`, `pathway`, `trait`, `nes`, `fdr` |
| `profiles` | e | a profile group and screen | `group`, `n_proteins`, `main_module`, `screen`, `mean_t` |
