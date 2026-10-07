# R

Functions the step scripts `source()`. Nothing here runs on its own.

| File | Provides | Called by |
|---|---|---|
| `helpers.R` | `fit_limpa()`, `subject_fit()` (limma on one value per person, with intervals), `spearman_by_row()` (rho with a Fisher-z interval), `storey_q()`, `for_design()` (a processing file as one grouping sees it), phenotype tables, heatmap pages, `write_workbook()` | every step |
| `column_dictionary.csv` | what every workbook column holds; `write_workbook()` stops on a column it does not describe | every step, through `helpers.R` |
| `design_contrasts.R` | `design_contrasts()`: the design, the nine contrasts and the groups | `01`/`02_Design_Contrasts` of each grouping |
| `protein_contrasts.R` | `protein_contrasts()`: limma on limpa, fry concordance, PCA, PERMANOVA | `Protein_Contrasts` of each grouping |
| `set_contrasts.R` | `set_contrasts()`: fgsea and fry on every gene set | `Set_Contrasts` of each grouping |
| `module_contrasts.R` | `module_contrasts()`: limma on eigengenes, fry on members | `Module_Contrasts` of each grouping |
| `classify.R` | `classify_step()` and its three levels: AUC with DeLong and bootstrap intervals | `06_Classification` of each grouping |
| `group_trait_associations.R` | `group_trait_associations()`: Spearman within group, limma group x trait | `07_Trait_Associations` of each grouping |
| `associate.R` | `associate_step()`: sample model, change score, two-predictor and interaction fits | `04_Continuous_no_group/01_Trait_Associations` |
| `figure_theme.R` | colours, the contrast palette, strip boxes, `trait_label()` | every figure script |
| `figure_panels.R` | panel builders and `save_figure()` | every figure script |
| `figure_overview.R` | `figure_overview()` and `panel_clustering()`: Figure 1 for any two-group folder | the F1 scripts of `02_HRvLR_groups` and `03_Cluster_groups` |

`03_Cluster_groups/Proteome/R/inputs.R` holds the proteome folder's own helpers: the five
clustering inputs, the validity tags and `stack_inputs()`.
