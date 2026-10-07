# R

Code shared by every stage. Stage scripts `source()` what they need from here; nothing in this
folder runs on its own.

| File | Provides | Sourced by |
|---|---|---|
| `helpers.R` | limpa and limma fits, `subject_fit()` (limma on one value per person, with intervals), `spearman_by_row()` (rho with a Fisher-z interval), Storey's q, phenotype tables, heatmap pages, `paginate()`, `write_workbook()` | every stage |
| `column_dictionary.csv` | what every workbook column holds; read by `write_workbook()` | every stage, through `helpers.R` |
| `classify.R` | `classify_step()`: AUC, DeLong and bootstrap intervals for the five between-group tasks | the `06_Classification` stages |
| `associate.R` | `associate_step()`: sample model, change score, two-predictor and interaction fits | `04_Continuous_no_group/01_Trait_Associations` |
| `figure_theme.R` | colours, the contrast palette, strip boxes, `trait_label()` | every figure script |
| `figure_panels.R` | panel builders (PCA, barcode, module and set panels) and `save_figure()` | every figure script |
| `figure_overview.R` | `figure_overview()` and `panel_clustering()`: Figure 1 for any two-group arc | the F1 scripts of `02_HRvLR_groups` and `03_Cluster_groups` |

`03_Cluster_groups/Proteome/R/inputs.R` holds the proteome arc's own helpers.
