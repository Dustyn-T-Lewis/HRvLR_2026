# Shared functions

These files define helpers; they do not run a complete analysis. Each walkthrough loads its
packages before sourcing the files it needs. The helpers are not a standalone R package.

| File | Purpose |
|---|---|
| `analysis_reporting_helpers.R` | Repeated statistical calculations, sample alignment, phenotype tables, report heatmaps and workbook writing |
| `workbook_columns.csv` | Column definitions used by the workbook writer |
| `figure_style.R` | Colours, labels and figure styling |
| `figure_panels.R` | Reusable figure panels and figure export |
| `group_overview_figure.R` | Group-overview panels and `build_group_overview()` |
| `proteome_inputs.R` | The five Proteome clustering inputs, `validity()` of a comparison against them, `stack_inputs()` |

Key helpers:

| Function | Purpose |
|---|---|
| `remove_blood_effect()` | Remove the blood covariate effect while protecting the supplied design |
| `fit_limpa_contrasts()` | Fit limpa protein contrasts with subject blocking and robust moderation |
| `subset_for_design()` | Subset processing objects and relabel metadata for a grouping's design |
| `values_at_timepoint()` | Extract one column per participant at a selected time point |
| `within_subject_change()` | Compute change for participants measured at both time points |
| `fit_subject_model()` | Fit a model to one value per participant |
| `observed_comparison_values()` | Extract nonmissing values for the two sides of a feature comparison |
| `test_feature_discrimination()` | Compute one feature's AUC, Wilcoxon p and DeLong interval |
| `results_to_matrix()` | Arrange long results into a feature-by-comparison matrix |
| `make_report_heatmap()` | Draw a report heatmap with shared conventions |
| `paginate_report_features()` | Select and paginate features for reports |
| `write_workbook()` | Write results, column descriptions, input checksums and package versions |

`analysis_reporting_helpers.R` also loads the column dictionary when sourced. Keep meaningful
checks here; stage-specific scientific decisions belong in the stage or named shared analysis.
