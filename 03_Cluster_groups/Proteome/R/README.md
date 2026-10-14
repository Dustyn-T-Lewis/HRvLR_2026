# R

Shared code for the proteome groupings: which clusterings exist, how far each comparison can be trusted
against them, and how their five workbooks become one.

`inputs.R` is sourced by stages `02_Design_Contrasts` to `07_Trait_Associations` and by
`Figures/F1_Overview`, after `00_Setup/Shared_functions/analysis_reporting_helpers.R`. `01_Clustering` and `08_Cluster_Phenotype` do not
source it.

| Object | What it provides |
|---|---|
| `inputs_used` | The five clustering inputs in run order: `T1`, `T2`, `T3`, `training`, `acute`. Every stage loops over it. `01_Clustering/a_script/clustering.R` keeps its own copy; change both together. |
| `time_points` | The time points behind each input (`training` is T1 and T2, `acute` T2 and T3), plus `all` for T1 to T3. |
| `uses(name)` | The time points a contrast, task or window name uses, read from its prefix: `T1`, `T2` and `T3` as named; `Baseline` or `baseline` T1; `Trained` and `post_level` T2; `Acute_C1vC2` and `post_bout_level` T3; `Training` or `training` the training window; `Acute` or `acute` the acute window; `all` and `subject mean` every time point. NA for any other name. |
| `validity(input, name)` | `circular` when the comparison uses exactly the input's time points, `caution` when it shares at least one, `clean` when it shares none, NA when `uses()` gives NA. Why each tag matters is in [`../README.md`](../README.md). |
| `stack_inputs(stage, leave_out)` | Reads the one workbook in each `<stage>/<input>/c_data/`, stacks every result sheet with an `input` column and, where a `contrast`, `task` or `window` column exists, a `validity` column, and writes `<stage>/c_data/<content>_all_inputs.xlsx` through `write_workbook()`. Sheets are read as text and typed after stacking, because an empty sheet otherwise reads as logical. `read_me`, `columns`, `input_manifest`, `package_versions`, `plot_index` and any sheet in `leave_out` are not stacked. |

`stack_inputs()` names the stacked file after the T1 workbook and takes the sheet list from it,
so every input's workbook must hold the same sheets.

```
R/
└── inputs.R
```
