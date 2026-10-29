# Shared functions

| File | Holds |
|---|---|
| `workbooks.R` | `write_workbook()`: writes a step's sheets with `read_me` and `columns` first and `input_manifest` (every file read, with its md5) and `package_versions` last |
| `workbook_columns.csv` | what every column in every workbook means; `write_workbook()` stops on a column it does not describe |
| `person_traits.R` | `person_values()` (one column per person: a level or a change), `trait_windows()` and `window_values()` (which feature value pairs with which trait), `load_features()` (proteins, set scores and eigengenes with their labels) |
| `set_tree.R` | `set_tree()` (overlapping gene sets grouped by Jaccard distance, average linkage, cut at three heights) and `most_specific()` (each cluster's representative) |

A row with a sheet pattern (`^run_map$`) applies to that sheet only and wins over a general row.
Steps source these after loading dplyr (and tidyr and purrr for `person_traits.R`).
