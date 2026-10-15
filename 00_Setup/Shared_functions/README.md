# Shared functions

| File | Holds |
|---|---|
| `workbooks.R` | `write_workbook()`: writes a step's sheets with `read_me` and `columns` first and `input_manifest` (every file read, with its md5) and `package_versions` last |
| `workbook_columns.csv` | what every column in every workbook means; `write_workbook()` stops on a column it does not describe |

A row with a sheet pattern (`^run_map$`) applies to that sheet only and wins over a general row.
Steps source `workbooks.R` after loading dplyr and tibble.
