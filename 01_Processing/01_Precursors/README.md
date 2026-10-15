# 01_Precursors

Reads DIA-NN's precursor report with limpa, keeps one injection per biopsy, and applies limpa's two
identity filters.

## Inputs

| File | Holds |
|---|---|
| `00_Setup/Input/report.parquet` | DIA-NN 2 precursor report: 57 runs from 18 people |
| `00_Setup/Input/metadata.csv` | the 48 biopsies analysed: sample, person, group, time point |

## Run

```sh
quarto render 01_Processing/01_Precursors/a_script/precursors.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/precursors.rds` | limma `EList`: log2 precursor intensities (`NA` where undetected), precursor annotation, biopsy table |
| `c_data/run_map.csv` | every report run: biopsy, injection, precursors detected, kept or not |
| `c_data/precursors.xlsx` | the run map, MBR and no-MBR q-filter counts, the filter log, detection per run |
| `b_reports/precursors.html` | the rendered walkthrough |

## Read with care

- The q filter uses the match-between-runs columns (`Lib.Q.Value`, `Lib.PG.Q.Value`). The report's
  library q-values differ from its global ones, which fits an MBR search; the search log was not
  checked.
- Each biopsy keeps its injection with the most precursors. Re-injections replace the first
  injection for S17-T1 and S29-T1.
- S28-T2 and S28-T3 are thin (6,699 and 4,048 precursors against a median of 13,102) and are kept.
