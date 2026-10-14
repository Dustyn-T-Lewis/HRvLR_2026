# 01_Protein_Filtering

Decides which proteins and runs are fit to analyse: contaminants by identity, outlier runs by a
three-of-four consensus, then proteins detected too rarely.

## Inputs

| File | Holds |
|---|---|
| `00_Setup/Input/HRvLR_raw.xlsx` | DIA-NN protein matrix, one column per run |
| `00_Setup/Input/metadata.csv` | run to person, group and time point |
| `references/blood_contaminants.csv` | 95 curated blood proteins, each with a reason |
| `references/HPA_annotations_full.tsv` | Human Protein Atlas annotation |
| `references/RBC_proteome_reference.tsv` | five published red-cell proteomes, reported only |

`references/README.md` gives the sources of the three reference tables.

## Run

```sh
quarto render 01_Processing/01_Protein_Filtering/a_script/protein_filtering.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/protein_filtering.html` | the walkthrough, with every check and figure |
| `c_data/DAList_filtered.rds` | filtered proteins and runs; read by `02_Normalisation` |
| `c_data/protein_filtering.xlsx` | verdict per protein, outlier flags per run, blood index, filter log; read by `05_Figures/Processing/S01_QC` |
