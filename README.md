# HRvLR Proteomics: trimmed scope, limpa from precursors

DIA mass-spectrometry proteomics of skeletal muscle from a resistance-training study. Sixteen
people were biopsied before training (T1), 96 h after the last of 16 sessions (T2), and 24 h after
one further session (T3). The study ranked 25 trained people on a composite hypertrophy score,
(2 × %ΔmCSA + %ΔfCSA) / 3, and profiled the top eight as higher responders (HR) and the bottom
eight as lower responders (LR).

This branch (`dustyn/trimmed-scope-limpa`) rebuilds the analysis for HR against LR only, quantified
with limpa from DIA-NN's precursor report. The eight-grouping version is on `main`. Each phase
returns here as it is rebuilt; every folder present runs.

## Tree

```
HRvLR_V2/
├── README.md, run_all.R, renv.lock
├── 00_Setup/
│   ├── Input/              report.parquet, metadata.csv, phenotype.csv
│   └── Shared_functions/   workbooks.R, workbook_columns.csv
└── 01_Processing/
    ├── 01_Precursors/      read the report, one injection per biopsy, identity filters
    └── 02_Quantification/  detection curve, protein values with standard errors, blood index,
                            contaminants, detection filter
```

Next, in order: gene sets, set scores and modules (in `01_Processing`); the HR against LR
comparisons; the continuous trait analysis; figures (F1 design, PCA, proteins and fry sets; F2
pathways; the continuous overview).

## Every step

```
NN_Step/
├── README.md     what the step does, its inputs, how to run it, its outputs
├── a_script/     the walkthrough (.qmd): code with the reasoning beside it
├── b_reports/    the rendered walkthrough (.html): every table and preview plot
├── c_data/       what later steps read (.rds) and the step's workbook (.xlsx)
└── references/   fixed files the step reads, where it has any
```

Every workbook opens with `read_me` (one line per sheet) and `columns` (every column defined, from
`00_Setup/Shared_functions/workbook_columns.csv`) and closes with `input_manifest` (every file
read, with its md5) and `package_versions`. A step reads only files earlier steps saved, so any
step reruns alone.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R                  # every step in order
Rscript run_all.R 01_Processing    # stages whose path starts with this prefix
```

Run from the project root, so Quarto's R sessions start with the project's renv library. The runner
stops at the first failure; logs go to the untracked `Scratch/runlogs/`.

## Terms

| Term | Meaning here |
|---|---|
| precursor | one peptide ion as DIA-NN measures it; limpa quantifies proteins from their precursors |
| MBR | match-between-runs: DIA-NN searches every run against one library built from all runs |
| detection probability curve | limpa's logistic model of how likely a precursor is detected given its abundance |
| standard error | limpa's uncertainty for each protein value in each biopsy; wide where few precursors were observed |
| blood index | the biopsy's mean log2 HBB and HBA1, centred; a covariate in the comparisons |

## Limits

- A needle biopsy holds a varying share of fibre, connective tissue and fat; no covariate adjusts
  for it.
- Three biopsies (S28-T3, S28-T2, S10-T2) detected few proteins, so most of their values come from
  the detection curve. They are kept and weighted down in the comparisons; their influence is
  checked there.
