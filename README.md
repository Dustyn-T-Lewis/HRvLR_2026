# HRvLR Proteomics: trimmed scope, limpa from precursors

DIA mass-spectrometry proteomics of skeletal muscle from a resistance-training study. Sixteen
people were biopsied before training (T1), 96 h after the last of 16 sessions (T2), and 24 h after
one further session (T3). The study ranked 25 trained people on a composite hypertrophy score,
(2 × %ΔmCSA + %ΔfCSA) / 3, and profiled the top eight as higher responders (HR) and the bottom
eight as lower responders (LR).

This branch (`dustyn/trimmed-scope-limpa`) rebuilds the analysis for HR against LR only, quantified
with limpa from DIA-NN's precursor report. The eight-grouping version is on `main`.

## An exploratory analysis

With 8 people per group, this is exploratory data analysis. `01_Processing` builds the values;
`02_HRvLR` and `03_Continuous` then examine the data at every level available (proteins, gene
sets, modules, variance, traits), each set against what chance gives. Those results informed the
figure plan in `04_Figures`: which panels became main figures, which went to supplements, and the
story the figures tell.

That story: both groups respond to training and to the extra bout, and HR's response reaches
significance in more gene sets than LR's (158 against 36 with training, 189 against 32 after the
bout, FDR < 0.05). The direct test of a larger HR response, the interaction, is weaker: 23 sets
for training and 58 for the bout at FDR < 0.05, and no single protein at that threshold.

## Tree

```
HRvLR_V2/
├── README.md, run_all.R, renv.lock
├── 00_Setup/
│   ├── Input/                     report.parquet, metadata.csv, phenotype.csv
│   └── Shared_functions/          workbooks.R, workbook_columns.csv
├── 01_Processing/                 from precursors to protein values, gene sets and modules
│   ├── 01_Precursors/             read the report, one injection per biopsy, identity filters
│   ├── 02_Quantification/         detection curve, protein values with standard errors, blood
│   │                              index, contaminants, detection filter
│   ├── 03_Blood_Correction/       blood slope removed, for methods without a covariate
│   ├── 04_Gene_Sets/              one protein per gene, MSigDB and GO Slim sets
│   ├── 05_Set_Scores/             singscore per set and biopsy
│   ├── 06_Module_Building/        WGCNA modules, eigengenes, kME
│   └── 07_Module_Characterisation/  gene-set ORA, top members, STRING test per module
├── 02_HRvLR/                      HR against LR; helpers.R holds the per-person and trait helpers
│   ├── 01_Protein_Contrasts/      limpa dpcDE, nine contrasts, F-tests, PCA, PERMANOVA
│   ├── 02_Set_Contrasts/          fgsea and fry per gene set
│   ├── 03_Module_Contrasts/       limma on eigengenes, fry on modules
│   ├── 04_Variance/               variance shares, DiffVar, gap statistic, PCs against covariates
│   ├── 05_Module_Preservation/    modules built in each group, preserved in the other
│   ├── 06_Classification/         single-feature AUC
│   └── 07_Group_Trait_Slopes/     trait slopes within each group, HR slope minus LR slope
├── 03_Continuous/                 traits as numbers, all 16 people pooled
│   ├── 01_Trait_Associations/     sample model, change score, two-predictor, interaction
│   ├── 02_Trait_Pathways/         fgsea on each trait ranking
│   └── 03_Profile_Clusters/       proteins grouped by their trait links
└── 04_Figures/                    main figures F1 to F5, supplements S1 to S5
    ├── style/, assets/            shared theme, palette, grids; the design schematic
    └── F1/ … F5/, S1/ … S5/       panel scripts, panel PDFs, composite PDF, one workbook
```

## Every step

```
NN_Step/
├── README.md     what the step does, its inputs, how to run it, its outputs
├── a_script/     the walkthrough (.qmd): code with the reasoning beside it
├── b_reports/    the rendered walkthrough (.html): every table and preview plot
├── c_data/       what later steps read (.rds) and the step's workbook (.xlsx)
└── references/   fixed files the step reads, where it has any
```

A figure folder differs: `a_script/` holds one R script per panel and a composite script,
`b_reports/` the panel PDFs and the composite PDF, `c_data/` one workbook with a sheet per panel
(see `04_Figures/README.md`).

Every workbook opens with `read_me` (one line per sheet) and `columns` (every column defined, from
`00_Setup/Shared_functions/workbook_columns.csv`) and closes with `input_manifest` (every file
read, with its md5) and `package_versions`. A step reads only files earlier steps saved, so any
step reruns alone.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R                  # every step in order, then the figures
Rscript run_all.R 01_Processing    # stages whose path starts with this prefix
Rscript run_all.R 02_HRvLR/03 04_Figures/F5   # several prefixes; earlier outputs must exist
```

Run from the project root, so Quarto's R sessions start with the project's renv library. A step's
`.qmd` is rendered into its `b_reports/`; a figure's composite `.R` runs with `Rscript`. The runner
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
