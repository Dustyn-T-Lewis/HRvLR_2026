# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance-training study. Sixteen
people were biopsied three times: at baseline (T1), after the training programme (T2), and after
one acute bout at the end of it (T3). Of 48 runs, 45 enter the analysis;
`01_Processing/01_Protein_Filtering` records which three are dropped and why.

The study labelled each person a High or Low Responder (HR, LR) by a median split of a composite
hypertrophy score. That split is mostly a fibre-area split, so the repository asks its questions
three ways:

| Arc | Who is compared | Question |
|---|---|---|
| [`02_HRvLR_groups`](02_HRvLR_groups/README.md) | the study's HR and LR | do the proteome and its response differ between responder groups? |
| [`03_Cluster_groups`](03_Cluster_groups/README.md) | groups found by clustering people on size, on strength, or on their proteome | do groups defined another way differ? |
| [`04_Continuous_no_group`](04_Continuous_no_group/README.md) | everyone pooled, each trait a number | which proteins, sets and modules track size and strength? |

All three read the same [`01_Processing`](01_Processing/README.md) and never each other. 01_Processing uses no
group label in any value it produces, with one exception: the detection filter in
`01_Protein_Filtering` keeps a protein when it is detected often enough within at least one of the
six HR/LR-by-time cells.

## Map

```
HRvLR_V2/
├── 00_Input/                        study data, no code
├── 01_Processing/                        shared steps every arc reads
│   ├── 01_Protein_Filtering/        contaminants, outlier runs, detection filter
│   ├── 02_Normalisation/            cyclic loess
│   ├── 03_Quantification/           limpa: a value and a standard error for every protein in every biopsy
│   ├── 04_Gene_Sets/                MSigDB and GO Slim sets mapped to the measured proteins
│   ├── 05_Set_Scores/               singscore: one score per set per biopsy
│   ├── 06_Module_Building/          WGCNA co-expression modules and eigengenes
│   └── 07_Module_Characterisation/  what each module is: enriched sets, hubs, STRING
├── 02_HRvLR_groups/
│   ├── 01_Design_Contrasts/         six cell means, nine contrasts
│   ├── 02_Protein_Contrasts/        limma on every protein
│   ├── 03_Set_Contrasts/            fgsea and fry on every gene set
│   ├── 04_Module_Contrasts/         limma on eigengenes, fry on module members
│   ├── 05_Module_Preservation/      are each group's modules found in the other?
│   ├── 06_Classification/           can one protein, set or module tell HR from LR?
│   │   ├── 01_Proteins/  02_Sets/  03_Modules/
│   ├── 07_Trait_Associations/       does a feature's link to a trait differ by group?
│   └── Figures/                     F1 to F3 and supplements S01 to S03
├── 03_Cluster_groups/
│   ├── Muscle_size/                 01_Clustering, then copies of the HR/LR stages on Hi and Lo
│   ├── Muscle_strength/             the same, clustered on strength
│   └── Proteome/                    clustered five ways on the proteome; 08_Cluster_Phenotype
├── 04_Continuous_no_group/
│   ├── 01_Trait_Associations/       limma, Spearman and bicor of every feature against each trait
│   ├── 02_Trait_Pathways/           fgsea on each trait's protein ranking
│   ├── 03_Profile_Clusters/         proteins grouped by their pattern of trait links
│   └── Figures/                     F1
├── R/                               shared code (below)
├── docs/                            decisions log and concept images
└── run_all.R                        runs every stage in order
```

## What every stage folder holds

```
NN_Stage_Name/
├── README.md        what the stage asks, how it works, how to read its outputs
├── a_script/        one R script, named for the stage (stage_name.R)
├── b_reports/       stage_name_report.pdf: overview pages, then one panel per notable feature
└── c_data/          stage_name.xlsx, plus any .rds a later stage reads
```

Every workbook has the same frame:

| Sheet | Holds |
|---|---|
| `read_me` | one line per sheet |
| `columns` | every column of every sheet, with what it holds (from `R/column_dictionary.csv`) |
| result sheets | named for what one row is, e.g. `protein_by_contrast` |
| `input_manifest` | every file the script read, with its md5 |
| `package_versions` | packages loaded at run time |

Columns use one vocabulary across the repository. An estimate is named for what it is (`logFC`,
`slope`, `slope_difference`, `interaction`, `rho`, `auc`) and its 95% interval is named after it
(`logFC_lower`, `logFC_upper`, `rho_lower`, ...). `p` is the nominal p-value, `fdr` the
Benjamini-Hochberg adjustment within the screen the `read_me` names, `qvalue` Storey's q over the
same screen.

## Reading order

1. `01_Processing/README.md`, then its stages in order: how a raw run becomes a protein value.
2. `02_HRvLR_groups/README.md` and its `Figures/F1_Overview`: the study's own question.
3. `02_HRvLR_groups/02_Protein_Contrasts/README.md`: the contrasts every later stage reuses.
4. `04_Continuous_no_group`: the same data without the split.
5. `03_Cluster_groups`: the split made other ways.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R            # every stage in order, figures last; logs in .runlogs/
Rscript run_all.R 01 02      # only the arcs whose folder starts with 01 or 02
```

`00_Input/README.md` has the download command for the STRING files
`01_Processing/07_Module_Characterisation` reads.

## Shared code in `R/`

| File | Holds |
|---|---|
| `helpers.R` | limpa and limma fits, `subject_fit()` (limma on one value per person, with intervals), `spearman_by_row()` (rho with a Fisher-z interval), Storey's q, phenotype tables, heatmap pages, `write_workbook()` |
| `column_dictionary.csv` | what every workbook column holds; a column it does not describe stops the run |
| `classify.R` | the classification step shared by the protein, set and module levels |
| `associate.R` | the continuous association step shared by the three levels |
| `figure_theme.R` | colours, the contrast palette, strip boxes, label helpers |
| `figure_panels.R`, `figure_overview.R` | panel builders; `figure_overview()` builds Figure 1 for any two-group arc |

## Conventions every stage shares

- Missing values: limpa (Li, Cobbold & Smyth 2025) fits a detection-probability curve, the chance
  a peptide is seen given its abundance, and gives every protein a value and a standard error in
  every biopsy. Protein contrasts use those standard errors as weights; every other test reads the
  complete matrix, and protein rows carry `observed_share`, the share of values measured rather
  than estimated.
- Blood: an index (mean log2 HBB and HBA1 against the sample median) enters every limma model as a
  covariate. Tests that take no covariate read a copy with each protein's blood slope removed
  (`limma::removeBatchEffect`, design = time point, so no group label enters).
- People: each person gives up to three biopsies, so person is a random block in biopsy-level
  models, and change scores carry one value per person.
- Reporting: each package's own statistic, nominal p and BH FDR, with Storey's q beside them.
  Nothing is filtered on a threshold before it is written; reports draw at most 24 features per
  screen and the workbooks keep every row.
- Identifiers: proteins by UniProt accession, gene sets by `collection|name`, modules by colour; every
  result table carries them as `feature` with a readable `label`. Module colours come from WGCNA and
  change if 01_Processing is rebuilt.

## Limits stated rather than tested

- Biopsy composition: a needle biopsy holds a varying share of fibre, connective tissue and fat.
  Large co-moving sarcomere changes could partly reflect that share; no covariate adjusts for it.
- The composite hypertrophy formula and the biopsy timing (programme length, time from bout to
  T3) are not in the repository.

## Dependencies

`renv.lock` pins every package. One comes from outside CRAN and Bioconductor's current releases:
`enrichVolcano`, pinned to GitHub commit `d3fed631f4` for the HR/LR volcano supplement.

The first version of this project, with its prediction screen, `functions/`, `tests/` and
`archive/`, is the git tag `legacy-v1`. Decisions and their reasons are logged in
`docs/decisions.md`.
