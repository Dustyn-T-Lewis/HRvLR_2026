# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance-training study. Sixteen
people were biopsied before training (T1), 96 h after the last of 16 sessions (T2), and 24 h after
one further standardised session (T3). Of 48 runs, 45 pass quality control, and 1,900 of 2,400
proteins are analysed.

Twenty-five people trained. The study ranked them on a composite hypertrophy score,
(2 × %ΔmCSA + %ΔfCSA) / 3, and profiled the top eight as higher responders (HR) and the bottom
eight as lower responders (LR); the middle nine were not profiled. The repository asks its
questions three ways: with those groups, with groups found by clustering, and with no groups at
all.

## Folders

| Folder | Holds |
|---|---|
| [`00_Setup`](00_Setup/README.md) | study inputs, shared helpers and explicitly named shared analyses |
| [`01_Processing`](01_Processing/README.md) | the shared steps: filtering, normalisation, quantification, gene sets, set scores, modules |
| [`03_Group_comparisons/HRvLR`](03_Group_comparisons/HRvLR/README.md) | HR against LR: contrasts, gene sets, modules, classification, trait links, figures |
| [`03_Group_comparisons`](03_Group_comparisons/README.md) | the same steps on groups found by clustering size, strength or the proteome |
| [`04_Continuous_no_group`](04_Continuous_no_group/README.md) | everyone pooled, each trait a number |

The three analysis folders read `01_Processing` and never each other, except the parked HR/LR
Figure 3, which also reads `04_Continuous_no_group`.

## Every step folder

```
NN_Step_Name/
├── README.md     what it reads and writes, how to run it, the key result
├── a_script/     one R script, named for the step
├── b_reports/    the step's PDF report
└── c_data/       the step's workbook, plus any .rds a later step reads
```

Every workbook opens with `read_me` (one line per sheet) and `columns` (every column defined, from
`00_Setup/Shared_functions/workbook_columns.csv`), and closes with `input_manifest` (every file read, with its md5) and
`package_versions`. A step reads only files earlier steps saved, so any step reruns alone.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R            # every analysis step in order
Rscript run_all.R 01         # processing only
```

Run commands from the repository root. Each stage runs in a separate R session; the runner
stops at the first failure. Logs go to the ignored `Scratch/runlogs/` directory. Selected
folders require their upstream outputs to exist already.

Stage READMEs are maintained manually. Reports and workbooks contain the analysis results;
the runner does not rewrite documentation or run development tests.

Study inputs are in `00_Setup/Input`. Shared helpers are in `00_Setup/Shared_functions`;
reused analysis implementations are in `00_Setup/Shared_analyses`. Stage entry scripts load
these explicitly and select their own inputs and output locations. Each processing stage keeps
its reference files in its own `references/` directory; no step downloads them. Random steps
are seeded.

## Terms

| Term | Meaning here |
|---|---|
| contrast | one comparison of group means, e.g. `Training_HR` (HR's T2 minus T1) or `Training_Interaction` (HR's training change minus LR's) |
| logFC, `_lower`, `_upper` | the contrast on the log2 scale (+1 is twice as much) and its 95% interval; every estimate's interval is named after it (`slope_lower`, `rho_upper`) |
| `p`, `fdr`, `qvalue` | nominal p; Benjamini-Hochberg within the screen the sheet names; Storey's q over the same screen |
| limpa | gives every protein a value and a standard error in every biopsy from a detection-probability curve, so missing values are modelled rather than filled |
| fgsea, fry | two gene-set tests: fgsea asks whether a set ranks high against other proteins; fry asks whether its proteins move consistently across people |
| module, eigengene, kME | a WGCNA group of co-varying proteins, its one-number summary per biopsy, and a protein's correlation with that summary |
| singscore | one score per gene set per biopsy, from the set's rank within that biopsy |
| AUC | the chance a random person from the first group has a higher value than one from the second |
| floor | the baseline contrast: the groups come from an outcome, so a T1 difference shows how far the split separates people before training |

## Conventions

- Blood: an index (mean log2 HBB and HBA1 against the sample median) is a covariate in every
  limma model. Tests that take no covariate read a copy with each protein's blood slope removed
  against time point only.
- People: each person gives up to three biopsies, so person is a random block in biopsy-level
  models; change scores carry one value per person.
- Reporting: nothing is filtered on a threshold before it is written. Reports draw at most 24
  features per screen; the workbooks keep every row.
- Identifiers: proteins by UniProt accession, gene sets by `collection|name`, modules by colour.
  Module colours change if `01_Processing` is rebuilt.

## Limits

- A needle biopsy holds a varying share of fibre, connective tissue and fat; no covariate adjusts
  for it.

## Dependencies

`renv.lock` pins every package the code loads. `enrichVolcano` comes from GitHub, pinned to commit
`d3fed631f4`. Stage READMEs describe the methods and their limitations.
