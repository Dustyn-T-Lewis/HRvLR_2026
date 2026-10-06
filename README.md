# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance training study. Sixteen
subjects were biopsied at baseline (T1), after the training programme (T2) and after an acute bout
at the end of it (T3). They were labelled High or Low Responder by a median split of a composite
hypertrophy score, eight each. 48 MS runs; three are dropped as outliers, leaving 45.

The primary comparison is the training interaction: whether the proteome changed differently over
training in HR than in LR. The study also asks which proteins, pathways and co-expression modules
separate the groups, and which track the size and strength phenotypes as continuous traits.

## Stages

| Stage | Substages | Writes for later stages |
|---|---|---|
| [`00_Input/`](00_Input/README.md) | study data, no code | |
| [`01_Preprocess/`](01_Preprocess/README.md) | `01_Filter`, `02_Normalize`, `03_Quantify` (limpa) | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`02_Proteins/`](02_Proteins/README.md) | `01_Design`, `02_Contrasts`, `03_Classify`, `04_Associate` | `design.rds`, `fit.rds`, `associate.rds` |
| [`03_Pathways/`](03_Pathways/README.md) | `01_Gene_Sets`, `02_Scores`, `03_Contrasts`, `04_Classify`, `05_Associate` | `gene_sets.rds`, `singscore.rds`, `set_tests.rds` |
| [`04_Modules/`](04_Modules/README.md) | `01_Build`, `02_Contrasts`, `03_Classify`, `04_Associate`, `05_Characterise`, `06_Preserve` | `modules.rds` |
| [`05_Integration/`](05_Integration/README.md) | `01_Chance_Line`, `02_Candidates`, `03_Panel`, `04_Continuous`, `05_Clustering` | `shuffled_t.rds` |
| [`06_Figures/`](06_Figures/README.md) | figures built from saved outputs only | |

Proteins, pathways and modules run the same three steps: `Contrasts` (the nine contrasts and a
baseline-adjusted ANCOVA), `Classify` (eight single-feature classification tasks) and `Associate`
(phenotype, in groups and as continuous traits). `05_Integration` works across the three levels.

Every substage holds `a_script/` (one script), `b_reports/` (`<step>_figures.pdf`, 11 x 8.5 in),
`c_data/` (one workbook, plus any `.rds` a later step reads) and a README naming what it reads and
writes. Every workbook opens on `read_me` and ends with `input_manifest` (inputs with md5) and
`package_versions`. Data passes through disk, so a substage reruns alone once its inputs exist.

Shared code lives in `R/`: `helpers.R` (blood correction, limpa and limma fits, ANCOVA, Spearman via
`cor.test`, Storey's q, phenotype tables and their shuffle, heatmap and workbook plumbing),
`classify.R` and `associate.R` (the `Classify` and `Associate` steps, called by each level),
`figure_theme.R` (colours and theme for every figure) and `figure_panels.R` (panel builders the
figures share).

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R            # every substage in order; logs in .runlogs/
Rscript run_all.R 05         # one stage
```

`00_Input/README.md` has the download command for the STRING files `01_Trunk/07_Characterise`
and `05_Integration` read.

## Statistical conventions

- Missing values: limpa (Li, Cobbold & Smyth 2025) fits a detection probability curve and gives
  every protein a value in every sample (`dpcQuantByRow`). The protein contrasts use `dpcDE` with
  the standard errors as precision weights. Every other test reads that complete matrix; each
  protein row carries `observed_share`, the share of its values measured rather than estimated.
- Blood: an index (mean log2 HBB and HBA1 against the sample median) enters every limma model as a
  covariate. Tests that take no covariate read a copy with each protein's blood slope removed
  (`limma::removeBatchEffect`).
- Subjects: HR and LR are different people, so subject is a random block in sample-level models;
  change scores carry one value per subject.
- Multiple testing: BH within each screen (contrast, task, window and outcome, trait and term), with
  Storey's q beside it. Nothing is gated on either; tables report every nominal result.
- Chance: shuffling whole phenotype rows across subjects keeps the traits' correlations and breaks
  only their link to the proteome. Each continuous screen reports its count and strongest hit
  against 200 shuffles and an empirical FDR (shuffle median over observed count); candidate counts,
  panels, pathway rankings and clusterings carry their own shuffle nulls.
- Identifiers: proteins by `uniprot_id`, gene sets by `set_id`, modules by colour; every result
  table carries them as `feature`, with a readable `label`.

## Dependencies

`renv.lock` pins every package. Three come from outside CRAN and Bioconductor's current releases:
`enrichVolcano` (GitHub commit `d3fed631f4`, the build that reproduces the project's set labels and
still exports `volcano_ring()`, `volcano_ring_theme()` and `ev_clean_label()`), `DGCA` (GitHub,
`andymckenzie/DGCA`, archived from CRAN) and `biclust` 2.0.3 (CRAN archive).

## Version 1 is the tag legacy-v1

The first version of this project, with its prediction screen (F06), `functions/`, `tests/` and
`archive/`, is the git tag `legacy-v1`.
