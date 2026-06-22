# HRvLR Proteomics Analysis

Label-free DIA-MS skeletal-muscle proteomics comparing High Responders (`HR`)
vs Low Responders (`LR`) in a 2x3 repeated-measures design:

- `T1`: baseline
- `T2`: 72 h post-training
- `T3`: 1 h acute post-bout

This repository now treats `A_YvO_2026` as the validated reference pipeline and
uses HRvLR-specific contrasts and metadata rather than copying YvO labels
forward.

## Design and Canonical Contrasts

- `Training_HR = HR_T2 - HR_T1`
- `Training_LR = LR_T2 - LR_T1`
- `Baseline_HRvLR = HR_T1 - LR_T1`
- `Training_Interaction = (HR_T2 - HR_T1) - (LR_T2 - LR_T1)`
- `Acute_HR = HR_T3 - HR_T2`
- `Acute_LR = LR_T3 - LR_T2`
- `Acute_Interaction = (HR_T3 - HR_T2) - (LR_T3 - LR_T2)`

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | Raw intensity matrix, metadata, phenotype table, HPA annotations |
| `01` | `01_normalization/` | YvO-style HPA filter, blood contaminant removal, UniProt deduplication, group-wise missingness filter, consensus outlier detection, `cycloess` normalization |
| `02` | `02_Imputation/` | Exploratory imputation arms (`imp4p`, MsCoreUtils hybrid, `missForest`), each writing a method-tagged `DAList_imputed_<method>.rds`; primary DEP remains non-imputed |
| `03` | `03_DEP/` | `limma + duplicateCorrelation`, 7 HRvLR contrasts, Pi-score summaries, robustness analyses |
| `04` | `04_Figures/` | Figure scripts and outputs; direct YvO analogues plus HRvLR-specific extensions |
| `05` | `05_WGCNA/` | Module-level repeated-measures network analysis feeding downstream figure assets |

## Canonical Run Order

### Core Stages

```sh
Rscript 01_normalization/a_script/01_run_normalization.R
Rscript 01_normalization/a_script/02_norm_reports.R

Rscript 03_DEP/a_script/01_run_dep.R
Rscript 03_DEP/a_script/02_dep_reports.R
Rscript 03_DEP/a_script/03_dep_robustness.R
```

The primary DEP runs on the non-imputed normalized matrix. Imputation is
exploratory and feeds only QC, sensitivity, and figure/WGCNA inputs. Each arm is
independent and writes a method-tagged `DAList_imputed_<method>.rds`; the
`missForest` arm is the one downstream figures and the DEP sensitivity check read
by default:

```sh
Rscript 02_Imputation/a_script/c_missforest.R    # default downstream arm
Rscript 02_Imputation/a_script/a_imp4p.R         # exploratory alternative
Rscript 02_Imputation/a_script/b_mscoreutils.R   # exploratory alternative
```

### Figure and Network Stages

Direct-transfer and validated analogue figure scripts should be rerun only after
stages `01` to `03` complete cleanly.

- `04_Figures/F01`: phenotype
- `04_Figures/F02`: proteome overview and QC
- `04_Figures/F03`: contrast-level summary and enrichment overview
- `04_Figures/F04`: volcano rings
- `04_Figures/F05`: between-group concordance
- `05_WGCNA/a_script/01_run_wgcna.R`
- `04_Figures/F06`: within-group training vs acute concordance

`F07`, `F08`, and later exploratory figure directories exist in the repository
but should be treated as HRvLR-specific extensions, not direct YvO manuscript
transfers, unless explicitly included in the active figure plan.

## Repository Conventions

- `a_script/`: scripts and optional narrative notebooks
- `b_reports/`: generated PDFs and figure renders
- `c_data/`: stage outputs used by downstream steps

## Reproducibility Rules

- Working directory resolution uses `rprojroot::find_rstudio_root_file()`
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed normalized matrix
- repeated-measures blocking uses authoritative `Subject_ID` metadata
- acute contrasts always mean `T3 - T2`

## Transfer Provenance

The stage-by-stage transfer decision log is maintained in:

- `docs/yvo_to_hrvlr_method_mapping.md`

That document records direct transfers, partial transfers, non-transfers,
validation checks, output paths, and rerun order.
