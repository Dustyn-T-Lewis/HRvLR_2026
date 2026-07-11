# HRvLR Proteomics Analysis

Label-free DIA-MS skeletal-muscle proteomics comparing High Responders (`HR`)
vs Low Responders (`LR`) in a 2x3 repeated-measures design:

- `T1`: baseline
- `T2`: 72 h post-training
- `T3`: 1 h acute post-bout

This repository now treats `A_YvO_2026` as the validated reference pipeline and
uses HRvLR-specific contrasts and metadata rather than copying YvO labels
forward.

## Start here

`HRvLR_pipeline.qmd` is the ground-up walkthrough and the map into everything else.
Each stage and figure also ships its own tutorial (`00_inputs.qmd`,
`01_filtering.qmd`, `02_normalization.qmd`, `03_dep.qmd`, `HRvLR_F01`–`F06.qmd`),
written to be read on its own: what it does, how, why that method and not another, and
how to read the output — including what a null looks like.

## What the pipeline found

The phenotype is real: HR and LR trained the same and grew apart (F01). **The proteome
does not separate them.** Every global test in F02 is null (PERMANOVA p = 0.62; RRPP
p ≥ 0.38; CAP fails to classify). HR and LR responses are only weakly concordant
(ρ = 0.20 training, 0.15 acute) and every fry rotation test is null. No WGCNA module
tracks the phenotype (best BH q = 0.24), and nothing predicts it out of sample.

The one genuinely FDR-controlled signal is pathway-level: 691 significant pathway ×
contrast tests, with HR's acute response enriching 217 pathways against LR's 33 (F03).

**Read the protein-level counts with care.** Most of them come from the π gate
(`p^|logFC| < 0.05`), which controls no error rate and admits proteins with raw p up to
0.269. There are 489 π-calls against 12 proteins at BH < 0.05. Treat π counts as a
ranking, not as discoveries.

Open issues that are documented but not fixed:
`docs/2026-07-11-documentation-pass-and-open-issues.md`.

## Design and Canonical Contrasts

DEP fits the means model `~ 0 + group` (one mean per `Group_Time` cell) with
`duplicateCorrelation` blocking on `Subject_ID`, and computes all 9 contrasts
below. Each is a linear combination of the six cell means, so all are estimable.

HR (within-responder):

- `Training_HR = HR_T2 - HR_T1`
- `Acute_HR = HR_T3 - HR_T2`

LR (within-responder):

- `Training_LR = LR_T2 - LR_T1`
- `Acute_LR = LR_T3 - LR_T2`

HRvLR (between-responder):

- `Baseline_HRvLR = HR_T1 - LR_T1`
- `Trained_HRvLR = HR_T2 - LR_T2`
- `Acute_HRvLR = HR_T3 - LR_T3`

Interaction (differential response):

- `Training_Interaction = (HR_T2 - HR_T1) - (LR_T2 - LR_T1)`
- `Acute_Interaction = (HR_T3 - HR_T2) - (LR_T3 - LR_T2)`

Figure defaults: volcano rings show all 9, grouped by the four families above.
Other figures (proteome overview onward) default to the 7-contrast set that
drops `Trained_HRvLR` and `Acute_HRvLR`.

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | Raw intensity matrix, metadata, phenotype table, HPA annotations |
| `01` | `01_Filtering/` | HPA presence filter, blood-concentration-gated myonuclei-rescue contaminant removal, UniProt deduplication, group-wise missingness filter, consensus outlier detection -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization of the filtered matrix; `imputation/` holds the four exploratory arms (`imp4p`, MsCoreUtils hybrid, `missForest`, Perseus MNAR), each writing a method-tagged `DAList_imputed_<method>.rds` |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, 9 HRvLR contrasts, Pi-score summaries. `b_imputed/`: exploratory DEP on the imputed matrices with logFC concordance |
| `04` | `04_Figures/` | Six figures: F01 phenotype atlas; F02 proteome overview + QC; F03 enrichVolcano ring-volcanoes, which also builds the shared fgsea cache; F04/F05 HR-vs-LR training/acute concordance, which read that cache; F06 WGCNA module-phenotype linkage |

## Canonical Run Order

### Core Stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R
```

Clustering is computed self-contained inside `04_Figures/F06` (see Figures);
WGCNA is the inferential engine for the module-phenotype linkage.

The primary DEP runs on the non-imputed normalized matrix. Imputation is
exploratory and feeds only QC and figure/WGCNA inputs. Each arm is independent
and writes a method-tagged `DAList_imputed_<method>.rds`; the `missForest` arm is
the one downstream figures and the imputed-DEP concordance check read by default:

```sh
Rscript 02_Normalization/imputation/a_script/c_missforest.R    # default downstream arm
Rscript 02_Normalization/imputation/a_script/a_imp4p.R         # exploratory alternative
Rscript 02_Normalization/imputation/a_script/b_mscoreutils.R   # exploratory alternative
Rscript 02_Normalization/imputation/a_script/d_perseus.R       # exploratory alternative (MNAR)

Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R         # exploratory imputed DEP, all four arms
```

### Figures

Six active figures plus one supplement, rendered into `b_reports`; rerun only
after stages `01` to `03` complete cleanly. F03 recomputes the fgsea enrichment
fresh each run and writes it for F04/F05 to read, so run F03 before F04/F05. The
`S_imputation` supplement reads the `03_DEP/b_imputed` concordance outputs, so run
it after the imputed DEP.

- `04_Figures/F01`: phenotype atlas
- `04_Figures/F02`: global proteome overview and QC
- `04_Figures/F03`: enrichVolcano ring-volcanoes (and the shared fgsea source data)
- `04_Figures/F04`: HR-vs-LR training-phase concordance
- `04_Figures/F05`: HR-vs-LR acute-phase concordance
- `04_Figures/F06`: WGCNA module-phenotype linkage (self-contained on the missForest-imputed proteome)
- `04_Figures/S_imputation`: imputation-method comparison supplement (non-imputed reference vs the four arms)

```sh
Rscript 04_Figures/F01/a_script/HRvLR_F01_run.R
Rscript 04_Figures/F02/a_script/HRvLR_F02_run.R
Rscript 04_Figures/F03/a_script/HRvLR_F03_run.R
Rscript 04_Figures/F04/a_script/HRvLR_F04_run.R
Rscript 04_Figures/F05/a_script/HRvLR_F05_run.R
Rscript 04_Figures/F06/a_script/HRvLR_F06_run.R
Rscript 04_Figures/S_imputation/a_script/HRvLR_S_imputation.R
```

## Repository Conventions

- `a_script/`: scripts and optional narrative notebooks
- `b_reports/`: generated PDFs and figure renders
- `c_data/`: stage outputs used by downstream steps
- `functions/`: reusable helpers that scripts source — cross-figure helpers in `04_Figures/functions/` and figure-specific helpers in each figure's `a_script/functions/` (e.g. F06's clustering loader)

## Reproducibility Rules

- Path resolution uses `here::here()` from the project root
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed normalized matrix
- repeated-measures blocking uses authoritative `Subject_ID` metadata
- acute contrasts always mean `T3 - T2`
