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
Most stages and figures also ship their own tutorial (`00_inputs.qmd`,
`01_filtering.qmd`, `02_normalization.qmd`, `03_dep.qmd`, and the figure narratives
`F01_phenotype.qmd`, `F02_proteome.qmd`, `F03_pathway.qmd`, `F04_modules.qmd`),
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
0.270. There are 506 π-calls against 15 proteins at BH < 0.05 (21 at BH < 0.10), and 70 of
the π-calls have raw p ≥ 0.05. Treat π counts as a ranking, not as discoveries.

The null is robust to imputation. On BH, missForest (MAR), MsCoreUtils (hybrid) and Perseus
(MNAR) each return zero significant proteins in all five HR-vs-LR and interaction contrasts,
matching the non-imputed arm. It is also robust to the blood filter: readmitting all 136
blood-tagged proteins leaves every one of those contrasts at zero.

Known limitations are stated on the page where the reader meets them: the π gate in
`HRvLR_pipeline.qmd`; the human-only search space with no contaminant FASTA and no decoys, so
reagent contaminants cannot be detected at all (`01_filtering.qmd`); the 34 proteins admitted
by the missingness filter that the model then cannot test; the circular module test in F04; and
the transductive eigengenes in F04 and `04_Figures/F03_pathway/supp`.

## Design and Canonical Contrasts

DEP fits the means model `~ 0 + group` (one mean per `Group_Time` cell) with
`duplicateCorrelation` blocking on `Subject_ID`, and computes all 9 contrasts
below. Each is a linear combination of the six cell means, and is estimable only
for proteins observed in every cell it touches. 34 proteins reach the model with
at least one empty cell, so the true tested-N is 1,897–1,912 per contrast, never
1,920 (`03_DEP/a_non_imputed/b_reports/untested_proteins.csv`).

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
| `04` | `04_Figures/` | Seven figures: F03_pathway/supp HR-vs-LR training/acute concordance, which reads F03's source data; F01 phenotype atlas; F02 proteome overview + QC; F03 enrichVolcano ring-volcanoes, which also builds the shared fgsea source data; F04 WGCNA module-phenotype linkage; F05 continuous training-response association; F06 out-of-sample prediction |

## Canonical Run Order

### Core Stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R
```

Clustering is computed self-contained inside `04_Figures/F04_modules` (see Figures);
WGCNA is the inferential engine for the module-phenotype linkage.

The primary DEP runs on the non-imputed normalized matrix. Imputation is
exploratory and feeds only QC and figure/WGCNA inputs. Each arm is independent
and writes a method-tagged `DAList_imputed_<method>.rds`; the `missForest` arm is
the one downstream figures and the imputed-DEP concordance check read by default:

```sh
Rscript 02_Normalization/imputation/a_script/impute_missforest.R   # default downstream arm
Rscript 02_Normalization/imputation/a_script/impute_imp4p.R        # exploratory alternative
Rscript 02_Normalization/imputation/a_script/impute_mscoreutils.R  # exploratory alternative
Rscript 02_Normalization/imputation/a_script/impute_perseus.R      # exploratory alternative (MNAR)

Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R         # exploratory imputed DEP, all four arms
Rscript 03_DEP/b_imputed/a_script/02_imp4p_circularity.R       # permutation control: why imp4p breaks the null
Rscript 02_Normalization/imputation/a_script/imputation_qc.R   # imputation QC figure (reads the imputed DEP)
```

On BH the null holds under every imputer that does not impute inside the tested factor.
missForest, MsCoreUtils and Perseus each return **zero** BH-significant proteins in all five
HR-vs-LR and interaction contrasts, matching the non-imputed arm. imp4p returns 117, because
`impute.mle` fits a separate EM within each `Group_Time` cell and the contrasts then test among
those same cells; re-imputing within random cells of equal size collapses it back to zero. Never
rank robustness on π-counts — π exponentiates |log2FC|, the quantity imputation distorts, and
Perseus tops the π table while being null on BH.

### Figures

Seven figures, rendered into `b_reports`; rerun only after stages `01` to `03`
complete cleanly. F03 writes `F03_pathway_source_data.xlsx`, which the two
F03_pathway/supp concordance leaves read, so run F03 before F03_pathway/supp.

- `04_Figures/F01_phenotype`: phenotype atlas
- `04_Figures/F02_proteome`: global proteome overview and QC
- `04_Figures/F03_pathway`: enrichVolcano ring-volcanoes (and the shared fgsea source data)
- `04_Figures/F03_pathway/supp/concordance_training`: HR-vs-LR training-phase concordance
- `04_Figures/F03_pathway/supp/concordance_acute`: HR-vs-LR acute-phase concordance
- `04_Figures/F03_pathway/supp/summary`: magnitude and concordance in one frame
- `04_Figures/F04_modules`: WGCNA module-phenotype linkage (self-contained on the missForest-imputed proteome)
- `04_Figures/F05_association`: continuous training-response association
- `04_Figures/F06_prediction`: out-of-sample HR/LR and continuous prediction, six leaves

```sh
Rscript 04_Figures/F01_phenotype/a_script/01_run_phenotype.R
Rscript 04_Figures/F02_proteome/a_script/01_run_proteome.R
Rscript 04_Figures/F03_pathway/a_script/01_run_volcanoes.R

Rscript 04_Figures/F03_pathway/supp/concordance_training/a_script/01_run_concordance_training.R
Rscript 04_Figures/F03_pathway/supp/concordance_acute/a_script/01_run_concordance_acute.R
Rscript 04_Figures/F03_pathway/supp/summary/a_script/01_run_summary.R

Rscript 04_Figures/F04_modules/a_script/01_run_modules.R
Rscript 04_Figures/F05_association/a_script/01_run_association.R

Rscript 04_Figures/F06_prediction/prediction_responder/baseline/a_script/01_run_baseline.R
Rscript 04_Figures/F06_prediction/prediction_responder/training/a_script/01_run_training.R
Rscript 04_Figures/F06_prediction/prediction_responder/acute/a_script/01_run_acute.R
Rscript 04_Figures/F06_prediction/prediction_continuous/baseline/a_script/01_run_baseline.R
Rscript 04_Figures/F06_prediction/prediction_continuous/training/a_script/01_run_training.R
Rscript 04_Figures/F06_prediction/prediction_continuous/acute/a_script/01_run_acute.R
```

## Repository Conventions

Every stage and figure unit is `a_script/` (code), `b_reports/` (renders), `c_data/`
(outputs the next step reads).

Shared helpers live by scope:

| Path | Contents |
| --- | --- |
| `functions/` | `shared_*` helpers used across stages and figures — `shared_style.R` (palettes, theme, sizing), `shared_pca.R` (sourced by stages 01–02), `shared_utils.R`, `shared_pathway_utils.R` (fgsea/ORA). |
| `04_Figures/functions/` | `f0N_*` helpers scoped to one figure — `f00_concordance.R` (the F03_pathway/supp driver) and `f00_concordance_panels.R` (its panel builders). |
| `04_Figures/shared/` | `references.bib` — the single bibliography every notebook cites. |
| `tests/` | The `testthat` suite. Run with `testthat::test_dir(here("tests", "testthat"))`. |

## Figures

Each figure is an `a_script/ b_reports/ c_data/` unit with its own run script. Most
ship a narrative `.qmd`; F05 and F03_pathway/supp/summary do not, and F06 carries its triad
and its two `a_narrative.qmd` files inside the six prediction leaves rather than at
the top level.

| Directory | Question | Engine |
| --- | --- | --- |
| `F03_pathway/supp/` | Where do HR and LR adapt alike over training and the acute bout, and where do they part? How large is the response and how concordant? | Quadrant ORA, `limma::fry`, pathway NES concordance, RRHO2, bootstrap CI; median/p90 \|logFC\| and Spearman ρ. |
| `F01_phenotype/` | The phenotype: matched training, divergent growth and strength. | Phenotype atlas + linear mixed models. |
| `F02_proteome/` | Global proteome overview and QC. | PCA, DEP counts, effect sizes, set overlaps, η². |
| `F03_pathway/` | Per-contrast enrichment. | enrichVolcano ring-volcanoes, fgsea, EnrichmentMap dedup. |
| `F04_modules/` | Which WGCNA modules track the phenotype? | Signed WGCNA on the missForest-imputed proteome, `limma::fry`, LOSO q². |
| `F05_association/` | Which proteins and pathways associate with the continuous training responses (ΔmCSA, strength, ΔfCSA)? | Mixed models on proteins and singscore pathway scores. Association only; no protein or pathway survives BH. |
| `F06_prediction/` | Can baseline, training-response, or acute features predict HR vs LR out of sample? | Elastic net (`glmnet`) + sparse PLS-DA (`mixOmics`), nested LOSO CV against a permutation null. Both arms null after BH. |

Prediction is scored against a permutation null, never zero; composite hypertrophy
stays out of any model carrying the HR/LR term (the groups were defined from it);
and fold-specific transforms stay train-only, with singscore's single-sample
scoring the one leakage-free exception. At n = 16 the null is the finding.

## Reproducibility Rules

- Path resolution uses `here::here()` from the project root
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed normalized matrix
- repeated-measures blocking uses authoritative `Subject_ID` metadata
- acute contrasts always mean `T3 - T2`
- packages are pinned in `renv.lock` and load from `renv/library`, not the system
  library; `.Rprofile` activates this on its own. `renv::restore()` rebuilds the
  library from the lockfile, `renv::status()` reports drift

## Working on enrichVolcano Without Disturbing This Pipeline

F03 depends on `enrichVolcano`, which is developed separately in
`D_Tools/enrichVolcano`. renv keeps the two apart: this project reads only
`renv/library`, and the system library is not on its search path. Installing a
work-in-progress build of the package the usual way puts it in the *system*
library, where this pipeline cannot see it, so iterate freely — the figures keep
rendering against the pinned commit.

The pin is deliberate, so adopting a new build is deliberate too:

```r
renv::install("Dustyn-T-Lewis/enrichVolcano@<sha>")
renv::snapshot()
```

Re-run F03 afterwards and diff the renders. Push the commit first: a sha that
exists only on one machine pins nothing.
