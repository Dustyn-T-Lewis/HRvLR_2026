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
`F01_phenotype.qmd`, `F02_proteome.qmd`, `F03_pathway.qmd`, `WGCNA.qmd`),
written to be read on its own: what it does, how, why that method and not another, and
how to read the output — including what a null looks like.

## What the pipeline found

The phenotype is real: HR and LR trained the same and grew apart (F01). **The proteome
does not separate them.** Every global test in F02 is null (PERMANOVA p = 0.62; RRPP
p ≥ 0.38; CAP fails to classify). HR and LR responses are only weakly concordant
(ρ = 0.20 training, 0.15 acute) and every fry rotation test is null.

The one genuinely FDR-controlled signal is pathway-level: 691 significant pathway ×
contrast tests, with HR's acute response enriching 217 pathways against LR's 33 (F03).

**F04-F06 report a screen, not a discovery set.** Each of their 1,365 cells
(F04 = 420, F05 = 153, F06 = 792) reports a metric, a raw permutation p, the
screen size, and a leakage label — no BH q anywhere in these three, unlike
F01-F03's limma BH FDR. A screen this size and this correlated (shared subjects,
overlapping feature spaces, nested timepoint configs) has no defensible multiple-
comparison correction; a raw p plus a screen-size denominator is the honest
alternative to a q that would imply independence the cells do not have. A cell
only counts as a **lead** when `perm_p < .05` *and* its metric beats the trivial
baseline (`q2 > 0` for prediction, `auc > 0.5` for classification) — collapsed
permutation nulls otherwise manufacture significance on their own.

F05 classification is a complete null: 0 leads of 153. F06 prediction clears
32 leads of 792 (4.0%, under the 5% chance rate), yet 26 of the 32 concentrate
on `d_mcsa`. The module arm does not predict: all 8 module leads fall below
zero, median drop 0.461. That drop decomposes as 0.428 from restricting the
feature space to the two modules reproducible in every fold, and 0.000 from
rebuilding the network in-fold. The leads therefore rested on the nine modules
that do not survive subsampling, not on circular network construction, which
contributes nothing measurable here. HR and LR share module architecture
(preservation strong in both directions, though the test uses cohort-defined
labels and so is not free of the same circularity), so the WGCNA modules
describe this cohort without
distinguishing the two arms. Contrast-specific networks (built on the training
or acute contrast alone) were tested and are not viable at this cohort size.

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
by the missingness filter that the model then cannot test; the module-prediction circularity
in `shared/WGCNA` that in-fold refitting exposed (see "What the pipeline found" above); and
the transductive eigengenes in `shared/WGCNA` and `04_Figures/F03_pathway/supp`.

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
| `04` | `04_Figures/` | The results layer in arc order: F01 phenotype atlas; F02 proteome overview + QC; F03_pathway enrichVolcano ring-volcanoes, which also builds the shared fgsea source data, with HR-vs-LR training/acute concordance as its `supp`; F04_association per-cell association screen; F05_classification HR/LR classification screen; F06_prediction continuous-adaptation prediction screen; `shared/WGCNA` builds the module eigengenes F04-F06 consume |

## Canonical Run Order

### Core Stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R
```

Clustering is computed self-contained inside `04_Figures/shared/WGCNA` (see Figures);
WGCNA builds the module eigengenes that feed the `modules` level of F04, F05 and F06.

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

The results-layer figures render into `b_reports`; rerun only after stages `01`
to `03` complete cleanly. F03_pathway writes `F03_pathway_source_data.xlsx`,
which its two `supp` concordance leaves read, so run F03_pathway before the
`supp` leaves.

- `04_Figures/F01_phenotype`: phenotype atlas
- `04_Figures/F02_proteome`: global proteome overview and QC
- `04_Figures/F03_pathway`: enrichVolcano ring-volcanoes (and the shared fgsea source data)
- `04_Figures/F03_pathway/supp/concordance_training`: HR-vs-LR training-phase concordance
- `04_Figures/F03_pathway/supp/concordance_acute`: HR-vs-LR acute-phase concordance
- `04_Figures/F03_pathway/supp/summary`: magnitude and concordance in one frame
- `04_Figures/shared/WGCNA`: builds the module eigengenes (missForest-imputed
  proteome) that F04, F05 and F06 read at the `modules` level; `loso_refit/` tests
  whether the modules survive leave-one-subject-out re-definition; `preservation/`
  tests whether HR and LR share module architecture; `contrast_networks/` tests
  whether a training- or acute-only network is viable
- `04_Figures/shared/reference`: 85 worked design references (one per stage x
  level x config, plus per-level heatmaps and raw-observation detail views)
- `04_Figures/F04_association`: per-cell association screen — 420 cells over
  `<level>/<config>/<phenotype>/<method>`, three levels (proteins, pathways,
  modules) x seven configs (T1, T2, T3, training, acute, total, trajectory) x
  the six adaptation deltas plus `group_diff`
- `04_Figures/F05_classification`: HR/LR classification screen — 153 cells over
  `<level>/<config>/HR_LR/<model>`, nested leave-one-subject-out against a
  permutation null
- `04_Figures/F06_prediction`: continuous-adaptation prediction screen — 792
  cells over `<level>/<config>/<phenotype>/<model>`, nested leave-one-subject-out
  against a permutation null

Each screen runs `run_*` (compute every leaf cell), then `split_*` (fan the
pre-split `<level>/<config>/<method>` leaf into per-phenotype panels), then
`rollup_*` (all three — pool the leaves into one workbook, write MANIFEST.xlsx, and render the
specification-curve figure), then `composite_*` (assemble the figure and write
`MANIFEST.xlsx`). `functions/sweep_grid.R`'s `leaf_done()` checks the pre-split
three-level path (`<level>/<config>/<method>/c_data/results.xlsx`) because the
runners write that shape and `split_*` converts it afterward; since the
three-level leaves are deleted once `split_*` has verified them, a killed
`run_*` cannot resume from where it left off and a re-run recomputes the whole
screen from scratch.

```sh
Rscript 04_Figures/F01_phenotype/a_script/01_run_phenotype.R
Rscript 04_Figures/F02_proteome/a_script/01_run_proteome.R
Rscript 04_Figures/F03_pathway/a_script/01_run_volcanoes.R

Rscript 04_Figures/F03_pathway/supp/concordance_training/a_script/01_run_concordance_training.R
Rscript 04_Figures/F03_pathway/supp/concordance_acute/a_script/01_run_concordance_acute.R
Rscript 04_Figures/F03_pathway/supp/summary/a_script/01_run_summary.R

Rscript 04_Figures/shared/WGCNA/a_script/01_run_modules.R
Rscript 04_Figures/shared/WGCNA/preservation/a_script/01_run_preservation.R
Rscript 04_Figures/shared/WGCNA/preservation/a_script/02_run_preservation_balanced.R
Rscript 04_Figures/shared/WGCNA/contrast_networks/a_script/01_run_contrast_stability.R

# F04 association: run every cell, split into per-phenotype panels, manifest,
# composite
Rscript 04_Figures/F04_association/a_script/run_F04_association.R
Rscript 04_Figures/F04_association/a_script/split_F04_association.R
Rscript 04_Figures/F04_association/a_script/rollup_F04_association.R
Rscript 04_Figures/F04_association/a_script/composite_F04_association.R

# F05 classification: run, split, roll up (spec curve + manifest), composite
Rscript 04_Figures/F05_classification/a_script/run_F05_classification.R
Rscript 04_Figures/F05_classification/a_script/split_F05_classification.R
Rscript 04_Figures/F05_classification/a_script/rollup_F05_classification.R
Rscript 04_Figures/F05_classification/a_script/composite_F05_classification.R

# F06 prediction: run, split, roll up (spec curve + manifest), composite
Rscript 04_Figures/F06_prediction/a_script/run_F06_prediction.R
Rscript 04_Figures/F06_prediction/a_script/split_F06_prediction.R
Rscript 04_Figures/F06_prediction/a_script/rollup_F06_prediction.R
Rscript 04_Figures/F06_prediction/a_script/composite_F06_prediction.R

# Module validation last: the refit reads the F05 and F06 manifests, so it must
# follow both.
Rscript 04_Figures/shared/WGCNA/loso_refit/a_script/01_run_loso_refit.R
Rscript 04_Figures/shared/WGCNA/preservation/a_script/01_run_preservation.R
Rscript 04_Figures/shared/WGCNA/preservation/a_script/02_run_preservation_balanced.R
Rscript 04_Figures/shared/WGCNA/contrast_networks/a_script/01_run_contrast_stability.R
```

## Repository Conventions

Every stage and figure unit is `a_script/` (code), `b_reports/` (renders), `c_data/`
(outputs the next step reads).

Shared helpers live by scope:

| Path | Contents |
| --- | --- |
| `functions/` | `shared_*` helpers used across stages and figures — `shared_style.R` (palettes, theme, sizing), `shared_pca.R` (sourced by stages 01–02), `shared_utils.R`, `shared_pathway_utils.R` (fgsea/ORA); `sweep_*` helpers run the F04-F06 screen — `sweep_grid.R` (leaf paths, `leaf_done()`), `sweep_assoc.R`/`sweep_assoc_leaf.R` (F04), `sweep_pred_leaf.R` (F05/F06), `sweep_split.R`, `sweep_rollup.R`, `sweep_manifest.R`, `sweep_cell_panel.R`, `sweep_composites.R`, `sweep_speccurve.R`, `sweep_drivers.R`. |
| `04_Figures/functions/` | `f0N_*` helpers scoped to one figure — `f00_concordance.R` (the F03_pathway/supp driver) and `f00_concordance_panels.R` (its panel builders). |
| `04_Figures/shared/` | `references.bib` — the single bibliography every notebook cites; `WGCNA/` — the module source for F04-F06; `reference/` — the 85 worked design references. |
| `tests/` | The `testthat` suite. Run with `testthat::test_dir(here("tests", "testthat"))`. |

## Figures

Each figure is an `a_script/ b_reports/ c_data/` unit with its own run script. Most
ship a narrative `.qmd`; F03_pathway/supp/summary does not. F04, F05 and F06 are
organized `<level>/<config>/<phenotype-or-HR_LR>/<method>`, with `run_*` computing
every cell, `split_*`/`rollup_*` pooling them, and `composite_*` assembling the
figure and writing `MANIFEST.xlsx`.

| Directory | Question | Engine |
| --- | --- | --- |
| `F03_pathway/supp/` | Where do HR and LR adapt alike over training and the acute bout, and where do they part? How large is the response and how concordant? | Quadrant ORA, `limma::fry`, pathway NES concordance, RRHO2, bootstrap CI; median/p90 \|logFC\| and Spearman ρ. |
| `F01_phenotype/` | The phenotype: matched training, divergent growth and strength. | Phenotype atlas + linear mixed models. |
| `F02_proteome/` | Global proteome overview and QC. | PCA, DEP counts, effect sizes, set overlaps, η². |
| `F03_pathway/` | Per-contrast enrichment. | enrichVolcano ring-volcanoes, fgsea, EnrichmentMap dedup. |
| `shared/WGCNA/` | Which WGCNA modules track the phenotype, and do they generalize? | Signed WGCNA on the missForest-imputed proteome; `loso_refit/` refits the network with each subject held out; `preservation/` cross-preserves HR- and LR-only networks; `contrast_networks/` builds training- and acute-only networks. |
| `F04_association/` | Where do features associate with HR/LR or with continuous adaptation, per level and config? | `limma`/`lm`/Spearman/Wilcoxon per `<level>/<config>/<phenotype>/<method>` cell; 420 cells, raw permutation p, no BH across the screen. |
| `F05_classification/` | Can the proteome classify HR vs LR out of sample? | Elastic net, lasso, ridge, sparse PLS-DA, PAM, RF, SVM (`glmnet`, `mixOmics`, `pamr`, `randomForest`, `e1071`) per `<level>/<config>/HR_LR/<model>` cell; 153 cells, nested LOSO against a permutation null. 0 leads. |
| `F06_prediction/` | Can the proteome predict continuous adaptation out of sample? | Elastic net, lasso, ridge, sPLS, RF, SVM per `<level>/<config>/<phenotype>/<model>` cell; 792 cells, nested LOSO against a permutation null. 32 leads (4.0%), 26 of them on `d_mcsa`; all 8 module leads fall below zero once restricted to the two reproducible modules; in-fold refitting adds nothing. |

A cell in F04-F06 reports a metric, a raw permutation p, the screen size, and a
leakage label; there is no BH q anywhere in these three, because a screen this
correlated (shared subjects, overlapping feature spaces, nested configs) has no
defensible multiple-comparison correction, unlike F01-F03's per-contrast limma
BH-FDR. A **lead** requires both `perm_p < .05` and the metric beating the
trivial baseline (`q2 > 0`, `auc > 0.5`) — a collapsed permutation null can
otherwise reach its p floor while predicting worse than the group mean.
Composite hypertrophy stays out of any model carrying the HR/LR term (the
groups were defined from it); fold-specific transforms stay train-only, with
singscore's single-sample scoring the one leakage-free exception; and the
module leads are additionally flagged because the eigengenes are built once on
the full cohort, so held-out subjects still shaped the network that scores
them. At n = 16 the null is the finding.

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

F03_pathway depends on `enrichVolcano`, which is developed separately in
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

Re-run F03_pathway afterwards and diff the renders. Push the commit first: a sha that
exists only on one machine pins nothing.
