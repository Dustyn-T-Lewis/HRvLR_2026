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

fgsea reports 1,004 significant pathway × contrast tests of 14,330, concentrated in
the acute contrasts (F03). `limma::fry` over the same Hallmark and GO Slim sets, which
rotates residuals and so carries inter-gene correlation instead of assuming it away,
returns **zero** in every HR-vs-LR and interaction contrast. Where the two disagree the
rotation test is the one whose null holds.

**No protein, pathway or module survives BH in any of the nine contrasts.** F04 fits
all nine at three feature levels and every cell is empty at q < .05. The smallest q
anywhere in the study is 0.0715, in `Acute_HR`. Earlier runs had eight survivors in
`Acute_LR`; six of them (SPTB, ANK1, STOM, CAT, BLVRB, SYNE2) were red-cell proteins
and the seventh, LCP1, is a leukocyte protein. All seven are now removed as blood by
`00_input/blood_contaminants.csv`.

**The T3 blood confound does not cancel in the interaction.** T3 biopsies carry
roughly twice the blood of T1 and T2, and the rise is not equal in the two arms: on
the log2 haemoglobin index the arm × T3 term is b = −1.21, p = 0.032, and p = 0.017 by
subject-label permutation, with LR rising 2.14 against HR's 0.57. A difference of
differences removes a constant offset, not a differential one.
`03_DEP/a_non_imputed/a_script/05_blood_adjusted.R` refits every contrast with the
index as a covariate.

**F05-F06 report a screen, not a discovery set.** Each of their 945 cells
(F05 = 153, F06 = 792) reports a metric, a permutation p over B = 200, the
screen size, and a leakage label — no BH q anywhere in these two, unlike
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
0.269. There are 464 π-calls against 0 proteins at BH < 0.05 (2 at BH < 0.10), and 68 of
the π-calls have raw p ≥ 0.05.

π is no longer the selection criterion, because it does not select. Shuffling the arm
label across subjects produces **more** π hits than the real labels do: 235 observed
against a permuted median of 274 across the five HR-vs-LR and interaction contrasts, and
no contrast reaches empirical p = 0.22
(`03_DEP/a_non_imputed/a_script/04_pi_permutation.R`). Never quote a π count without that
null beside it.

The null is robust to imputation. On BH, missForest (MAR), MsCoreUtils (hybrid) and Perseus
(MNAR) each return zero significant proteins in all five HR-vs-LR and interaction contrasts,
matching the non-imputed arm. It is also robust to the blood filter: readmitting every
blood-tagged protein leaves each of those contrasts at zero.

Known limitations are stated on the page where the reader meets them: the π gate in
`HRvLR_pipeline.qmd`; the human-only search space with no contaminant FASTA and no decoys, so
reagent contaminants cannot be detected at all (`01_filtering.qmd`); the 34 proteins admitted
by the missingness filter that the model then cannot test; the module-prediction circularity
in `04_Features/modules` that in-fold refitting exposed (see "What the pipeline found" above); and
the transductive eigengenes in `04_Features/modules` and `05_Figures/F03_pathway/supp`.

## Design and Canonical Contrasts

DEP fits the means model `~ 0 + group` (one mean per `Group_Time` cell) with
`duplicateCorrelation` blocking on `Subject_ID`, and computes all 9 contrasts
below. Each is a linear combination of the six cell means, and is estimable only
for proteins observed in every cell it touches. 34 proteins reach the model with
at least one empty cell, so the true tested-N is 1,877–1,892 per contrast, never
1,900 (`03_DEP/a_non_imputed/b_reports/bh_denominators.csv`).

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
| `01` | `01_Filtering/` | curated blood + handling contaminant list, blood-concentration-gated myonuclei-rescue HPA removal, UniProt deduplication, group-wise missingness filter, consensus outlier detection -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization of the filtered matrix; `imputation/` holds the four exploratory arms (`imp4p`, MsCoreUtils hybrid, `missForest`, Perseus MNAR), each writing a method-tagged `DAList_imputed_<method>.rds` |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, 9 HRvLR contrasts, Pi-score summaries. `b_imputed/`: exploratory DEP on the imputed matrices with logFC concordance |
| `04` | `04_Features/` | The derived-feature layer. `proteins/` carries stage 03's tested protein matrix forward; `pathways/` scores singscore pathway sets; `modules/` builds the signed WGCNA modules and tests whether they generalize (`loso_refit/`, `preservation/`, `contrast_networks/`). Each arm fits the same nine contrasts stage 03 fits, so all three feature levels share one estimator, and each writes its own QC report and results table |
| `05` | `05_Figures/` | The results layer in arc order: F01 phenotype atlas; F02 proteome overview + QC; F03_pathway enrichVolcano ring-volcanoes, which also builds the shared fgsea source data, with HR-vs-LR training/acute concordance as its `supp`; F04_association the HR-vs-LR contrast heatmaps; F05_classification HR/LR classification screen; F06_prediction continuous-adaptation prediction screen; F07_keepers the best cell each model reached against its own null |

## Canonical Run Order

### Core Stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R
```

Clustering is computed self-contained inside `04_Features/modules` (see Figures);
WGCNA builds the module eigengenes that feed the `modules` level of F05 and F06.

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

- `05_Figures/F01_phenotype`: phenotype atlas
- `05_Figures/F02_proteome`: global proteome overview and QC
- `05_Figures/F03_pathway`: enrichVolcano ring-volcanoes (and the shared fgsea source data)
- `05_Figures/F03_pathway/supp/concordance_training`: HR-vs-LR training-phase concordance
- `05_Figures/F03_pathway/supp/concordance_acute`: HR-vs-LR acute-phase concordance
- `05_Figures/F03_pathway/supp/summary`: magnitude and concordance in one frame
- `04_Features/modules`: builds the module eigengenes (missForest-imputed
  proteome) that F05 and F06 read at the `modules` level; `loso_refit/` tests
  whether the modules survive leave-one-subject-out re-definition; `preservation/`
  tests whether HR and LR share module architecture; `contrast_networks/` tests
  whether a training- or acute-only network is viable
- `05_Figures/shared/reference`: 85 worked design references (one per stage x
  level x config, plus per-level heatmaps and raw-observation detail views)
- `archive/pooled_association_2026-07-29`: the retired pooled-association sweep
  (the old 420-cell F04 and the F07 synthesis heatmaps). Pooled rho partly
  re-expresses the HR-vs-LR contrast — the HR/LR label correlates +0.755 with
  `d_fcsa_I` — so the screen answered a question the study is not asking
- `05_Figures/F04_association`: the nine stage 03 contrasts at protein, pathway
  and module level, in three column families (HR-LR by timepoint, the two
  interactions, the four within-arm changes). One panel per level plus a stacked
  composite; logFC fill on a scale per family, BH within each contrast
- `05_Figures/F05_classification`: HR/LR classification screen — 153 cells over
  `<level>/<config>/HR_LR/<model>`, nested leave-one-subject-out against a
  permutation null
- `05_Figures/F06_prediction`: continuous-adaptation prediction screen — 792
  cells over `<level>/<config>/<phenotype>/<model>`, nested leave-one-subject-out
  against a permutation null

Each screen runs `run_*` (compute every leaf cell), then `split_*` (fan the
pre-split `<level>/<config>/<method>` leaf into per-phenotype panels), then
`rollup_*` (all three — pool the leaves into one workbook, write MANIFEST.xlsx, and render the
specification-curve figure), then `composite_*` (assemble the figure and write
`MANIFEST.xlsx`). `functions/sweep_grid.R`'s `leaf_done()` checks the pre-split
three-level path (`<level>/<config>/<method>/c_data/results.xlsx`) because the
runners write that shape and `split_*` converts it afterward.

Each runner writes a leaf as soon as it finishes it, so a killed `run_*`
**resumes** — relaunch it and `leaf_done()` skips every completed leaf, costing
at most the one cell that was in flight. `split_*` copies the pre-split tree
and never deletes it, so resume keeps working until that tree is removed by
hand. Verified 2026-07-28: a killed F05 run resumed at leaf 119 of 153.

Both runners take `[levels] [max B]`. **Pass `"" 200`.** `bmax` defaults to
`1000`, which expands the grid to `c(0, 200, 1000)` and computes a
1,000-permutation null instead of 200 — five times the work for p resolution
(1/1001 vs 1/201) that nothing in this screen uses. Launching without arguments
has cost this project two long detours; see `docs/decisions.md` 2026-07-26.

```sh
Rscript 05_Figures/F01_phenotype/a_script/01_run_phenotype.R
Rscript 05_Figures/F02_proteome/a_script/01_run_proteome.R
Rscript 05_Figures/F03_pathway/a_script/01_run_volcanoes.R

Rscript 05_Figures/F03_pathway/supp/concordance_training/a_script/01_run_concordance_training.R
Rscript 05_Figures/F03_pathway/supp/concordance_acute/a_script/01_run_concordance_acute.R
Rscript 05_Figures/F03_pathway/supp/summary/a_script/01_run_summary.R

Rscript 04_Features/modules/a_script/01_run_modules.R
Rscript 04_Features/modules/preservation/a_script/01_run_preservation.R
Rscript 04_Features/modules/preservation/a_script/02_run_preservation_balanced.R
Rscript 04_Features/modules/contrast_networks/a_script/01_run_contrast_stability.R

# The nine contrasts on each feature level. Proteins reshape stage 03's fit and
# assert a refit reproduces it; modules and pathways are fitted here.
Rscript 04_Features/proteins/a_script/01_run_proteins.R
Rscript 04_Features/pathways/a_script/01_run_pathways.R
Rscript 04_Features/modules/a_script/02_run_module_contrasts.R

# F04 contrast heatmaps: one per level, then the stacked all-levels sheet.
# The composite re-runs each level, so running it alone is enough.
Rscript 05_Figures/F04_association/a_script/composite_F04_association.R

# F05 classification: run, split, roll up (spec curve + manifest), composite
Rscript 05_Figures/F05_classification/a_script/run_F05_classification.R
Rscript 05_Figures/F05_classification/a_script/split_F05_classification.R
Rscript 05_Figures/F05_classification/a_script/rollup_F05_classification.R
Rscript 05_Figures/F05_classification/a_script/composite_F05_classification.R

# F06 prediction: run, split, roll up (spec curve + manifest), composite
Rscript 05_Figures/F06_prediction/a_script/run_F06_prediction.R
Rscript 05_Figures/F06_prediction/a_script/split_F06_prediction.R
Rscript 05_Figures/F06_prediction/a_script/rollup_F06_prediction.R
Rscript 05_Figures/F06_prediction/a_script/composite_F06_prediction.R

# Module validation last: the refit reads the F05 and F06 manifests, so it must
# follow both.
Rscript 04_Features/modules/loso_refit/a_script/01_run_loso_refit.R
Rscript 04_Features/modules/preservation/a_script/01_run_preservation.R
Rscript 04_Features/modules/preservation/a_script/02_run_preservation_balanced.R
Rscript 04_Features/modules/contrast_networks/a_script/01_run_contrast_stability.R
```

## Repository Conventions

Every stage and figure unit is `a_script/` (code), `b_reports/` (renders), `c_data/`
(outputs the next step reads).

Shared helpers live by scope:

| Path | Contents |
| --- | --- |
| `functions/` | `shared_*` helpers used across stages and figures — `shared_style.R` (palettes, theme, sizing), `shared_pca.R` (sourced by stages 01–02), `shared_utils.R`, `shared_pathway_utils.R` (fgsea/ORA); `f00_concordance.R` and `f00_concordance_panels.R` drive F03_pathway/supp; `sweep_*` helpers run the F05-F06 screen — `sweep_grid.R` (leaf paths, `leaf_done()`), `sweep_pred_leaf.R` (F05/F06), `sweep_split.R`, `sweep_rollup.R`, `sweep_manifest.R`, `sweep_cell_panel.R`, `sweep_composites.R`, `sweep_speccurve.R`, `sweep_drivers.R`. |
| `05_Figures/shared/` | `references.bib` — the single bibliography every notebook cites; `WGCNA/` — the module source for F05-F06; `reference/` — the worked design references. |
| `tests/` | The `testthat` suite. Run with `testthat::test_dir(here("tests", "testthat"))`. |

## Figures

Each figure is an `a_script/ b_reports/ c_data/` unit with its own run script. Most
ship a narrative `.qmd`; F03_pathway/supp/summary does not. F05 and F06 are
organized `<level>/<config>/<phenotype-or-HR_LR>/<method>`, with `run_*` computing
every cell, `split_*`/`rollup_*` pooling them and writing `MANIFEST.xlsx`, and
`composite_*` assembling the figure.

| Directory | Question | Engine |
| --- | --- | --- |
| `F03_pathway/supp/` | Where do HR and LR adapt alike over training and the acute bout, and where do they part? How large is the response and how concordant? | Quadrant ORA, `limma::fry`, pathway NES concordance, RRHO2, bootstrap CI; median/p90 \|logFC\| and Spearman ρ. |
| `F01_phenotype/` | The phenotype: matched training, divergent growth and strength. | Phenotype atlas + linear mixed models. |
| `F02_proteome/` | Global proteome overview and QC. | PCA, DEP counts, effect sizes, set overlaps, η². |
| `F03_pathway/` | Per-contrast enrichment. | enrichVolcano ring-volcanoes, fgsea, EnrichmentMap dedup. |
| `04_Features/modules/` | Which WGCNA modules track the phenotype, and do they generalize? | Signed WGCNA on the missForest-imputed proteome; `loso_refit/` refits the network with each subject held out; `preservation/` cross-preserves HR- and LR-only networks; `contrast_networks/` builds training- and acute-only networks. |
| `F04_association/` | How do high responders differ from low responders, per feature level? | The nine stage 03 contrasts read from `04_Features`; logFC fill with one scale per contrast family, stars for nominal p, black box for BH q < .05 within a contrast. Zero survivors in all nine contrasts at all three levels; smallest q anywhere is 0.0715. |
| `F05_classification/` | Can the proteome classify HR vs LR out of sample? | Elastic net, lasso, ridge, sparse PLS-DA, PAM, RF, SVM (`glmnet`, `mixOmics`, `pamr`, `randomForest`, `e1071`) per `<level>/<config>/HR_LR/<model>` cell; 153 cells, nested LOSO against a permutation null. 0 leads. |
| `F06_prediction/` | Can the proteome predict continuous adaptation out of sample? | Elastic net, lasso, ridge, sPLS, RF, SVM per `<level>/<config>/<phenotype>/<model>` cell; 792 cells, nested LOSO against a permutation null. 32 leads (4.0%), 26 of them on `d_mcsa`; all 8 module leads fall below zero once restricted to the two reproducible modules; in-fold refitting adds nothing. |

A cell in F05-F06 reports a metric, a permutation p, the screen size, and a
leakage label; there is no BH q anywhere in these two, because a screen this
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
