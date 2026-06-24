# YvO to HRvLR Method Mapping

Updated: 2026-05-02

## Scope

This document maps the validated `A_YvO_2026` proteomics pipeline onto
`A_HRvLR_2026`.

The goal is not to clone YvO code. The goal is to preserve validated analysis
logic where the HRvLR design supports it, re-derive all HRvLR factor levels and
contrasts from HRvLR metadata, and explicitly reject YvO stages whose biology
does not transfer.

## Target Design Summary

### YvO Reference Design

- Biological question:
  - how age modifies baseline proteome state and chronic resistance-training
    response
- Design:
  - 2 groups (`Young`, `Old`)
  - 2 timepoints (`Pre`, `Post`)
  - repeated measures within subject
- Canonical contrasts:
  - `Training_Young = Young_Post - Young_Pre`
  - `Training_Old = Old_Post - Old_Pre`
  - `Aging = Old_Pre - Young_Pre`
  - `Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)`

### HRvLR Target Design

- Biological question:
  - how responder status (`HR`, `LR`) differs at baseline, after chronic
    training, and after an acute post-training bout
- Design:
  - 2 groups (`HR`, `LR`)
  - 3 timepoints (`T1`, `T2`, `T3`)
  - repeated measures within subject
  - 16 subjects total, 8 per responder group, 3 samples per subject
- Metadata columns confirmed:
  - `Col_ID`
  - `Subject_ID`
  - `Group`
  - `Timepoint`
  - `Group_Time`
- Canonical HRvLR contrasts supported by the design:
  - `Training_HR = HR_T2 - HR_T1`
  - `Training_LR = LR_T2 - LR_T1`
  - `Baseline_HRvLR = HR_T1 - LR_T1`
  - `Trained_HRvLR = HR_T2 - LR_T2`
  - `Acute_HRvLR = HR_T3 - LR_T3`
  - `Training_Interaction = (HR_T2 - HR_T1) - (LR_T2 - LR_T1)`
  - `Acute_HR = HR_T3 - HR_T2`
  - `Acute_LR = LR_T3 - LR_T2`
  - `Acute_Interaction = (HR_T3 - HR_T2) - (LR_T3 - LR_T2)`

### Design Implications

- YvO repeated-measures logic transfers directly to HRvLR.
- YvO baseline group comparison transfers directly.
- YvO chronic training comparison transfers directly.
- HRvLR adds a valid acute-response axis that YvO does not have.
- YvO aging-specific or reversal-specific interpretation does not transfer.

## Core Pipeline Mapping

### Stage 00: Inputs, Metadata, Sample IDs

| YvO source | YvO question | YvO assumptions | HRvLR target | Transfer status | HRvLR decision |
| --- | --- | --- | --- | --- | --- |
| `00_input/YvO_raw.xlsx`, `00_input/YvO_meta.xlsx` | Are sample annotations authoritative before any filtering? | Sample columns must match metadata `Col_ID`; repeated-measures subject ID exists | `00_input/HRvLR_raw.xlsx`, `00_input/HRvLR_meta.csv`, `00_input/Phenotype.xlsx` | Direct | Keep the same authoritative metadata-first contract. Use HRvLR `Group`, `Timepoint`, `Group_Time`, `Subject_ID`; never derive factors from raw column names unless validating. |

Required HRvLR validation checks:

- `setequal(raw sample columns, metadata$Col_ID)`
- required metadata columns exist before any stage runs
- `Subject_ID` is present and non-missing for all repeated-measures analyses
- `Group_Time` levels are re-derived from HRvLR metadata, not hard-coded from YvO

Primary downstream outputs and claims supported:

- every stage
- every figure
- all contrast directionality

### Stage 01: Normalization and QC

| YvO source | YvO question | YvO assumptions | HRvLR target | Transfer status | HRvLR decision |
| --- | --- | --- | --- | --- | --- |
| `01_normalization/a_script/01_normalize.R` | How do we convert raw protein intensities into a muscle-specific, contamination-filtered, normalized matrix suitable for downstream statistics? | HPA skeletal-muscle filter, blood contaminant removal, UniProt deduplication, group-wise missingness filter, consensus outlier removal, `cycloess` normalization | `01_normalization/a_script/01_run_normalization.R`, `02_norm_reports.R` | Direct with threshold adaptation | Preserve the YvO stage order and QC logic. Keep `cycloess`, the same contaminant logic, deterministic ordering, and consensus outlier framework. Adapt the missingness threshold to HRvLR cell size: `>=5/8` per `Group_Time` is the HRvLR analogue of YvO `>=10` in larger cells. |

Required HRvLR metadata columns:

- `Col_ID`
- `Subject_ID`
- `Group`
- `Timepoint`
- `Group_Time`

Primary outputs and claims supported:

- `01_normalization/c_data/02_normalized.csv`
- `01_normalization/c_data/03_DAList_normalized.rds`
- pre/post normalization QC
- sample retention/exclusion
- stage 02 and stage 03 inputs

Validation checks:

- sample order in the intensity matrix exactly matches metadata order
- outlier consensus uses only retained HRvLR metadata
- no raw files overwritten

### Stage 02: Missingness Classification and Imputation

| YvO source | YvO question | YvO assumptions | HRvLR target | Transfer status | HRvLR decision |
| --- | --- | --- | --- | --- | --- |
| `02_Imputation/a_script/01_impute.R` | Which proteins are MAR vs MNAR, and can we generate a deterministic complete matrix for QC-sensitive or complete-data downstream tasks without using imputation in the primary DEP model? | 3-method MAR/MNAR consensus, within-matrix missingness profiling, `missForest`, deterministic row ordering before stochastic imputation, unreliability flag at `>50%` missing | Exploratory arms: `a_imp4p.R`, `b_mscoreutils.R`, `c_missforest.R`, each writing `DAList_imputed_<method>.rds` | Superseded by 2026 audit | Stage 02 restructured to CvH-style lean per-method arms (imp4p / MsCoreUtils / missForest). Imputation stays out of the primary DEP model; the `missForest` arm is what downstream QC/figures read. The YvO-style 3-method MAR/MNAR consensus, per-protein classification, and the 12-method benchmark were removed. |

Transferred assumptions:

- primary DEP uses non-imputed normalized data
- imputed matrix exists for sensitivity analyses and complete-data tools
- stochastic imputation must be seeded
- deterministic ordering must be enforced before `missForest`

Non-transferred assumptions:

- YvO two-timepoint group structure is not copied directly
- HRvLR-specific benchmark-first workflow is not the canonical transferred stage

Primary outputs and claims supported:

- `02_Imputation/c_data/01_imputed.csv`
- `02_Imputation/c_data/01_DAList_imputed.rds`
- `02_Imputation/c_data/02_mar_mnar_classification.csv`
- `02_Imputation/c_data/07_imputation_mask.csv`
- `02_Imputation/c_data/08_mnar_imputation_audit.csv`
- `02_Imputation/c_data/09_imputation_summary.txt`
- `02_Imputation/c_data/02_imputation.xlsx`
- missingness/imputation QC figures

Validation checks:

- zero `NA` values remain in the imputed matrix
- imputed matrix dimensions match normalized matrix dimensions
- imputation reports read the same top-level output paths written by the canonical run script

### Stage 03: Differential Abundance

| YvO source | YvO question | YvO assumptions | HRvLR target | Transfer status | HRvLR decision |
| --- | --- | --- | --- | --- | --- |
| `03_DEP/a_script/01_run_dep.R` | Which proteins change within group, differ at each biological state, or show response interactions, under repeated-measures blocking? | non-imputed normalized input, `limma`, `duplicateCorrelation`, explicit contrast table, empirical Bayes shrinkage, raw `p <= 0.10` plus `pi_score < 0.05` summary logic | `03_DEP/a_script/01_run_dep.R`, `02_dep_reports.R`, `03_dep_robustness.R` | Direct with expanded HRvLR contrast set | Keep YvO’s repeated-measures `limma + duplicateCorrelation` framework. Use HRvLR `Subject_ID` as the blocking variable. Retain YvO-style effect-size and Pi-score outputs. Expand the contrast table to the 9 HRvLR contrasts supported by the 2x3 design. |

Direct YvO analogues:

- `Aging` -> `Baseline_HRvLR`
- `Training_Young` -> `Training_HR`
- `Training_Old` -> `Training_LR`
- `Interaction` -> `Training_Interaction`

HRvLR-specific additions:

- `Trained_HRvLR`
- `Acute_HRvLR`
- `Acute_HR`
- `Acute_LR`
- `Acute_Interaction`

Primary outputs and claims supported:

- `03_DEP/c_data/01_limma_DAList.rds`
- `03_DEP/c_data/03_combined_results.csv`
- per-contrast tables and summaries
- static volcano and mean-difference plots
- downstream figure inputs

Validation checks:

- blocking must use authoritative HRvLR subject metadata, not regex-stripped sample IDs when `Subject_ID` is available
- acute contrasts must consistently mean `T3 - T2`
- contrast direction must be explicit in both tables and figure subtitles

### Stage 03 Robustness and Supplementary DEP Summaries

| YvO source | YvO question | YvO assumptions | HRvLR target | Transfer status | HRvLR decision |
| --- | --- | --- | --- | --- | --- |
| `03_DEP/a_script/03_run_robustness.R` | Are the core DEP patterns stable to resampling and alternative summaries? | same design matrix and contrast definitions as primary DEP | `03_DEP/a_script/03_dep_robustness.R` | Partial direct transfer | Keep robustness logic only where it remains interpretable for HRvLR. Leave-one-subject-out, effect-size CI summaries, and imputation sensitivity are valid. HRvLR-specific response profiling comparing `HR` vs `LR` effect-size distributions is a valid analogue, not a direct YvO transfer. |

## Figure Mapping

### Direct or Partial Figure Transfers

| YvO figure | YvO question | HRvLR target | Transfer status | HRvLR interpretation |
| --- | --- | --- | --- | --- |
| `F00` pipeline QC | Did filtering, outlier handling, MAR/MNAR classification, imputation, and model setup behave as expected? | Optional stitched `04_Figures/F00` composite plus canonical stage-local QC in `01_normalization`, `02_Imputation`, and `03_DEP` | Implemented optional composite | Keep QC generation inside each stage subdir as the canonical truth. The stitched `F00` should remain a manuscript-style composite over those stage-local outputs, not a replacement for them. |
| `F01` phenotype | Are phenotype outcomes coherent with the biological grouping? | `04_Figures/F01` | Direct with expansion | HRvLR `F01` should become a full phenotype atlas using all available proteomics-linked phenotype endpoints from `00_input/Phenotype.xlsx`, not a trimmed outcome subset. |
| `F02` proteome overview | What is the global overview of QC, variance, effect distributions, and proteome structure? | `04_Figures/F02` | Partial direct transfer | Keep the YvO overview feel, but define HRvLR `F02` explicitly as QC plus proteome overview. CV and variability summaries are valid additions and should be integrated rather than split into a separate figure. |
| `F03` DEP and pathway overview | How many DEPs and enriched pathways are present, and how do they distribute across contrasts? | `04_Figures/F03` | Partial direct transfer | Keep DEP-count, overlap, barcode, and pathway-summary logic, but aggregate across the HRvLR main contrast set rather than forcing the YvO 2x2 contrast layout. |
| `F03` volcano rings | How do per-contrast protein- and pathway-level effects look in a consistent visual system? | `04_Figures/F04` | Direct with expansion | Same figure type transfers. HRvLR volcano rings show all 9 contrasts grouped by family (HR, LR, HRvLR, Interaction). Other figures default to the 7-contrast set that drops `Trained_HRvLR` and `Acute_HRvLR`. Acute subtitles must use `T3 - T2`, not `T3 - T1`. |
| `F04` training concordance | Are between-group responses concordant? | `04_Figures/F05` | Direct with HRvLR extension | `Training_HR` vs `Training_LR` is the direct HRvLR analogue. `Acute_HR` vs `Acute_LR` is a valid HRvLR-specific main-panel extension and should remain in the main figure. |
| `F05` reversal/concordance analogue | Does HRvLR support a non-aging analogue of reversal-style concordance? | `04_Figures/F06` | Analogue only | Within-group `Training` vs `Acute` concordance is valid for HRvLR, but it is supplementary support for response-shape interpretation, not a reversal analysis. |
| `F06` WGCNA module-trait | Do network modules relate to traits or group/time effects? | `05_WGCNA` outputs and downstream figure assets | Direct with design adaptation | Repeated-measures network modeling is valid for HRvLR. The module-level interpretation must reference responder status and chronic/acute contrasts, not aging. |

### Non-Transfer or HRvLR-Specific Analogues

| YvO figure | Why it does not transfer directly | HRvLR-specific analogue |
| --- | --- | --- |
| `F05` aging reversal | HRvLR has no age axis and no biologically meaningful “reversal toward young baseline” interpretation | Within-group training-versus-acute concordance is valid for HRvLR and is currently represented by `04_Figures/F06`. It is an analogue, not a reversal analysis. |
| `F07` phenotype prediction | YvO prediction logic depends on a specific age/training framing and a validated LOSO objective | HRvLR should only keep prediction or classifier figures if a responder-status prediction question is explicitly specified and the train/test contract is validated. This is not part of the current canonical transfer. |

### Figure Style Provenance Rules

Transferred unchanged where valid:

- shared visual language
- font sizing hierarchy
- legend placement conventions
- export size discipline
- panel label style

Re-derived for HRvLR:

- all contrast labels
- all timepoint subtitles
- all biological interpretations
- all significance annotations tied to the acute axis

Current figure-specific issues still remaining after the current implementation pass:

- `04_Figures/F02` and `04_Figures/F03` now have explicit run entrypoints and
  canonical `MAIN_F0X_composite` outputs, but their panel internals still
  retain some older layout assumptions and should only be refactored
  cautiously.
- `04_Figures/F04` now follows the agreed main/supp split, but its export tree
  still carries some legacy panel files from the earlier all-main layout and
  should be cleaned cautiously rather than overwritten blindly.
- `F01` now covers the phenotype endpoint families present in
  `HRvLR_meta.csv`; if the manuscript needs full session-by-session volume-load
  trajectories from `00_input/Phenotype.xlsx`, that should be added as an
  explicit HRvLR-only expansion rather than assumed to be part of the direct
  YvO transfer.

## Non-Canonical or Exploratory HRvLR Components

These exist in the repository but should not be treated as validated YvO
transfers without separate review:

- the benchmark-first stage-02 workflow under `02_Imputation/a_script/benchmark/`
- figure directories whose biological role is not yet tied back to the YvO
  reference manuscript structure
- any responder-prediction workflow not explicitly validated against an HRvLR
  outcome definition

## Validation Checklist

Every canonical rerun should verify:

- required metadata columns exist
- intensity matrix columns match metadata `Col_ID`
- `Subject_ID` exists and is used for repeated-measures blocking where claimed
- factor levels come from HRvLR metadata
- `Training_*` contrasts use `T2 - T1`
- `Acute_*` contrasts use `T3 - T2`
- `Baseline_HRvLR` uses `HR_T1 - LR_T1`
- `Trained_HRvLR` uses `HR_T2 - LR_T2`
- `Acute_HRvLR` uses `HR_T3 - LR_T3`
- stochastic steps are seeded with `42`
- row and column ordering is deterministic before imputation or label selection
- no raw input files are overwritten
- top-level README run commands point to real scripts

## Canonical Rerun Order

### Core Pipeline

1. `01_normalization/a_script/01_run_normalization.R`
2. `01_normalization/a_script/02_norm_reports.R`
3. `03_DEP/a_script/01_run_dep.R`
4. `03_DEP/a_script/02_dep_reports.R`
5. `03_DEP/a_script/03_dep_robustness.R`

Imputation is exploratory (not part of the primary DEP path). Run the arm(s) you
need before figures/WGCNA; `c_missforest.R` is the default downstream arm:

- `02_Imputation/a_script/c_missforest.R`
- `02_Imputation/a_script/a_imp4p.R`
- `02_Imputation/a_script/b_mscoreutils.R`

### Direct-Transfer Figures

1. `04_Figures/F00/*` (optional stitched composite only after stage-local QC exists)
2. `04_Figures/F01/*`
3. `04_Figures/F02/*`
4. `04_Figures/F03/*`
5. `04_Figures/F04/*`
6. `04_Figures/F05/*`
7. `04_Figures/F06/*`

### Valid HRvLR-Specific Analogues

1. `05_WGCNA/a_script/01_run_wgcna.R`
2. `04_Figures/F06/*`
3. `04_Figures/F07/*` only after confirming the enrichment question is the
   intended manuscript figure
4. `04_Figures/F08/*` only after confirming the acute-enrichment figure is the
   intended manuscript figure

## Immediate Cleanup Priorities From This Audit

1. Decide whether `F01` should stay anchored to `HRvLR_meta.csv` only, or gain
   an HRvLR-specific expansion for full per-session volume-load trajectories
   from `00_input/Phenotype.xlsx`.
2. Clean legacy figure exports cautiously where older filenames from prior
   layouts still coexist beside the canonical refreshed outputs.
3. If needed for manuscript polish, refactor the remaining older internal panel
   layout code in `F02` and `F03` without changing the refreshed output
   contract.
