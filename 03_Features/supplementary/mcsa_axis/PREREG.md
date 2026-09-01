# mCSA axis — pre-registration

Written 2026-08-19, before any protein model was fitted. These rules do not
change after results are seen. This stage closes the one empty cell in the
phenotype grid: whether whole-muscle CSA change has a protein-level correlate
of its own. It is not a search for a way to overturn the null. A hit that
fails the gates below is a bug in this stage, not biology.

Seed 42 before every stochastic step.

## What is already settled, and is not re-tested here

Verified against `00_input/c_data/phenotype.csv` and
`04_Figures/supp_screens/c_data/pred_cells_summary.csv` on 2026-08-19.

`d_mcsa` is an ingredient of the composite the arms were defined on, not a
separate phenotype. Regressing `comp_hypertrophy` on the z-scores of the three
CSA measures returns R² = 0.964 with weights 4.98 (fCSA I), 2.41 (fCSA II),
4.12 (mCSA); dropping `d_mcsa` takes R² to 0.737. `d_mcsa` alone separates the
arms at AUC 0.83, Wilcoxon p = 0.028. What is true instead: `d_mcsa` carries a
third of the composite's formula and almost none of its ordering. The galamm
measurement model weights the six indicators from the data and returns mCSA at
0.27 (z = 0.9, 8% of variance) against fCSA at 1.03 and 1.02, and its factor
scores correlate 0.94 with the composite and 0.32 with `d_mcsa`.

LR_S14 carries much of that orthogonality: the study's largest whole-muscle
gain (z = +1.63) with near-worst fibre change (z = −1.75, −1.80). Removing it
moves `d_mcsa` vs fCSA I from 0.09 to 0.32 and vs the composite from 0.49 to
0.69. This stage treats that as a live measurement-error hypothesis and
reports every protein result twice, with and without the subject.

The HR/LR split has no middle to trim (LR tops at 0.85, HR starts at 8.06).
No analysis here subsets to extreme responders. No analysis here pools the six
traits into one "response".

Already answered elsewhere and not repeated: proteins against the continuous
latent response axis (`04_galamm_pilot`, min BH q = 0.117, zero survivors);
all six outcomes at module, pathway and delta-cluster level
(`05_phenotype_modules`, nothing beats its consistency null); the nine arm ×
timepoint contrasts (`F04_association`, zero BH survivors at protein and
module level).

## Fixed inputs

- Protein subset: the complete-case rows of
  `02_Normalization/c_data/DAList_normalized.rds`, observed in all 45
  analysed samples, nothing imputed. `pilot_data()` in
  `03_Features/04_galamm_pilot/a_script/pilot_helpers.R` applies the gate and
  asserts n = 931. Reused, not rebuilt.
- Blood index: `functions/blood_index_model.R::blood_index_data()`, joined by
  `pilot_data()`, no NA after the join. It enters every model as a fixed
  covariate; the T3 confound does not cancel (arm × T3 b = −1.21, p = 0.032).
- Phenotype: `00_input/c_data/phenotype.csv`, 16 subjects.
- Sample table: `functions/shared_hlm.R::hlm_meta()`. Not rebuilt.

S28 is T1-only and S29 is T2/T3-only, so a single-timepoint config holds 15
subjects and the training delta holds 14. n = 16 is unreachable for anything
per-subject and is not claimed.

## The model

One sample per subject within a config, so a mixed model buys nothing.
Per config: `limma::lmFit` on the 931 × n matrix with
`model.matrix(~ d_mcsa + blood_index)`, `eBayes`, coefficient 2, BH across the
931 within that config.

Four configs, each labelled by what it can claim:

| Config | Samples | n | Claim |
| --- | --- | --- | --- |
| `T1` | baseline | 15 | the only baseline forecast |
| `T2` | 72 h post-training | 15 | concurrent |
| `T3` | 1 h acute | 15 | concurrent, carries the blood confound |
| `delta` | T2 − T1 per subject | 14 | concurrent |

`d_mcsa` is the change in area over the period ending at T2, so only `T1`
precedes its outcome. A survivor appearing in `T2`, `T3` or `delta` and not in
`T1` is a concurrent correlate, and the figure says so on its face.

## Permutation

`d_mcsa` is shuffled across subjects, B = 200, seed 42, refit, recount. Subject
is the unit of randomisation because the phenotype is a subject property, the
scheme of `03_Features/01_Proteins/a_script/pi_permutation.R`. Empirical p
carries the +1 correction, `(n_ge + 1) / (B + 1)`.

The permutation arms only behind a BH survivor, as `03_q2_protein.R` does. With
zero survivors it never runs and the stage reports that plainly.

## Decision rules

1. A **survivor** is BH q < 0.05 within its config across the 931.
2. A survivor must also clear empirical p < 0.05 against the B = 200
   permutation of its own config.
3. Zero survivors is written as a plain negative. No trend language, no
   lowered threshold, no fallback to a π gate — π controls no error rate and
   this project has already shown shuffling beats it (permuted median 274
   against 235 observed).
4. Every config is run twice, with LR_S14 and without. The second is reported
   as a parallel column, never as a gate, and neither run is preferred after
   the fact.
5. The blood index enters every model. No exceptions.

## Expected outcome

A negative. The module level put `d_mcsa`'s best cell at consistency p = 0.057
on a module that fails leave-one-subject-out rebuilding, and the pathway level
put its `d_mcsa` hits inside their permuted distribution. The scripts are
written so that a negative reads as the result rather than as a failure to
find one.

## Deliverables

`c_data/01_phenotype_geometry.csv` and `01_composite_weights.csv`;
`c_data/02_protein_mcsa.csv` (per protein × config × subject-set) and
`02_permutation.csv` if it arms; `b_reports/F_mcsa_axis.{png,pdf}` with
`F_mcsa_axis_legend.md`; `tests/testthat/test-mcsa-axis.R`; `FINDINGS.md`
carrying a "Deviations from PREREG.md" section and a proposed README
paragraph left uninserted.
