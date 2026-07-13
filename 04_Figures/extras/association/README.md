# phenotype mapping (association)

Which proteins and pathways associate with the continuous training responses
(mCSA, strength, fibre CSA). Association only; prediction is out of scope here
because n=16 cannot support a generalizable classifier. Rationale and citations:
`docs/phenotype_mapping_methods_review.md`.

## Method

- Collapse the normalized proteome to one value per subject per phase (baseline
  level, or the training/acute delta), then regress each feature on each trait
  with limma empirical-Bayes moderation. One observation per subject makes it
  leakage-free; moderation stabilises variance at small n.
- Pathways use singscore (`multiScore`), scored per sample so the score cannot
  leak across samples, then the same limma association on the score matrix.

## Layout

- `a_script/functions/associate.R` — `phase_subject_matrix`, `associate_traits`,
  `score_pathways`.
- `a_script/functions/plots.R` — association volcano, pathway scatter, and
  singscore's own `plotRankDensity` / `plotDispersion`.
- `a_script/run.R` — runs both levels, writes c_data tables and lead figures.
- Tests: `tests/testthat/test-phenotype-association.R` (includes the singscore
  sample-independence check).

## Result

At BH<0.05, no protein and no pathway associates with any trait at any phase.
This null is consistent with the sample size and the exercise-omics literature;
report it as such rather than reaching for prediction.
