# galamm pilot — pre-registration

Written 2026-08-18, before any model was fitted. These rules do not change
after results are seen. The pilot tests whether three model-specification
choices are holding signal down; it is not a search for a way to overturn
the null. A new hit that fails the checks below is a bug in the pilot, not
biology.

Engine: galamm 0.4.0 (Sørensen et al. 2023, Psychometrika 88:456-486,
doi:10.1007/s11336-023-09910-z), installed via `renv::install("galamm")`
and recorded in `renv.lock`. Q1 uses `nlme::lme` first because it answers
the same question at lower cost. Seed 42 before every stochastic step.

## Fixed inputs

- Protein matrix: `02_Normalization/c_data/DAList_normalized.rds`,
  1,900 proteins x 45 analysed samples (the three partial-roster biopsies
  never entered the pipeline and do not enter here).
- Protein subset, all questions except Q3: the complete-case set — rows
  observed in all 45 samples, nothing imputed. Verified n = 931 before
  writing this document; every script asserts the count and stops if it
  drifts. This is the same set F06's survivor rescore used.
- Sample table: `functions/shared_hlm.R::hlm_meta()`. Not rebuilt.
- Blood index: `functions/blood_index_model.R::blood_index_data()`, joined
  on `Col_ID`; assert no NA after the join. Not rebuilt. The index enters
  every model in this pilot as a fixed covariate — the T3 confound does
  not cancel in the interaction (arm x T3 b = -1.21, p = 0.032).
- Phenotype: `00_input/c_data/phenotype.csv`, 16 subjects, six indicators
  (`comp_hypertrophy`, `d_fcsa_I`, `d_fcsa_II`, `d_mcsa`,
  `d_1rm_legpress`, `d_1rm_ext`; the last has one NA, which drops that
  single indicator row in long format and nothing else).

## Q1 — per-timepoint residual variance

Model, per protein: `y ~ group * timepoint + blood_index`, random
intercept per subject, fitted twice with `nlme::lme` (REML):
once homoscedastic, once with `weights = varIdent(~ 1 | timepoint)`.
Same fixed effects in both, so the REML likelihood-ratio test on 2 df is
valid for the variance structure.

Readouts, decided in advance:

1. Per-protein variance ratios sigma_T2/sigma_T1 and sigma_T3/sigma_T1;
   report medians and IQRs over the 931.
2. The six contrasts of `shared_hlm.R` (`group_main`, `T1`, `T2`, `T3`,
   `training`, `acute`) via `emmeans` on the lme fits (emmeans supports
   lme; the galamm caveat about hand-built contrasts applies only if
   galamm is reached). BH within contrast across the 931, q < 0.05.
3. Side-by-side with the existing lmerTest fit on the same proteins,
   stating that lme uses containment df and lmerTest Satterthwaite, so
   part of any p-value difference is the df approximation, not the
   variance model.

Close rule: Q1 is negative if both median ratios lie inside [0.8, 1.25]
and no contrast gains a BH q < 0.05 protein that the homoscedastic fit
lacked. Escalate to galamm only if lme cannot express something needed
(none anticipated).

## Q2 — latent hypertrophic-response factor

Order is fixed: measurement model first, on phenotypes alone; proteins
only if it passes; the d_mcsa comparison before any protein is read.

1. Standardize all six indicators to z-scores across the 16 subjects
   (units differ by three orders of magnitude). Anchor:
   `comp_hypertrophy` loading fixed to 1 — it is the a-priori composite;
   after standardization the anchor sets scale only.
2. Fit the one-factor measurement model in galamm on the six indicators
   in long format (95 rows: 16 x 6 minus the one NA). Report every
   loading with its SE and the share of each indicator's variance the
   factor explains.
   Stop rule: if the model does not converge, is singular, or every
   free loading has |lambda|/SE < 2, Q2 closes as "not identifiable at
   n = 16" and no protein model is fitted. That is the honest answer.
3. Extract empirical-Bayes factor scores; report Spearman rho against
   `d_mcsa`. If rho > 0.9, the factor is d_mcsa, the pilot has explained
   the F06 concentration rather than extended it, and Q2 closes there —
   protein fits are then optional confirmation, not discovery.
4. Only if the measurement model passes and rho <= 0.9: per-protein joint
   model in galamm — the protein's abundance (45 samples) enters as an
   additional Gaussian response loading on the latent factor, with
   `group * timepoint + blood_index` fixed effects and `(1 | subject)`
   on the protein response. Test statistic: the Wald z of the protein's
   loading. BH across the 931, q < 0.05. galamm inference is Wald from
   the marginal likelihood; say so wherever it is compared to the
   Satterthwaite-based existing results.
5. Fallback, pre-declared: if per-protein joint fits fail wholesale
   (> 20% non-convergence), downgrade to two-step (EB factor scores as
   outcome, lmerTest per protein) and label every resulting p as
   two-step, which understates uncertainty in the scores.

Timing gate: time one joint fit before looping. If it exceeds 5 s, cut
the subset to the intersection of the 931 with the 37 F06 lead proteins
plus a random 100-protein calibration sample (seed 42), and say exactly
that in FINDINGS.md.

## Q3 — joint intensity + detection (stretch)

Attempted only after Q1 and Q2 are finished and written up. Scope: the
untested-protein list from
`03_Features/01_Proteins/a_script/03_untested_proteins.R` (the ~34
proteins the missingness filter admits but the cell-means model cannot
test). Per protein: log-intensity (Gaussian) and detected/not-detected
(binomial) share one latent variable in galamm, `(1 | subject)`,
blood index on the Gaussian response. Deliverable is feasibility —
convergence count and interval widths — not discovery. Any protein
promoted from this model is subject to rules 3-5 below unchanged.

## Confirmation rules for any BH survivor, all questions

1. Arm-label permutation across subjects, B = 200, reusing the
   subject-as-unit scheme of
   `03_Features/01_Proteins/a_script/pi_permutation.R::permute_arm_labels`
   (the label is a subject property; timepoint travels with the sample).
   For Q2, where the model carries no arm label, the analog is fixed now:
   permute the six-indicator phenotype block across subjects, leaving the
   protein data untouched. The observed statistic must beat its own
   permuted null at empirical p < 0.05, computed as (n_ge + 1)/(B + 1).
2. A hit the existing lmerTest fit ranked outside its top 200 (by minimum
   p across the six contrasts; for Q2 additionally by |Spearman rho| with
   d_mcsa) is a specification artefact until the coefficient that moved
   is named, with its size in both fits.
3. Every non-converged, singular, or boundary fit is counted and the
   denominator reported next to any count of successes.

## Deliverables and hygiene

Everything new lives under `03_Features/04_galamm_pilot/` (`a_script/`,
`b_reports/`, `c_data/`). No existing script, figure, or result is
modified. Scripts follow the project header convention and `.lintr`;
any function computing a statistic gets a test under `tests/testthat/`.
FINDINGS.md answers Q1, Q2, Q3 in order, each explicit negative stated
as a negative. The README paragraph is written but not inserted.
