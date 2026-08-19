# galamm pilot — findings

Run 2026-08-19 under the rules fixed in PREREG.md. One figure summarises
everything: `b_reports/F_galamm_pilot.{png,pdf}`. galamm 0.4.0, seed 42,
931 complete-case proteins, blood index in every model.

## Q1 — per-timepoint residual variance: negative

The homoscedastic fit was fine. Across 929 converged proteins of 931
(2 failed at the homoscedastic stage; both reported in
`c_data/01_q1_variance.csv`), the median residual-SD ratios are
sigma_T2/sigma_T1 = 0.87 and sigma_T3/sigma_T1 = 0.86 — both inside the
pre-declared [0.8, 1.25] window, and below 1 rather than above it: with
the blood index in the mean model, T3 shows no residual variance
inflation left to absorb. 167 proteins reject the equal-variance LRT at
nominal p < 0.05 (18%, above the 5% chance rate), so heteroscedasticity
exists protein-by-protein, but it changes no conclusion: zero BH q < 0.05
proteins in all six contrasts under both error models, and the
per-contrast p-values sit on the diagonal (figure panel B). The single-
sigma spec in `shared_hlm.R` was not holding signal down.

## Q2 — latent hypertrophic-response factor: estimable, and negative on proteins

The measurement model is identifiable at n = 16, which was not a given.
With `comp_hypertrophy` anchored at 1, the fibre-CSA indicators carry the
factor (d_fcsa_I 1.03, SE 0.34; d_fcsa_II 1.02, SE 0.34; both |z| = 3.0;
~55% of variance explained each). `d_1rm_legpress` loads at 0.57
(SE 0.31, z = 1.8, 28%), `d_mcsa` at 0.27 (SE 0.28, z = 0.9, 8%), and
`d_1rm_ext` at 0.00. The per-item residual variance
(`dispformula = ~ (1 | item)`) failed at its starting values, so the
common-sigma model is reported; the attempt is logged by the script.

The pre-declared d_mcsa question resolves the F06 reading, in the
direction nobody scripted: the factor is **not** d_mcsa (Spearman
rho = 0.32 between EB scores and d_mcsa) — it is the composite/fCSA axis
(rho = 0.94 vs `comp_hypertrophy`). So the pilot neither explains the F06
d_mcsa concentration (the factor ranks subjects differently) nor extends
it. Whole-muscle CSA and fibre-level hypertrophy are two axes in this
cohort, and the six-indicator trait is the fibre one.

On proteins the answer is a clean negative. All 931 joint fits converged
(0.05 s median, timing gate untouched); the smallest BH q for the
protein-on-factor loading is 0.117; **zero survivors**, so the
permutation stage never armed. The top nominal proteins (NIT2, CAND2,
GSTP1, ECHS1, LANCL1) are nominal only and are shown as such in panel E.
The dichotomy is not discarding protein-level signal: modelling the
response as a continuous latent trait, within arm, finds nothing the
median split missed.

## Q3 — joint intensity + detection: not attempted

Not attempted in this run. The prereg scopes it as a stretch after Q1
and Q2; both are now closed, so the ~34 untestable proteins remain an
open, pre-registered follow-up. No claim is made about them.

## Deviations from PREREG.md

- The per-item dispersion attempt in the measurement model fell back to
  a common residual variance (declared fallback, recorded above).
- Nothing else: subset, covariate, thresholds, gates and stop rules ran
  as written.

## Proposed README paragraph (not inserted)

> A galamm pilot (`03_Features/04_galamm_pilot/`, pre-registered in its
> PREREG.md) tested whether three specification choices were holding
> signal down: a single residual variance across timepoints, the HR/LR
> median split, and Satterthwaite-vs-marginal inference. None was. Per-
> timepoint residual variances sit at median ratios 0.87 (T2/T1) and
> 0.86 (T3/T1) with zero BH survivors under either error model; a latent
> hypertrophic-response factor is estimable from the six phenotype
> indicators but is the composite/fCSA axis, not d_mcsa (Spearman
> rho = 0.32), and associates with none of the 931 complete-case
> proteins at BH q < 0.05 (minimum q = 0.117, all fits converged). The
> null survives.
