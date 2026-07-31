# Method ranking

Written 2026-07-31, against commit `9d6d1e50`, **before any alternative method was
fitted to this dataset**. That ordering is the point. A ranking written after the
alternatives have been run is a ranking fitted to the answer, which is the failure
the `d_mcsa` criterion was fixed in advance to avoid (`docs/decisions.md`,
2026-07-30).

Each section ranks candidates on design arguments and published benchmark
evidence. Whichever candidate ranks first becomes the primary analysis, whatever
it yields. Where the incumbent method wins, it wins on merit and not by being
already installed.

## The decision rule

1. Rank 1 is primary. Ranks 2 and 3 are pre-declared sensitivity arms, reported
   whether they agree or not.
2. A method ranked below the incumbent is not run.
3. If rank 1 is not the incumbent, the incumbent is refitted alongside it and both
   are reported, with the ranking's date and commit quoted.
4. Nothing in this file is revised after a result is seen. Revisions get a new
   dated section stating what new evidence forced them.

## 1. Differential abundance estimator

1. **limma + `duplicateCorrelation`** (incumbent, retained). Empirical Bayes
   variance moderation is the best-evidenced small-n intervention available
   (Phipson 2016, `10.1214/16-AOAS920`; Kammers 2015,
   `10.1016/j.euprot.2015.02.002`), and the independent `prolfqua` benchmark
   recommends fixed-effect linear models with variance moderation over mixed
   models when sample sizes are small (`10.1021/acs.jproteome.2c00441`).
2. **msqrob2**, as the primary sensitivity arm. The only candidate offering both a
   per-protein random effect and a missingness component in one fit
   (`10.1074/mcp.RA119.001624`, `10.1021/acs.analchem.9b04375`).
3. **proDA**, as the MNAR-specific arm (`10.1101/661496`, preprint only).
4. DEqMS, a better variance prior but no change to blocking or missingness.
5. MSstats, not run. The only method in the `prolfqua` benchmark flagged for
   anti-conservative FDR at small nominal thresholds, which is the one property a
   null result cannot afford.

**Why the incumbent survives a real challenge.** Hoffman & Roussos
(`10.1093/bioinformatics/btaa687`) is the strongest published case against a
single shared correlation: methods without a per-feature random effect produce
reproducible false positives driven by between-individual variance. That paper's
setting is large cohorts with many samples per subject. Ours is 16 subjects with
two or three observations each, where a per-protein variance component is barely
identified. In the `prolfqua` benchmark the mixed model estimated the fewest
contrasts of any method (94.3%) and lost ranking accuracy to unstable denominator
degrees of freedom. Convergence failure is not random across proteins, so it would
bias any comparison built on it.

**The finding that does not favour the incumbent.** Neither of the two arms we
currently run models missingness. Complete-case contrast estimates are biased
under non-ignorable missingness (O'Brien 2018, `10.1214/18-aoas1144`), and
missForest is an ad-hoc imputer whose accuracy degrades as the MNAR fraction rises
(Jin 2021, `10.1038/s41598-021-81279-4`) and which was validated under essentially
MCAR conditions (`10.1093/bioinformatics/btr597`). Running msqrob2 and proDA is
therefore not optional polish; it tests the one assumption the primary arm cannot.

## 2. Multiplicity

1. **BH within each contrast, contrasts pre-declared, with Benjamini–Bogomolov
   adjustment when contrasts are reported selectively** (`10.1111/rssb.12028`).
2. **BH within each contrast as-is** (incumbent), with an explicit statement that
   the guarantee is per-contrast and that expected false discoveries across all
   nine run up to nine times the per-contrast budget.
3. A contrast-level graphical scheme gating BH inside each opened contrast
   (`10.1002/sim.3495`).
4. IHW within each contrast (`10.1038/nmeth.3885`). Not run: 1900 hypotheses is
   thin for data-driven weight learning, and the most tempting covariate,
   completeness, is entangled with MNAR and so not null-independent.
5. Storey q within each contrast. Not run: it gains power only when π₀ is well
   below 1, which is the opposite of this dataset.
6. **BH pooled across all nine. Not run, and not merely because it is
   conservative.** Efron (`10.1214/07-AOAS141`) shows that pooling families with
   different null distributions and signal densities distorts the empirical null
   and buries sparse-family signal. Within-arm, between-arm and interaction
   contrasts are demonstrably not exchangeable here. Pooling is the wrong choice,
   not the safe one.

The incumbent is exactly right on the separate-versus-pooled question and
incomplete on selective reporting. The gap is one adjustment and one sentence.

## 3. Set-level testing

1. **`limma::fry`** over the detected-protein background, same design matrix and
   subject block as the DE fit. The only option testing a null defensible at n=16
   that also respects the blocking, and rotation needs no permutation budget the
   unbalanced cells cannot supply (`10.1093/bioinformatics/btq401`).
2. **`camera`** alongside it as the competitive complement, the only competitive
   test with a published inter-gene-correlation correction (`10.1093/nar/gks461`).
3. `mroast` if the directional decomposition is wanted.
4. **fgsea, demoted from inference to ranking and display.** Its null permutes
   gene labels and assumes proteins vary independently. Goeman & Bühlmann
   (`10.1093/bioinformatics/btm051`) call such p-values wildly anti-conservative;
   Gatti (`10.1186/1471-2164-11-574`) measured very high false positive rates on
   randomised labels, with correlated sets most likely to be called; Tamayo
   (`10.1177/0962280212460441`) gives the variance-inflation mechanism.

**The set-size floor gets a citation and a correction.** Reimand
(`10.1038/s41596-018-0103-9`) sets 10–15 minimum and 200–500 maximum members, and
GSEA's own default minimum is 15. The floor must be applied to members **detected
in our 1900**, not annotated in the ontology (Wijesooriya, `10.1371/journal.pcbi.1009935`;
Timmons, `10.1186/s13059-015-0761-7`). A Hallmark set with 9 of 200 members
measured is a 9-gene set. This retires the current practice of filtering at
display time: the floor belongs at set construction, which is where the pipeline
does not yet apply it.

## 4. Prediction and classification screens

1. **Nested leave-one-subject-*pair*-out for AUC, with a max-statistic
   subject-level permutation null over the whole screen.**
2. Repeated stratified 5-fold grouped by subject, full distribution reported.
3. Nested LOSO retained but with AUC replaced by a metric well defined under
   pooling.
4. `.632+` bootstrap over subjects.
5. **Nested LOSO with pooled AUC, the incumbent. Do not keep.**

**The 144 cells at AUC exactly 0 are a named artifact, not a curiosity.** AUC is a
pairwise statistic and a leave-one-out fold holds one observation, so the only way
to obtain a number is to pool out-of-fold predictions, which Airola
(`10.1016/j.csda.2010.11.018`, `10.1007/s10618-018-00607-x`) shows carries a
substantial negative bias. Parker (`10.1186/1471-2105-8-326`) names the mechanism:
with balanced classes, removing one observation tilts the training fold away from
the held-out class, 8/8 becomes 8/7, so a degenerate intercept-only model predicts
the training majority and is systematically wrong. Pooled AUC then goes to exactly
0 rather than 0.5. **An AUC of 0 here means no signal plus a degenerate model. It
has never meant inverted signal, and the figure must not imply otherwise.**

**Our permutation null is not a valid null for a selected maximum.** Nichols &
Holmes (`10.1002/hbm.1058`) require permuting labels, rerunning the entire screen
including model selection, taking the maximum within each permutation, and
comparing the observed maximum to that distribution. We built per-cell nulls and
then looked at the maximum across cells. That tests each cell, not the argmax, so
every screen-maximum claim in F06 and F07 currently lacks its null. The exchangeable
unit is the subject, which the existing permutation code already respects.

**Riley** (`10.1002/sim.7992`) is to be quoted in the write-up: at 16 subjects and
8 events no minimum-sample-size criterion for a prediction model is met by any
margin. That is the reason this is framed as a screen and never as a model.

## 5. Imputation and leakage

1. **In-fold imputation** for anything cross-validated, refitted inside every
   training fold (Kaufman, `10.1145/2382577.2382579`).
2. **Complete-case restriction** as a pre-registered outcome-blind co-primary for
   the prediction screen. Bourgon (`10.1073/pnas.0914005107`) makes a marginal,
   group-blind completeness filter valid and power-increasing.
3. Modelled missingness (proDA, msqrob2 hurdle) for the differential-abundance
   arms, replacing imputation rather than tuning it.
4. Cohort-wide missForest as a clearly labelled non-cross-validated sensitivity
   arm only.
5. **Cohort-wide missForest before cross-validation. Do not use.**

**The defence we have been leaning on is false.** The argument that response-blind
preprocessing is safe before cross-validation is exactly what Moscovich & Rosset
(`10.1111/rssb.12537`) disprove: unsupervised preprocessing introduces substantial
bias into cross-validation estimates, **of either sign**, with small n and high
dimension named as the worst case. So "any leakage would only have made our null
more conservative" does not hold, and cannot be used. Rosenblatt
(`10.1038/s41467-024-46150-w`) adds that small datasets exacerbate leakage and that
it shifts coefficients, not just scores.

Two aggravating features specific to this design, argued rather than cited:
missForest learns cross-sample structure, so a cohort-wide fit uses the held-out
subject to impute others and others to impute the held-out subject; and with three
timepoints per subject, a subject's own other timepoints are its nearest
neighbours in feature space. Cohort-wide imputation quietly reconstructs a
held-out sample from its own subject, which is subject-level leakage wearing an
imputation costume.

**On the complete-case rescore already run.** It qualifies as rank 2 and its
validity holds only because completeness was computed marginally across all 45
samples, blind to arm and timepoint. A per-cell completeness filter would be
group-dependent and would forfeit the Bourgon guarantee. The cost is coverage, not
Type I error: 931 of 1900 proteins is enriched for abundant proteins and depleted
of the low-abundance range where detection-limit truncation lives, so "nothing in
931" is a weaker statement than "nothing in 1900" and must be reported as such.

## What this ranking changes

| Area | Incumbent | Verdict |
|---|---|---|
| DA estimator | limma + dupCor | **Retained as rank 1.** Add msqrob2 and proDA as pre-declared arms, since neither current arm models missingness. |
| Multiplicity | BH within contrast | **Retained as rank 2, promoted to rank 1 by adding Benjamini–Bogomolov.** Pooling is ruled out on Efron, not on conservatism. |
| Set test | fry beside fgsea | **fry confirmed rank 1; fgsea demoted to display.** Move the detected-member floor from display time to set construction. |
| Screen metric | pooled AUC under LOSO | **Ranked last.** Replace with leave-one-subject-pair-out, or with a metric defined under pooling. |
| Screen null | per-cell permutation | **Invalid for a selected maximum.** Needs a max-statistic null that reruns selection inside each permutation. |
| Imputation | cohort-wide missForest | **Ranked last for anything cross-validated.** In-fold or complete-case only. |

Four of the six move. The two that do not, within-contrast BH and fry, were the
two most likely to be dismissed as over-conservative, and the literature says both
are correct for reasons stronger than the ones the pipeline currently gives.
